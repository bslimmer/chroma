import csv
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]


def child_update_no(csv_path):
    update_suffix = csv_path.name.rsplit('_cfg_', 1)[1]
    return int(update_suffix.split('.gluecor', 1)[0])


def load_series(data_folder, glob_pattern):
    csv_files = sorted(
        Path(data_folder).glob(glob_pattern),
        key=lambda csv_path: (child_update_no(csv_path), csv_path.name),
    )

    if not csv_files:
        raise ValueError("No matching gluecor CSV files found in " + str(data_folder))

    return np.stack(
        [np.loadtxt(csv_file, delimiter=',', skiprows=1, usecols=1) for csv_file in csv_files]
    ), csv_files


def load_child_series(data_folder, child_label):
    child_series, csv_files = load_series(data_folder, child_label + "_hmc.*.gluecor.csv")

    if not len(csv_files):
        raise ValueError(
            "No " + child_label + "_hmc gluecor CSV files found in " + str(data_folder)
        )

    return child_series


def parse_int_vector(text):
    return [int(value) for value in text.split()]


def load_sidecar_maps(sidecar_path):
    root = ET.parse(sidecar_path).getroot()
    parent_nrow = parse_int_vector(root.findtext('./Parent/nrow'))
    parent_temporal_extent = parent_nrow[-1]

    child_maps = {}
    for child_tag in ('Child0', 'Child1'):
        local_to_global = parse_int_vector(root.findtext(f'./{child_tag}/local_to_global_t'))
        child_maps[child_tag.lower()] = local_to_global

    return parent_temporal_extent, child_maps


def periodic_min_separation(t0, t1, temporal_extent):
    forward_delta = (t1 - t0) % temporal_extent
    reverse_delta = (t0 - t1) % temporal_extent
    return min(forward_delta, reverse_delta)


def plot_annotated_matrix(
    matrix,
    title,
    output_path,
    use_heatmap=False,
    x_label='child1 timeslice index',
    y_label='child0 timeslice index',
    x_tick_labels=None,
    y_tick_labels=None,
):
    fig, ax = plt.subplots(figsize=(8, 8))

    if use_heatmap:
        ax.imshow(matrix, cmap='bwr', alpha=0.2, origin='lower')
    else:
        ax.imshow(matrix, cmap='bwr', alpha=0.0, origin='lower')

    nrows, ncols = matrix.shape
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xticks(np.arange(ncols))
    ax.set_yticks(np.arange(nrows))
    if x_tick_labels is None:
        x_tick_labels = [str(col) for col in range(ncols)]
    if y_tick_labels is None:
        y_tick_labels = [str(row) for row in range(nrows)]
    ax.set_xticklabels(x_tick_labels)
    ax.set_yticklabels(y_tick_labels)

    if any(len(label) > 1 for label in x_tick_labels):
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    for row in range(nrows):
        for col in range(ncols):
            ax.text(
                col,
                row,
                f"{matrix[row, col]:.3e}",
                ha='center',
                va='center',
                fontsize=9,
                color='black',
            )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def jackknife_std_from_estimates(jackknife_estimates):
    nsamples = jackknife_estimates.shape[0]
    if nsamples < 2:
        return np.zeros_like(jackknife_estimates[0])

    jackknife_mean = np.mean(jackknife_estimates, axis=0)
    return np.sqrt(
        #(nsamples - 1) / nsamples
        (nsamples - 1)
        * np.sum((jackknife_estimates - jackknife_mean) ** 2, axis=0)
    )


def jackknife_child_pair_product(child_pair_outer_products, child0_avg_vectors, child1_avg_vectors):
    nsamples = child_pair_outer_products.shape[0]

    raw_child_pair_product_avg = np.mean(child_pair_outer_products, axis=0)
    child0_onepoint_avg = np.mean(child0_avg_vectors, axis=0)
    child1_onepoint_avg = np.mean(child1_avg_vectors, axis=0)

    # Connected two-point function:
    # <O0(t0) O1(t1)> - <O0(t0)> <O1(t1)>
    child_pair_product_avg = (
        raw_child_pair_product_avg
        - np.outer(child0_onepoint_avg, child1_onepoint_avg)
    )


    if nsamples < 2:
        return (
            child_pair_product_avg,
            np.zeros_like(child_pair_product_avg),
            np.stack([child_pair_product_avg]),
        )

    raw_sum = np.sum(child_pair_outer_products, axis=0)
    child0_sum = np.sum(child0_avg_vectors, axis=0)
    child1_sum = np.sum(child1_avg_vectors, axis=0)

    jackknife_estimates = []
    for sample_index in range(nsamples):
        loo_raw_avg = (raw_sum - child_pair_outer_products[sample_index]) / (nsamples - 1)
        loo_child0_avg = (child0_sum - child0_avg_vectors[sample_index]) / (nsamples - 1)
        loo_child1_avg = (child1_sum - child1_avg_vectors[sample_index]) / (nsamples - 1)

        loo_connected_avg = loo_raw_avg - np.outer(loo_child0_avg, loo_child1_avg)
        jackknife_estimates.append(loo_connected_avg)

    jackknife_estimates = np.stack(jackknife_estimates)
    jackknife_std = jackknife_std_from_estimates(jackknife_estimates)

    return child_pair_product_avg, jackknife_std, jackknife_estimates


def jackknife_mean_and_std_from_measurements(measurements):
    measurements = np.asarray(measurements, dtype=float)
    mean = np.mean(measurements)

    if measurements.size < 2:
        return mean, 0.0

    total = np.sum(measurements)
    jackknife_estimates = np.array(
        [
            (total - measurements[measurement_index]) / (measurements.size - 1)
            for measurement_index in range(measurements.size)
        ],
        dtype=float,
    )
    stddev = jackknife_std_from_estimates(jackknife_estimates)
    return mean, stddev


def group_measurement_cube_by_periodic_separation(
    measurement_cube,
    child0_global_t,
    child1_global_t,
    temporal_extent,
):
    values_by_separation = {}

    for sample_index in range(measurement_cube.shape[0]):
        for child0_local_t, child0_parent_t in enumerate(child0_global_t):
            for child1_local_t, child1_parent_t in enumerate(child1_global_t):
                separation = periodic_min_separation(
                    child0_parent_t,
                    child1_parent_t,
                    temporal_extent,
                )
                values_by_separation.setdefault(separation, []).append(
                    measurement_cube[sample_index, child0_local_t, child1_local_t]
                )

    separations = np.array(sorted(values_by_separation), dtype=int)
    grouped_measurements = [
        np.array(values_by_separation[separation], dtype=float) for separation in separations
    ]
    counts = np.array(
        [grouped_values.size for grouped_values in grouped_measurements],
        dtype=int,
    )
    return separations, grouped_measurements, counts


def reduce_periodic_separation_bins(
    measurement_cube,
    child0_global_t,
    child1_global_t,
    temporal_extent,
):
    separations, grouped_measurements, counts = group_measurement_cube_by_periodic_separation(
        measurement_cube,
        child0_global_t,
        child1_global_t,
        temporal_extent,
    )

    averages = []
    stddevs = []
    for measurements in grouped_measurements:
        mean, stddev = jackknife_mean_and_std_from_measurements(measurements)
        averages.append(mean)
        stddevs.append(stddev)

    return (
        separations,
        np.array(averages, dtype=float),
        np.array(stddevs, dtype=float),
        counts,
    )


def write_pair_map_csv(
    output_path,
    connected_avg_matrix,
    connected_std_matrix,
    child0_global_t,
    child1_global_t,
    temporal_extent,
    row_name='child0',
    col_name='child1',
):
    with output_path.open('w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            [
                row_name + '_local_t',
                col_name + '_local_t',
                row_name + '_global_t',
                col_name + '_global_t',
                'forward_delta_t_mod_T',
                'periodic_min_delta_t',
                'connected_avg',
                'connected_jackknife_std',
            ]
        )

        for child0_local_t, child0_parent_t in enumerate(child0_global_t):
            for child1_local_t, child1_parent_t in enumerate(child1_global_t):
                forward_delta = (child1_parent_t - child0_parent_t) % temporal_extent
                periodic_delta = periodic_min_separation(
                    child0_parent_t,
                    child1_parent_t,
                    temporal_extent,
                )
                writer.writerow(
                    [
                        child0_local_t,
                        child1_local_t,
                        child0_parent_t,
                        child1_parent_t,
                        forward_delta,
                        periodic_delta,
                        connected_avg_matrix[child0_local_t, child1_local_t],
                        connected_std_matrix[child0_local_t, child1_local_t],
                    ]
                )


def write_separation_csv(output_path, separations, averages, stddevs, counts):
    with output_path.open('w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            [
                'periodic_delta_t',
                'pair_count',
                'connected_avg',
                'connected_jackknife_std',
            ]
        )

        for separation, count, average, stddev in zip(separations, counts, averages, stddevs):
            writer.writerow([separation, count, average, stddev])


def plot_separation_series(separations, averages, stddevs, output_path, title):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(separations[1:], averages[1:], yerr=stddevs[1:], fmt='o', linestyle='none', capsize=4)
    ax.set_title(title)
    ax.set_xlabel('periodic time separation on combined lattice')
    ax.set_ylabel('connected two-point function')
    ax.set_xticks(separations)
    ax.grid(True, alpha=0.3)
    #ax.set_yscale('log')
    #ax.set_ylim(1e-09, 1e-02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_combined_separation_series(series_specs, output_path, title):
    fig, ax = plt.subplots(figsize=(8, 5))

    for series_spec in series_specs:
        separations = series_spec['separations']
        averages = series_spec['averages']
        stddevs = series_spec['stddevs']

        if len(separations) == 0:
            continue

        ax.errorbar(
            separations[1:],
            averages[1:],
            yerr=stddevs[1:],
            fmt=series_spec.get('fmt', 'o'),
            linestyle=series_spec.get('linestyle', 'none'),
            capsize=4,
            label=series_spec['label'],
        )

    ax.set_title(title)
    ax.set_xlabel('periodic time separation on combined lattice')
    ax.set_ylabel('connected two-point function')

    all_separations = [
        separation
        for series_spec in series_specs
        for separation in series_spec['separations']
    ]
    if all_separations:
        ax.set_xticks(sorted(set(all_separations)))

    ax.grid(True, alpha=0.3)
    #ax.set_yscale('log')
    #ax.set_ylim(1e-09, 1e-02)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def build_outer_product_cube(source_vectors, sink_vectors):
    return np.stack(
        [
            np.outer(source_vectors[sample_index], sink_vectors[sample_index])
            for sample_index in range(source_vectors.shape[0])
        ]
    )



def two_level_data_analysis(vol):
    child_pair_outer_products = []
    child_pair_labels = []
    child0_avg_vectors = []
    child1_avg_vectors = []
    parent_temporal_extent = None
    child0_global_t = None
    child1_global_t = None

    for parent_stream in range(10):
        for parent_substream in range(10, 101, 10):
            if parent_stream == 0:
                file_prefix = 'hmc_parent_cfg_' + str(parent_substream)
            else:
                file_prefix = 'hmc_parent' + str(parent_stream) + '_cfg_' + str(parent_substream)

            child_run_dir = SCRIPT_DIR / 'child_hmc_runs' / file_prefix
            completion_file = child_run_dir / 'gluecor_output' / ('child1_hmc.' + file_prefix + '_cfg_1000.gluecor.csv')

            # Check if the full child HMC stream is done; skip if not.
            if completion_file.exists():
                data_folder = child_run_dir / 'gluecor_output'
                sidecar_path = child_run_dir / (file_prefix + '.sidecar.xml')

                current_parent_temporal_extent, current_child_maps = load_sidecar_maps(sidecar_path)
                current_child0_global_t = current_child_maps['child0']
                current_child1_global_t = current_child_maps['child1']

                if parent_temporal_extent is None:
                    parent_temporal_extent = current_parent_temporal_extent
                    child0_global_t = current_child0_global_t
                    child1_global_t = current_child1_global_t
                else:
                    if parent_temporal_extent != current_parent_temporal_extent:
                        raise ValueError(
                            'Inconsistent parent temporal extent between child-run sidecars'
                        )
                    if child0_global_t != current_child0_global_t:
                        raise ValueError('Inconsistent child0 local-to-global time map')
                    if child1_global_t != current_child1_global_t:
                        raise ValueError('Inconsistent child1 local-to-global time map')

                # Each array has shape (n_measurements, n_timeslices).
                # Convert raw timeslice sums into spatial averages.
                child0_data = load_child_series(data_folder, 'child0') / vol
                child1_data = load_child_series(data_folder, 'child1') / vol

                child0_avg = np.mean(child0_data, axis=0)
                child1_avg = np.mean(child1_data, axis=0)

                # Store the raw per-parent connected-building block first.
                child_pair_outer_products.append(np.outer(child0_avg, child1_avg))
                child_pair_labels.append(file_prefix)
                child0_avg_vectors.append(child0_avg)
                child1_avg_vectors.append(child1_avg)

    if child_pair_outer_products:
        # Shape: (n_parent_child_sets, child0_timeslice, child1_timeslice)
        child_pair_outer_products = np.stack(child_pair_outer_products)
        child0_avg_vectors = np.stack(child0_avg_vectors)
        child1_avg_vectors = np.stack(child1_avg_vectors)

        child_pair_product_avg, child_pair_product_std, child_pair_product_jackknife = (
            jackknife_child_pair_product(
                child_pair_outer_products,
                child0_avg_vectors,
                child1_avg_vectors,
            )
        )

        periodic_separations, periodic_averages, periodic_stddevs, periodic_pair_counts = (
            reduce_periodic_separation_bins(
                child_pair_outer_products,
                child0_global_t,
                child1_global_t,
                parent_temporal_extent,
            )
        )

        # Convert to log for the stddev matrix plot.
        child_pair_product_std_plot = np.full_like(child_pair_product_std, np.nan, dtype=float)
        positive_std_mask = child_pair_product_std > 0.0
        child_pair_product_std_plot[positive_std_mask] = np.log10(
            child_pair_product_std[positive_std_mask]
        )

    else:
        child_pair_outer_products = np.empty((0, 0, 0))
        child_pair_product_avg = np.empty((0, 0))
        child_pair_product_std = np.empty((0, 0))
        child_pair_product_std_plot = np.empty((0, 0))
        periodic_separations = np.empty((0,), dtype=int)
        periodic_averages = np.empty((0,), dtype=float)
        periodic_stddevs = np.empty((0,), dtype=float)
        periodic_pair_counts = np.empty((0,), dtype=int)

    output_dir = REPO_ROOT / 'cfgs/two_level_glueball_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    print('loaded child-pair sets:', len(child_pair_labels))
    if child_pair_product_avg.size > 0:
        plot_annotated_matrix(
            child_pair_product_avg,
            'Child Pair Product Average',
            output_dir / 'child_pair_product_avg.png',
            use_heatmap=True
        )
        plot_annotated_matrix(
            child_pair_product_std_plot,
            'Child Pair Product Standard Deviation log10',
            output_dir / 'child_pair_product_std_log10.png',
            use_heatmap=True,
        )
        plot_annotated_matrix(
            child_pair_product_avg,
            'Child Pair Product Average (Parent Global Timeslices)',
            output_dir / 'child_pair_product_avg_global_time.png',
            use_heatmap=True,
            x_label='child1 parent-global timeslice',
            y_label='child0 parent-global timeslice',
            x_tick_labels=[str(t) for t in child1_global_t],
            y_tick_labels=[str(t) for t in child0_global_t],
        )
        plot_annotated_matrix(
            child_pair_product_std_plot,
            'Child Pair Product Standard Deviation log10 (Parent Global Timeslices)',
            output_dir / 'child_pair_product_std_log10_global_time.png',
            use_heatmap=True,
            x_label='child1 parent-global timeslice',
            y_label='child0 parent-global timeslice',
            x_tick_labels=[str(t) for t in child1_global_t],
            y_tick_labels=[str(t) for t in child0_global_t],
        )
        write_pair_map_csv(
            output_dir / 'child_pair_product_pair_map.csv',
            child_pair_product_avg,
            child_pair_product_std,
            child0_global_t,
            child1_global_t,
            parent_temporal_extent,
        )
        write_separation_csv(
            output_dir / 'child_pair_product_by_periodic_dt.csv',
            periodic_separations,
            periodic_averages,
            periodic_stddevs,
            periodic_pair_counts,
        )
        plot_separation_series(
            periodic_separations,
            periodic_averages,
            periodic_stddevs,
            output_dir / 'child_pair_product_by_periodic_dt.png',
            'Child Pair Product Aggregated By Periodic Time Separation',
        )

    return periodic_separations, periodic_averages, periodic_stddevs


def parent_lattice_data_analysis(vol):
    data_folder = SCRIPT_DIR / 'parent_files' / 'gluecor_output'

    try:
        parent_data, csv_files = load_series(data_folder, 'hmc_parent*.gluecor.csv')
    except ValueError:
        parent_data = np.empty((0, 0))
        csv_files = []

    if parent_data.size > 0:
        parent_data = parent_data / vol
        parent_labels = [csv_file.stem.replace('.gluecor', '') for csv_file in csv_files]

        parent_pair_outer_products = build_outer_product_cube(parent_data, parent_data)
        parent_pair_product_avg, parent_pair_product_std, parent_pair_product_jackknife = (
            jackknife_child_pair_product(
                parent_pair_outer_products,
                parent_data,
                parent_data,
            )
        )

        parent_temporal_extent = parent_data.shape[1]
        parent_global_t = list(range(parent_temporal_extent))
        periodic_separations, periodic_averages, periodic_stddevs, periodic_pair_counts = (
            reduce_periodic_separation_bins(
                parent_pair_outer_products,
                parent_global_t,
                parent_global_t,
                parent_temporal_extent,
            )
        )

        parent_pair_product_std_plot = np.full_like(parent_pair_product_std, np.nan, dtype=float)
        positive_std_mask = parent_pair_product_std > 0.0
        parent_pair_product_std_plot[positive_std_mask] = np.log10(
            parent_pair_product_std[positive_std_mask]
        )
    else:
        parent_labels = []
        parent_pair_outer_products = np.empty((0, 0, 0))
        parent_pair_product_avg = np.empty((0, 0))
        parent_pair_product_std = np.empty((0, 0))
        parent_pair_product_jackknife = np.empty((0, 0, 0))
        parent_pair_product_std_plot = np.empty((0, 0))
        parent_temporal_extent = 0
        parent_global_t = []
        periodic_separations = np.empty((0,), dtype=int)
        periodic_averages = np.empty((0,), dtype=float)
        periodic_stddevs = np.empty((0,), dtype=float)
        periodic_pair_counts = np.empty((0,), dtype=int)

    output_dir = REPO_ROOT / 'cfgs/two_level_glueball_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    print('loaded parent measurement sets:', len(parent_labels))
    if parent_pair_product_avg.size > 0:
        plot_annotated_matrix(
            parent_pair_product_avg,
            'Parent Pair Product Average',
            output_dir / 'parent_pair_product_avg.png',
            use_heatmap=True
        )
        plot_annotated_matrix(
            parent_pair_product_std_plot,
            'Parent Pair Product Standard Deviation log10',
            output_dir / 'parent_pair_product_std_log10.png',
            use_heatmap=True,
        )
        plot_annotated_matrix(
            parent_pair_product_avg,
            'Parent Pair Product Average (Parent Global Timeslices)',
            output_dir / 'parent_pair_product_avg_global_time.png',
            use_heatmap=True,
            x_label='sink parent-global timeslice',
            y_label='source parent-global timeslice',
            x_tick_labels=[str(t) for t in parent_global_t],
            y_tick_labels=[str(t) for t in parent_global_t],
        )
        plot_annotated_matrix(
            parent_pair_product_std_plot,
            'Parent Pair Product Standard Deviation log10 (Parent Global Timeslices)',
            output_dir / 'parent_pair_product_std_log10_global_time.png',
            use_heatmap=True,
            x_label='sink parent-global timeslice',
            y_label='source parent-global timeslice',
            x_tick_labels=[str(t) for t in parent_global_t],
            y_tick_labels=[str(t) for t in parent_global_t],
        )
        write_pair_map_csv(
            output_dir / 'parent_pair_product_pair_map.csv',
            parent_pair_product_avg,
            parent_pair_product_std,
            parent_global_t,
            parent_global_t,
            parent_temporal_extent,
            row_name='source_parent',
            col_name='sink_parent',
        )
        write_separation_csv(
            output_dir / 'parent_pair_product_by_periodic_dt.csv',
            periodic_separations,
            periodic_averages,
            periodic_stddevs,
            periodic_pair_counts,
        )
        plot_separation_series(
            periodic_separations,
            periodic_averages,
            periodic_stddevs,
            output_dir / 'parent_pair_product_by_periodic_dt.png',
            'Parent Pair Product Aggregated By Periodic Time Separation',
        )

    return periodic_separations, periodic_averages, periodic_stddevs


def main():
    vol = 8.0**3
    two_level_periodic_separations, two_level_periodic_averages, two_level_periodic_stddevs = (
        two_level_data_analysis(vol)
    )
    parent_periodic_separations, parent_periodic_averages, parent_periodic_stddevs = (
        parent_lattice_data_analysis(vol)
    )

    output_dir = REPO_ROOT / 'cfgs/two_level_glueball_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_combined_separation_series(
        [
            {
                'label': 'Two-level child pair',
                'separations': two_level_periodic_separations,
                'averages': two_level_periodic_averages,
                'stddevs': two_level_periodic_stddevs,
                'fmt': 'o',
                'linestyle': 'none',
            },
            {
                'label': 'Parent lattice',
                'separations': parent_periodic_separations,
                'averages': parent_periodic_averages,
                'stddevs': parent_periodic_stddevs,
                'fmt': 's',
                'linestyle': 'none',
            },
        ],
        output_dir / 'combined_pair_product_by_periodic_dt.png',
        'Two-Level And Parent Periodic Time-Separation Comparison',
    )
    


if __name__ == "__main__":
    main()
