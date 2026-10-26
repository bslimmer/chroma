# Edge-SVD Boundary Projectors for the Factorized Perambulator

Status: two-measurement design implemented 2026-10-25 (working tree, not committed).
Compiles; the shared code passes a standalone round-trip test; not yet run inside Chroma
with real solvers. The 2026-10-22 version computed the basis inside the prop measurement;
that computation now lives in its own measurement, which writes the basis to a file, and
the prop measurement reads that file. Reason: the 2 h wall-time cap. A 32^3x64 prop run
with 128 colorvecs already takes about 2 h, and the basis adds about 864 slab solves at
R = 128 (k = 144, 2 power iterations), or about 3200 at R = 512.
Target files: new `lib/meas/inline/hadron/inline_edge_svd_projectors_superb_w.{cc,h}`;
`inline_prop_and_matelem_distillation_superb_w.{cc,h}`; `lib/util/ferm/factorized_distillation.{h,cc}`;
`inline_hadron_aggregate.cc`; `lib/Makefile.am`, `lib/CMakeLists.txt`.
References: `specs/hier/factorized_prop.md`, `specs/hier/factorized_prop_overlapping.md`;
Ce, Giusti, Schaefer, arXiv:1601.04587; Halko, Martinsson, Tropp, arXiv:0909.4061;
toys `fact_prop_testing/toy_overlap_vs_nonoverlap/{basis_cmp,check_rsvd,tune}.py`;
walkthrough doc "Edge-SVD Projectors: Code Walkthrough".

## Goal

Replace the vectors that define the boundary projector `Pi` in the factorized prop
measurement. Today `Pi` is a per-spin projector onto Laplacian eigenvectors
(`distillation_vecs`). The new option (`singular_vectors`) uses the dominant right
singular vectors of the edge-to-edge block M that carries the quark across the frozen
region. They are spin-colour vectors that depend only on frozen-region operators.

Two measurements:

1. `EDGE_SVD_PROJECTORS_SUPERB` (new): computes the basis for one configuration with a
   randomized SVD and writes it, with its singular values and metadata, to a file.
2. `PROP_AND_MATELEM_DISTILLATION_SUPERB` (existing): with
   `projector_type = singular_vectors`, reads the basis from that file and uses its first
   R vectors. The chain, cut fields, production solvers, signs, slices, phasings and output
   format are unchanged.

The basis is ordered by decreasing singular value. A file computed at a large rank
R_max therefore serves every R <= R_max, so a rank scan needs one basis run per config.
Taking the first R columns is the best rank-R choice within the computed subspace, and at
least as accurate as a run at rank R. The basis depends only on frozen links, so in a
two-level scheme it is computed once per frozen configuration.

Toy motivation (4^3x16, clover, rough field, per-quark ratio vs. the same chain with
`Pi = 1`, rank R on a 768-dim slice): at R = 32, Laplacian 0.37 / 0.41, edge-SVD
0.71 / 0.77 (non-overlap / overlap).

## Blocks (unchanged from the 2026-10-21 spec)

Notation: `R_t` = restriction to slice(s) t, `D` = full Dirac operator (`fLinOp`),
`D^dag` = `fLinOp` with `MINUS`, `g5` = `Gamma(Ns*Ns-1)`, `(X^-1)^dag = g5 X^-1 g5`.
The basis is the top-R right singular vectors on the projector slice p of `M : p -> q`.
Two frozen intervals `F0 = [0, fw-1]`, `F1 = [c1, c1+fw-1]`; both intervals are handled
in the same solves (M is block-diagonal over them).

Non-overlapping (`M = T^dag`, slab solves with the `u_sep` solver, cuts rooted on
`t_start` and `t_end - 1` of each interval):

```text
M   y = R_q D^dag R_{e_s} g5 S_00 g5 R_{e_p} D^dag y      (y on p)
M^dag x = R_p D R_{e_p} S_00 R_{e_s} D x                  (x on q)
```

| interval | p (Pi) | e_p | e_s | q |
|---|---|---|---|---|
| F1 | c1+fw-1 | c1+fw-2 | c1+1 | c1 |
| F0 | 0 | 1 | fw-2 | fw-1 |

Overlapping (`M = G_F`, slab solves with `D_F` on `u` with t links cut on
`{c1, c1+fw-1, fw-2, Lt-1}`):

```text
M   x = R_q S_F R_{w1} D x                    (x on p = w2)
M^dag y = R_p D^dag R_{w1} g5 S_F g5 R_q y      (y on q)
```

| interval | p = w2 (Pi) | w1 | q |
|---|---|---|---|
| F1 | c1 | c1+1 | c1+fw-1 |
| F0 | fw-1 | fw-2 | 0 |

Example c1 = 32, fw = 4: non-overlap p = {35, 0}, q = {32, 3}; overlap p = {32, 3},
q = {35, 0}. Requires `F0` to start at 0 and widths >= 3. The overlapping mode requires
`SIMPLE_FERM_STATE`, as the overlapping prop path already does. Testing uses simple
fermion states throughout.

## Algorithm (unchanged)

k = R + oversampling, qp = power iterations.

1. Save the QDP RNG state, seed with `seed`, draw k Gaussian fields, keep their values on
   the q slices (`Omega`), restore the RNG state.
2. `Y = M^dag Omega`.
3. qp times: `Y <- orth(Y)`; `Y <- M^dag (M Y)`.
4. `Q = orth(Y)` (throw if fewer than R columns survive).
5. `Z = M Q`; per slice `Z^dag Z = W Sigma^2 W^dag` (`SB::svd`, k x k); `phi = Q W[:, :R]`.

`orth` = Cholesky-QR through the Gram matrix, twice, dropping eigenvalues below
100 eps lambda_max. Cost: k (2 + 2 qp) slab solves, batched by `max_rhs`, shared by both
intervals. Default qp = 2: on the toy it matches the exact-SVD basis, while qp = 1 lost
up to 0.06 per quark at small R.

## Measurement 1: EDGE_SVD_PROJECTORS_SUPERB (new)

### XML

```xml
<elem>
  <Name>EDGE_SVD_PROJECTORS_SUPERB</Name>
  <Frequency>1</Frequency>
  <Param>
    <rank>512</rank>                                   <!-- R_max stored in the file -->
    <oversampling>16</oversampling>                     <!-- optional, default 16 -->
    <power_iterations>2</power_iterations>              <!-- optional, default 2 -->
    <seed>11</seed>                                     <!-- optional, default 11 -->
    <overlapping_factorization>true</overlapping_factorization>  <!-- default false -->
    <max_rhs>4</max_rhs>                                <!-- optional, default 8 -->
    <decay_dir>3</decay_dir>
    <Propagator>                                        <!-- same block as the prop measurement -->
      <FermionAction>...</FermionAction>
      <InvertParam>...</InvertParam>
    </Propagator>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <frozen_boundary_sidecar_file>...</frozen_boundary_sidecar_file>
    <frozen_boundary_child_id>0</frozen_boundary_child_id>
    <svd_basis_file>...</svd_basis_file>
  </NamedObject>
</elem>
```

- `rank` must satisfy 0 < rank <= 6 V3 (the spin half that the hop passes).
- `Propagator` is read with the existing `ChromaProp_t` reader; only `FermionAction` and
  `InvertParam` are used.

### Steps

1. Read the gauge field from the named object map; write the gauge observables
   (`MesPlq`) to the output XML, for information only.
2. Load the geometry with `loadFactorizedPropGeometry(sidecar, child_id, decay_dir)`.
3. Get the slice lists p, q and the two intermediate lists from the shared helper
   `edgeSVDSlices(geometry, overlapping, Lt)` (below). It runs the same layout checks
   as the prop measurement.
4. Build `fLinOp` from the uncut field (as the prop measurement does), and the slab
   solver: non-overlap `u_sep` via the shared helper `nonOverlapCutField`; overlap
   `u_F` with the four cuts above. One `SB::ChimeraSolver`.
5. Build `M`, `M^dag` from the formulas above and call `computeEdgeSVDBasis`.
6. Write the file (Storage format, below). Write a summary (rank, time, sigma_1,
   sigma_R/sigma_1 and orthonormality error per slice) to the output XML.

### File format

`SB::StorageTensor<Nd + 4, SB::ComplexD>`, order `cxyzXsnt` (the code's even-odd layout, the same as
`SB::asTensorView(LatticeFermion)`; x = Lx/2, X = 2), dimensions `latticeSize` with `n = rank` and `t = Lt`, `Sparse`
storage, `BlockChecksum`, `SharedFSFile`. Only the two projector slices are written:
column j of slice p at `t = p`, `n = j`. Metadata (XML, root `EdgeSVDBasis`, broadcast
from the master node like the perambulator metadata):

- `order`, `version` (1), `lattSize`, `decay_dir`
- `overlapping_factorization`, `frozen_width`, `cut1`
- `projector_slices`, `far_slices` (in matching order)
- `rank`, `oversampling`, `power_iterations`, `seed`
- `num_singular_values` (k) and `singular_values` (all k Ritz values per projector slice,
  slice-major), `orthonormality_error` (per projector slice)
- `Config_info` (gauge XML, for reference; not checked)
- `Propagator` (the `ChromaProp_t` block: FermionAction and InvertParam). The prop
  measurement reads the FermionAction back with `readXMLGroup` and compares it to its own
  with all whitespace removed.

At 32^3, one column is 12 V3 complex doubles = 6.3 MB per slice: 6.4 GB per config at
R_max = 512, 1.6 GB at 128.

## Measurement 2: PROP_AND_MATELEM_DISTILLATION_SUPERB (changes)

### XML

```xml
<projector_type>singular_vectors</projector_type>     <!-- distillation_vecs (default) | singular_vectors -->
<num_intermediate_projectors>128</num_intermediate_projectors>  <!-- R used, <= rank in the file -->
<svd_basis_file>...</svd_basis_file>                  <!-- required for singular_vectors -->
```

- Default (`projector_type` absent or `distillation_vecs`): unchanged, bit for bit.
- `svd_oversampling`, `svd_power_iterations`, `svd_seed` are removed from this
  measurement (they belong to measurement 1). There is no in-line computation:
  `singular_vectors` without `svd_basis_file` throws.
- Checks kept from the in-line version: R > 0 (non-overlap), zero
  `intermediate_projector_phase`, R <= 6 V3.

### Reading the basis

Replaces the in-line computation block (lambdas, `PP_F`, `computeEdgeSVDBasis` call).
Still before the `t_source` loop:

1. Open the file with `SB::StorageTensor<Nd + 4, SB::ComplexD>(file, true, "/EdgeSVDBasis/order")`.
2. Parse the metadata and check, throwing with a clear message on a mismatch:
   - `lattSize`, `decay_dir`, `overlapping_factorization`, `frozen_width`, `cut1`;
   - `projector_slices` equal the slices this run projects on (from `edgeSVDSlices`);
   - `rank` >= `num_intermediate_projectors`;
   - the `FermionAction` XML equals this run's, compared with all whitespace removed.
3. For each projector slice p: copy `n = 0..R-1`, `t = p` from the file into a slice
   tensor (`cxyztXns`, t extent 1, `SB::Complex`), stored in `sv_basis[p]`.
4. Write the file name, the R used and sigma_R/sigma_1 per slice to the output XML.

`projectOnSlicesSV` and the two call sites are unchanged.

## Code organisation

Move to `factorized_distillation.{h,cc}` (namespace `FactorizedDistillation`), shared by
both measurements:

- `EdgeSVDSlices edgeSVDSlices(const FactorizedPropGeometry&, bool overlapping, int Lt)`:
  returns p, a, b, q (one entry per interval, F1 then F0 for non-overlap, F0 then F1 for
  overlap, as now) and runs the geometry checks.
- `multi1d<LatticeColorMatrix> nonOverlapCutField(u, geometry, decay_dir)`: the existing
  `u_sep` loop. The prop measurement calls it too (same cuts, no behaviour change), so the
  two measurements cannot drift apart.
- `computeEdgeSVDBasis` and its helpers (`edgeSVDApply`, `edgeSVDOrthonormalize`,
  `edgeSVDGram`, `edgeSVDGramEig`, `edgeSVDTimesSmall`, `newEdgeSVDSliceTensor`), moved
  unchanged out of the prop `.cc`.
- `writeEdgeSVDBasis(file, basis, info, config_xml, prop)`, `readEdgeSVDBasisInfo(file)`,
  `checkEdgeSVDBasisInfo(info, slices, overlapping, decay_dir, R, fermact_xml)` and
  `readEdgeSVDBasis(file, slices, R)` (which also checks the layout dimensions and that the
  vectors read are orthonormal to 1e-5): the format lives in one place. `EdgeSVDBasisInfo`
  holds the metadata.

`projectOnSlicesSV` stays in the prop `.cc` (only that measurement applies `Pi`).

New measurement files follow the pattern of
`inline_factorized_baryon_tensor_distillation_superb_w.{cc,h}`: `Params` with
`read`/`write`, `InlineMeas::func`, `registerAll`, guarded by `BUILD_SB`; registered in
`inline_hadron_aggregate.cc`; added to `lib/Makefile.am` and `lib/CMakeLists.txt`.

## Known caveats

- `D_F` and `D_00` are built like the other cut operators (links zeroed before the clover
  term). This is fine for a basis (it only has to span the right directions).
- The basis is a deterministic truncation, so a bias remains. Check it against the overlap
  `Pi = 1` run.
- The file is tied to one configuration, one fermion action and one geometry. The checks
  above enforce the lattice, geometry, slices, rank and fermion action. Whether the file
  comes from the same gauge configuration is not checked (2026-10-26); matching basis files
  to configurations is up to the job scripts.

## Verification done

2026-10-25 (two-measurement version):

1. Compile: `factorized_distillation.cc`, `inline_edge_svd_projectors_superb_w.cc`,
   `inline_prop_and_matelem_distillation_superb_w.cc` and `inline_hadron_aggregate.cc` compile
   to object files (g++ / OpenMPI, QDP++ `eloy/localbinarydb` scalar, superbblas master,
   `BUILD_SB`, `-Wall`); no new warnings beyond patterns already in the file.
2. `fact_prop_testing/edge_svd_runtime_test/t_edge_svd_io.cc` (4^3x16 scalar, synthetic M),
   all passed:
   - `edgeSVDSlices` for both paths and `nonOverlapCutField` (exactly the links rooted on
     t_start and t_end - 1 of each interval);
   - write -> read: metadata exact (singular values to 1e-13), vectors exact for R and R/2,
     the file holds only the two projector slices (size = 2 x R x 12 V3 x 16 B + header);
   - checks: pass for the matching run (FermionAction with different whitespace); throw on
     a different mass, a file rank below R, the wrong path, and a read of more vectors than
     stored. (Re-run 2026-10-26 after removing the plaquette check: all passed.)

2026-10-22 (in-line version; the moved SVD code is unchanged):

3. Standalone test of the SVD helpers: full-rank sketch matches an exact `SB::svd` to 1e-14;
   Pi Pi = Pi; zero off-slice.
4. Toy: `M`/`M^dag` adjoint to 1e-15 in both paths; full rank reproduces `Pi = 1`; the
   default settings match the exact-SVD basis.

Not verified: multi-rank runs, the real ChimeraSolver slab solves, GPU builds.

## Validation to do

1. Default type: prop SDB identical to the current code on one config, both paths.
2. Basis job on a small lattice, then the prop run reading it (MPI, real solvers).
3. Full-rank check on a small lattice (R = 6 V3): the prop result equals `Pi = 1` to solver
   precision.
4. 32^3: basis-job time and slab iteration counts; singular-value decay; prop ratio vs R
   (64, 128, 256, 512) from one R_max = 512 file, against `distillation_vecs` at equal total
   rank and against overlap `Pi = 1`.

## Decisions

1. Overlap block: frozen-only `G_F`, not `G` through `PP_snk`.
2. Both paths supported. Non-overlap uses `M = T^dag`.
3. `num_intermediate_projectors` = spin-colour vectors per projector slice used by the
   prop run; the file can hold more.
4. (2026-10-25) The basis is computed by its own measurement and written to a file; the
   prop measurement reads it. File size is not a concern.
5. Randomized SVD with Gram-matrix orthonormalization; dense steps use `SB::svd`.
6. Default power iterations = 2.
7. (2026-10-25) Stout fermion states are out of scope; tests use simple fermion states.
8. (2026-10-25) The prop measurement only reads the basis from a file; no in-line fallback.
9. (2026-10-25) Any metadata mismatch throws, including the fermion action.
10. (2026-10-25) The file stores the code's even-odd layout (`cxyzXsnt`).
11. (2026-10-26) No plaquette / gauge-configuration check: the plaquette is neither stored
    in the file nor compared. The other metadata checks stay.

## Open

- XML tag name `projector_type` (the request said `Projector_Type`); lower snake case
  kept to match the other `NamedObject` tags.
