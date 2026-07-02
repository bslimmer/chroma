#include "chroma.h"
#include "util/gauge/gauge_subdomain_split.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace Chroma;

namespace
{
  const std::string MODE_PERIODIC_ALL_T = "PERIODIC_ALL_T";
  const std::string MODE_INTERIOR_NOWRAP_CHILD = "INTERIOR_NOWRAP_CHILD";
  const std::string MODE_PARENT_WINDOW = "PARENT_WINDOW";
  const std::string MODE_TWO_LEVEL_CROSS_DOMAIN = "TWO_LEVEL_CROSS_DOMAIN";
  const std::string MODE_TWO_LEVEL_OUTER_ENSEMBLE = "TWO_LEVEL_OUTER_ENSEMBLE";

  struct ExpectedCount
  {
    int index;
    int value;

    ExpectedCount() : index(-1), value(-1) {}
  };

  struct MeasurementResult
  {
    int stream_id;
    unsigned long update_no;
    multi1d<Double> timeslice_operator;
    multi1d<Double> correlator;
    multi1d<int> num_sources;

    MeasurementResult() : stream_id(0), update_no(0) {}
  };

  struct ParentTimePair
  {
    int delta_t_parent;
    int child0_local_t;
    int child1_local_t;
    int parent_t0;
    int parent_t1;

    ParentTimePair()
      : delta_t_parent(-1),
        child0_local_t(-1),
        child1_local_t(-1),
        parent_t0(-1),
        parent_t1(-1)
    {
    }
  };

  struct LocalCorrelatorSummary
  {
    multi1d<int> delta_t_child_set;
    multi1d<int> num_sources;
    multi1d<Double> mean;
    multi1d<Double> stddev;
    multi1d<Double> stderr;
  };

  struct ChildSummary
  {
    std::string outer_sample_id;
    int child_id;
    multi1d<int> child_nrow;
    multi1d<int> child_local_to_global_t;
    multi1d<GaugeSubdomainSplitInterval> frozen_local_intervals;
    int support_guard;
    multi1d<int> safe_slices;
    int stream_count;
    multi1d<int> stream_ids;
    multi1d<int> retained_measurements_per_stream;
    int retained_measurement_count;
    multi1d<Double> conditional_mean_timeslice_operator;
    LocalCorrelatorSummary child_local_summary;
    std::vector<MeasurementResult> measurements;

    ChildSummary()
      : child_id(-1),
        support_guard(0),
        stream_count(0),
        retained_measurement_count(0)
    {
    }
  };

  struct ParentWindowSummary
  {
    std::string outer_sample_id;
    unsigned long retained_update_no;
    multi1d<int> nrow;
    int decay_dir;
    multi1d<int> delta_t_parent_set;
    multi1d<int> pair_counts;
    multi1d<Double> timeslice_operator;
    multi1d<Double> correlator;
    std::vector<ParentTimePair> pairs;

    ParentWindowSummary() : retained_update_no(0), decay_dir(-1) {}
  };

  struct OuterSampleSummary
  {
    std::string outer_sample_id;
    multi1d<int> delta_t_parent_set;
    multi1d<int> pair_counts;
    multi1d<Double> two_level_correlator;
    multi1d<Double> parent_window_correlator;
    multi1d<Double> delta_correlator;
    std::vector<ParentTimePair> pairs;
    LocalCorrelatorSummary child0_local_summary;
    LocalCorrelatorSummary child1_local_summary;
  };

  struct OuterEnsembleSummary
  {
    int n_outer_samples;
    multi1d<int> delta_t_parent_set;
    multi1d<Double> mean_two_level_correlator;
    multi1d<Double> stderr_two_level_correlator;
    multi1d<Double> mean_parent_window_correlator;
    multi1d<Double> stderr_parent_window_correlator;
    multi1d<Double> mean_delta_correlator;
    multi1d<Double> stderr_delta_correlator;
    LocalCorrelatorSummary child0_local_summary;
    LocalCorrelatorSummary child1_local_summary;
    bool pass;

    OuterEnsembleSummary() : n_outer_samples(0), pass(false) {}
  };

  struct CheckerInput
  {
    std::vector<std::string> input_files;
    std::vector<std::string> outer_sample_summary_files;
    std::string summary_file;
    std::string csv_file;
    multi1d<int> nrow;
    int decay_dir;
    std::string mode;
    int support_guard;
    int max_dt;
    std::string sidecar_file;
    int child_id;
    int discard_updates;
    std::string outer_sample_id;
    unsigned long retained_update_no;
    multi1d<int> compare_delta_t_parent_set;
    std::string child0_summary_file;
    std::string child1_summary_file;
    std::string parent_summary_file;
    int min_outer_samples;
    int expected_measurements;
    int expected_timeslices;
    multi1d<ExpectedCount> expected_num_sources;
    multi1d<int> expected_safe_slices;
    multi1d<ExpectedCount> expected_pairs;

    CheckerInput()
      : decay_dir(-1),
        support_guard(0),
        max_dt(-1),
        child_id(-1),
        discard_updates(0),
        retained_update_no(0),
        min_outer_samples(-1),
        expected_measurements(-1),
        expected_timeslices(-1)
    {
    }
  };

  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_qactden_0mp_corr: " << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  bool hasPath(XMLReader& xml, const std::string& path)
  {
    return xml.count(path) != 0;
  }

  std::string makeOuterSampleId(unsigned long update_no)
  {
    std::ostringstream os;
    os << "outer_" << update_no;
    return os.str();
  }

  unsigned long parseOuterSampleUpdateNo(const std::string& outer_sample_id)
  {
    const std::string prefix = "outer_";
    check(outer_sample_id.substr(0, prefix.size()) == prefix,
          "outer_sample_id must start with outer_");

    std::istringstream is(outer_sample_id.substr(prefix.size()));
    unsigned long update_no = 0;
    is >> update_no;
    check(!is.fail(), "failed to parse update number from outer_sample_id");
    return update_no;
  }

  template <class T>
  multi1d<T> toMulti1d(const std::vector<T>& values)
  {
    multi1d<T> out(values.size());
    for (int i = 0; i < out.size(); ++i)
      out[i] = values[i];
    return out;
  }

  std::vector<int> toStdVector(const multi1d<int>& values)
  {
    std::vector<int> out(values.size());
    for (int i = 0; i < values.size(); ++i)
      out[i] = values[i];
    return out;
  }

  int indexOf(const multi1d<int>& values, int target)
  {
    for (int i = 0; i < values.size(); ++i)
    {
      if (values[i] == target)
        return i;
    }
    return -1;
  }

  bool sameIntList(const multi1d<int>& a, const multi1d<int>& b)
  {
    if (a.size() != b.size())
      return false;
    for (int i = 0; i < a.size(); ++i)
    {
      if (a[i] != b[i])
        return false;
    }
    return true;
  }

  bool sameIntervalList(const multi1d<GaugeSubdomainSplitInterval>& a,
                        const multi1d<GaugeSubdomainSplitInterval>& b)
  {
    if (a.size() != b.size())
      return false;
    for (int i = 0; i < a.size(); ++i)
    {
      if (a[i].t_start != b[i].t_start || a[i].t_end != b[i].t_end)
        return false;
    }
    return true;
  }

  int latticeVolume(const multi1d<int>& nrow)
  {
    int volume = 1;
    for (int mu = 0; mu < nrow.size(); ++mu)
      volume *= nrow[mu];
    return volume;
  }

  void read(XMLReader& xml, const std::string& path, ExpectedCount& expected)
  {
    XMLReader top(xml, path);

    if (hasPath(top, "dt"))
      read(top, "dt", expected.index);
    else if (hasPath(top, "delta_t_parent"))
      read(top, "delta_t_parent", expected.index);
    else if (hasPath(top, "delta_t_child"))
      read(top, "delta_t_child", expected.index);
    else if (hasPath(top, "index"))
      read(top, "index", expected.index);
    else
      fail("expected count entry is missing dt/delta_t_parent/delta_t_child/index");

    read(top, "value", expected.value);
  }

  void writeExpectedCountList(XMLWriter& xml,
                              const std::string& path,
                              const multi1d<int>& indices,
                              const multi1d<int>& values,
                              const std::string& index_tag)
  {
    check(indices.size() == values.size(),
          "writeExpectedCountList size mismatch");

    push(xml, path);
    for (int i = 0; i < indices.size(); ++i)
    {
      push(xml, "elem");
      write(xml, index_tag, indices[i]);
      write(xml, "value", values[i]);
      pop(xml);
    }
    pop(xml);
  }

  void readExpectedCountList(XMLReader& xml,
                             const std::string& path,
                             multi1d<ExpectedCount>& counts)
  {
    if (!hasPath(xml, path))
      return;

    XMLReader top(xml, path);
    counts.resize(top.count("elem"));
    for (int i = 0; i < counts.size(); ++i)
    {
      std::ostringstream xpath;
      xpath << "elem[" << (i + 1) << "]";
      read(top, xpath.str(), counts[i]);
    }
  }

  void readStringList(XMLReader& xml,
                      const std::string& path,
                      std::vector<std::string>& values)
  {
    if (!hasPath(xml, path))
      return;

    XMLReader top(xml, path);
    const int count = top.count("elem");
    values.resize(count);
    for (int i = 0; i < count; ++i)
    {
      std::ostringstream xpath;
      xpath << "elem[" << (i + 1) << "]";
      read(top, xpath.str(), values[i]);
    }
  }

  void readCheckerInput(XMLReader& xml, const std::string& path, CheckerInput& input)
  {
    XMLReader input_xml(xml, path);

    if (hasPath(input_xml, "Input"))
    {
      XMLReader input_block(input_xml, "Input");
      if (hasPath(input_block, "file"))
      {
        std::string file;
        read(input_block, "file", file);
        input.input_files.push_back(file);
      }
      readStringList(input_block, "files", input.input_files);
      readStringList(input_block, "outer_sample_summaries",
                     input.outer_sample_summary_files);
    }

    if (hasPath(input_xml, "Output"))
    {
      XMLReader output_xml(input_xml, "Output");
      if (hasPath(output_xml, "summary_file"))
        read(output_xml, "summary_file", input.summary_file);
      if (hasPath(output_xml, "csv_file"))
        read(output_xml, "csv_file", input.csv_file);
    }

    if (hasPath(input_xml, "Geometry"))
    {
      XMLReader geometry_xml(input_xml, "Geometry");
      if (hasPath(geometry_xml, "nrow"))
        read(geometry_xml, "nrow", input.nrow);
      if (hasPath(geometry_xml, "decay_dir"))
        read(geometry_xml, "decay_dir", input.decay_dir);
    }

    if (hasPath(input_xml, "Reducer"))
    {
      XMLReader reducer_xml(input_xml, "Reducer");
      read(reducer_xml, "mode", input.mode);
      if (hasPath(reducer_xml, "support_guard"))
        read(reducer_xml, "support_guard", input.support_guard);
      if (hasPath(reducer_xml, "max_dt"))
        read(reducer_xml, "max_dt", input.max_dt);
      if (hasPath(reducer_xml, "sidecar_file"))
        read(reducer_xml, "sidecar_file", input.sidecar_file);
      if (hasPath(reducer_xml, "child_id"))
        read(reducer_xml, "child_id", input.child_id);
      if (hasPath(reducer_xml, "discard_updates"))
        read(reducer_xml, "discard_updates", input.discard_updates);
      if (hasPath(reducer_xml, "outer_sample_id"))
        read(reducer_xml, "outer_sample_id", input.outer_sample_id);
      if (hasPath(reducer_xml, "retained_update_no"))
        read(reducer_xml, "retained_update_no", input.retained_update_no);
      if (hasPath(reducer_xml, "compare_delta_t_parent_set"))
        read(reducer_xml, "compare_delta_t_parent_set",
             input.compare_delta_t_parent_set);
      if (hasPath(reducer_xml, "min_outer_samples"))
        read(reducer_xml, "min_outer_samples", input.min_outer_samples);
    }

    if (hasPath(input_xml, "OuterSample"))
    {
      XMLReader outer_xml(input_xml, "OuterSample");
      if (hasPath(outer_xml, "id"))
        read(outer_xml, "id", input.outer_sample_id);
      if (hasPath(outer_xml, "sidecar_file"))
        read(outer_xml, "sidecar_file", input.sidecar_file);
    }

    if (hasPath(input_xml, "Child0"))
    {
      XMLReader child0_xml(input_xml, "Child0");
      if (hasPath(child0_xml, "summary_file"))
        read(child0_xml, "summary_file", input.child0_summary_file);
    }

    if (hasPath(input_xml, "Child1"))
    {
      XMLReader child1_xml(input_xml, "Child1");
      if (hasPath(child1_xml, "summary_file"))
        read(child1_xml, "summary_file", input.child1_summary_file);
    }

    if (hasPath(input_xml, "Parent"))
    {
      XMLReader parent_xml(input_xml, "Parent");
      if (hasPath(parent_xml, "summary_file"))
        read(parent_xml, "summary_file", input.parent_summary_file);
    }

    if (hasPath(input_xml, "Checks"))
    {
      XMLReader checks_xml(input_xml, "Checks");
      if (hasPath(checks_xml, "expected_measurements"))
        read(checks_xml, "expected_measurements", input.expected_measurements);
      if (hasPath(checks_xml, "expected_timeslices"))
        read(checks_xml, "expected_timeslices", input.expected_timeslices);
      if (hasPath(checks_xml, "expected_safe_slices"))
        read(checks_xml, "expected_safe_slices", input.expected_safe_slices);
      if (hasPath(checks_xml, "min_outer_samples"))
        read(checks_xml, "min_outer_samples", input.min_outer_samples);
      readExpectedCountList(checks_xml, "expected_num_sources",
                            input.expected_num_sources);
      readExpectedCountList(checks_xml, "expected_pairs",
                            input.expected_pairs);
    }
  }

  void checkExpectedCounts(const multi1d<ExpectedCount>& expected,
                           const multi1d<int>& counts,
                           const std::string& label)
  {
    for (int i = 0; i < expected.size(); ++i)
    {
      check(expected[i].index >= 0 && expected[i].index < counts.size(),
            label + " requested an out-of-range index");
      if (counts[expected[i].index] != expected[i].value)
      {
        std::ostringstream os;
        os << label << " check failed at index " << expected[i].index
           << ": expected " << expected[i].value
           << ", got " << counts[expected[i].index];
        fail(os.str());
      }
    }
  }

  void checkExpectedCountsForLabels(const multi1d<ExpectedCount>& expected,
                                    const multi1d<int>& labels,
                                    const multi1d<int>& counts,
                                    const std::string& value_label)
  {
    check(labels.size() == counts.size(),
          value_label + " label/count size mismatch");

    for (int i = 0; i < expected.size(); ++i)
    {
      const int idx = indexOf(labels, expected[i].index);
      check(idx >= 0,
            value_label + " requested a label not present in the reduced set");
      if (counts[idx] != expected[i].value)
      {
        std::ostringstream os;
        os << value_label << " check failed at label " << expected[i].index
           << ": expected " << expected[i].value
           << ", got " << counts[idx];
        fail(os.str());
      }
    }
  }

  void validateGeometry(const multi1d<int>& nrow, int decay_dir)
  {
    check(nrow.size() == Nd, "nrow size must match Nd");
    check(decay_dir >= 0 && decay_dir < Nd,
          "decay_dir must be in [0, Nd)");
  }

  void validateInput(const CheckerInput& input)
  {
    check(!input.mode.empty(), "Reducer/mode is required");

    if (input.mode == MODE_PERIODIC_ALL_T)
    {
      check(input.input_files.size() == 1,
            "PERIODIC_ALL_T expects exactly one Input/file");
      validateGeometry(input.nrow, input.decay_dir);
      if (input.max_dt < 0)
        fail("PERIODIC_ALL_T requires max_dt");
      check(input.max_dt <= input.nrow[input.decay_dir] / 2,
            "PERIODIC_ALL_T max_dt must be <= floor(T/2)");
    }
    else if (input.mode == MODE_INTERIOR_NOWRAP_CHILD)
    {
      check(!input.input_files.empty(),
            "INTERIOR_NOWRAP_CHILD expects Input/file or Input/files");
      validateGeometry(input.nrow, input.decay_dir);
      if (input.max_dt < 0)
        fail("INTERIOR_NOWRAP_CHILD requires max_dt");
      check(input.max_dt <= input.nrow[input.decay_dir] - 1,
            "INTERIOR_NOWRAP_CHILD max_dt must be <= T-1");
      check(!input.sidecar_file.empty(),
            "INTERIOR_NOWRAP_CHILD requires Reducer/sidecar_file");
      check(input.child_id == 0 || input.child_id == 1,
            "INTERIOR_NOWRAP_CHILD child_id must be 0 or 1");
      check(!input.outer_sample_id.empty(),
            "INTERIOR_NOWRAP_CHILD requires Reducer/outer_sample_id");
    }
    else if (input.mode == MODE_PARENT_WINDOW)
    {
      check(input.input_files.size() == 1,
            "PARENT_WINDOW expects exactly one Input/file");
      validateGeometry(input.nrow, input.decay_dir);
      check(!input.sidecar_file.empty(),
            "PARENT_WINDOW requires Reducer/sidecar_file");
      check(input.retained_update_no > 0,
            "PARENT_WINDOW requires Reducer/retained_update_no > 0");
      check(input.support_guard >= 0,
            "PARENT_WINDOW requires support_guard >= 0");
    }
    else if (input.mode == MODE_TWO_LEVEL_CROSS_DOMAIN)
    {
      check(!input.sidecar_file.empty(),
            "TWO_LEVEL_CROSS_DOMAIN requires OuterSample/sidecar_file");
      check(!input.child0_summary_file.empty(),
            "TWO_LEVEL_CROSS_DOMAIN requires Child0/summary_file");
      check(!input.child1_summary_file.empty(),
            "TWO_LEVEL_CROSS_DOMAIN requires Child1/summary_file");
      check(!input.parent_summary_file.empty(),
            "TWO_LEVEL_CROSS_DOMAIN requires Parent/summary_file");
    }
    else if (input.mode == MODE_TWO_LEVEL_OUTER_ENSEMBLE)
    {
      check(!input.outer_sample_summary_files.empty(),
            "TWO_LEVEL_OUTER_ENSEMBLE requires Input/outer_sample_summaries");
      if (input.min_outer_samples >= 0)
      {
        check(int(input.outer_sample_summary_files.size()) >= input.min_outer_samples,
              "outer sample summary list is shorter than min_outer_samples");
      }
    }
    else
    {
      fail("unsupported reducer mode: " + input.mode);
    }
  }

  void checkHmcGeometry(XMLReader& hmc_xml, const multi1d<int>& expected_nrow)
  {
    multi1d<int> hmc_nrow;
    read(hmc_xml, "/hmc/Input/Params/HMCTrj/nrow", hmc_nrow);
    check(hmc_nrow.size() == expected_nrow.size(), "HMC nrow size mismatch");

    for (int mu = 0; mu < expected_nrow.size(); ++mu)
    {
      if (hmc_nrow[mu] != expected_nrow[mu])
      {
        std::ostringstream os;
        os << "HMC nrow mismatch at direction " << mu
           << ": expected " << expected_nrow[mu]
           << ", got " << hmc_nrow[mu];
        fail(os.str());
      }
    }
  }

  GaugeSubdomainSplitPlan readSplitPlan(const std::string& sidecar_file)
  {
    XMLReader sidecar_in(sidecar_file);
    GaugeSubdomainSplitPlan plan;
    read(sidecar_in, "/GaugeSubdomainSplitInfo", plan);
    return plan;
  }

  const multi1d<int>& childNrow(const GaugeSubdomainSplitPlan& plan, int child_id)
  {
    return (child_id == 0) ? plan.child0_nrow : plan.child1_nrow;
  }

  const multi1d<int>& childLocalToGlobalT(const GaugeSubdomainSplitPlan& plan, int child_id)
  {
    return (child_id == 0) ? plan.child0_local_to_global_t
                           : plan.child1_local_to_global_t;
  }

  const multi1d<GaugeSubdomainSplitInterval>& childFrozenIntervals(const GaugeSubdomainSplitPlan& plan,
                                                                   int child_id)
  {
    return (child_id == 0) ? plan.child0_frozen_local_intervals
                           : plan.child1_frozen_local_intervals;
  }

  bool isFrozenLocalTime(int local_t,
                         const multi1d<GaugeSubdomainSplitInterval>& intervals)
  {
    for (int i = 0; i < intervals.size(); ++i)
    {
      if (local_t >= intervals[i].t_start && local_t <= intervals[i].t_end)
        return true;
    }
    return false;
  }

  multi1d<int> computeSafeSlices(const multi1d<GaugeSubdomainSplitInterval>& intervals,
                                 int child_extent,
                                 int support_guard)
  {
    std::vector<int> safe;

    for (int t = 0; t < child_extent; ++t)
    {
      bool excluded = false;
      for (int i = 0; i < intervals.size(); ++i)
      {
        const int start = std::max(0, intervals[i].t_start - support_guard);
        const int end = std::min(child_extent - 1, intervals[i].t_end + support_guard);
        if (t >= start && t <= end)
        {
          excluded = true;
          break;
        }
      }

      if (!excluded)
        safe.push_back(t);
    }

    return toMulti1d(safe);
  }

  std::vector<ParentTimePair>
  buildPairs(const GaugeSubdomainSplitPlan& plan, int support_guard)
  {
    const multi1d<int> safe0 =
      computeSafeSlices(plan.child0_frozen_local_intervals,
                        plan.child0_nrow[plan.param.t_dir],
                        support_guard);
    const multi1d<int> safe1 =
      computeSafeSlices(plan.child1_frozen_local_intervals,
                        plan.child1_nrow[plan.param.t_dir],
                        support_guard);

    std::vector<ParentTimePair> pairs;
    for (int i = 0; i < safe0.size(); ++i)
    {
      const int t0 = safe0[i];
      const int parent_t0 = plan.child0_local_to_global_t[t0];
      for (int j = 0; j < safe1.size(); ++j)
      {
        const int t1 = safe1[j];
        const int parent_t1 = plan.child1_local_to_global_t[t1];
        const int delta_t_parent = parent_t1 - parent_t0;
        if (delta_t_parent <= 0)
          continue;

        ParentTimePair pair;
        pair.delta_t_parent = delta_t_parent;
        pair.child0_local_t = t0;
        pair.child1_local_t = t1;
        pair.parent_t0 = parent_t0;
        pair.parent_t1 = parent_t1;
        pairs.push_back(pair);
      }
    }

    std::sort(pairs.begin(), pairs.end(),
              [](const ParentTimePair& a, const ParentTimePair& b) {
                if (a.delta_t_parent != b.delta_t_parent)
                  return a.delta_t_parent < b.delta_t_parent;
                if (a.child0_local_t != b.child0_local_t)
                  return a.child0_local_t < b.child0_local_t;
                return a.child1_local_t < b.child1_local_t;
              });

    return pairs;
  }

  multi1d<int> deriveDeltaSet(const std::vector<ParentTimePair>& pairs)
  {
    std::vector<int> deltas;
    for (std::size_t i = 0; i < pairs.size(); ++i)
    {
      if (deltas.empty() || deltas.back() != pairs[i].delta_t_parent)
        deltas.push_back(pairs[i].delta_t_parent);
    }
    return toMulti1d(deltas);
  }

  std::vector<ParentTimePair> filterPairs(const std::vector<ParentTimePair>& pairs,
                                          const multi1d<int>& delta_set)
  {
    if (delta_set.size() == 0)
      return pairs;

    std::vector<ParentTimePair> filtered;
    for (std::size_t i = 0; i < pairs.size(); ++i)
    {
      if (indexOf(delta_set, pairs[i].delta_t_parent) >= 0)
        filtered.push_back(pairs[i]);
    }
    return filtered;
  }

  multi1d<int> countPairsByDelta(const std::vector<ParentTimePair>& pairs,
                                 const multi1d<int>& delta_set)
  {
    multi1d<int> counts(delta_set.size());
    for (int i = 0; i < counts.size(); ++i)
      counts[i] = 0;

    for (std::size_t i = 0; i < pairs.size(); ++i)
    {
      const int idx = indexOf(delta_set, pairs[i].delta_t_parent);
      if (idx >= 0)
        counts[idx] += 1;
    }

    return counts;
  }

  void write(XMLWriter& xml, const std::string& path, const ParentTimePair& pair)
  {
    push(xml, path);
    write(xml, "delta_t_parent", pair.delta_t_parent);
    write(xml, "child0_local_t", pair.child0_local_t);
    write(xml, "child1_local_t", pair.child1_local_t);
    write(xml, "parent_t0", pair.parent_t0);
    write(xml, "parent_t1", pair.parent_t1);
    pop(xml);
  }

  void read(XMLReader& xml, const std::string& path, ParentTimePair& pair)
  {
    XMLReader top(xml, path);
    read(top, "delta_t_parent", pair.delta_t_parent);
    read(top, "child0_local_t", pair.child0_local_t);
    read(top, "child1_local_t", pair.child1_local_t);
    read(top, "parent_t0", pair.parent_t0);
    read(top, "parent_t1", pair.parent_t1);
  }

  void writePairList(XMLWriter& xml,
                     const std::string& path,
                     const std::vector<ParentTimePair>& pairs)
  {
    push(xml, path);
    for (std::size_t i = 0; i < pairs.size(); ++i)
      write(xml, "elem", pairs[i]);
    pop(xml);
  }

  void readPairList(XMLReader& xml,
                    const std::string& path,
                    std::vector<ParentTimePair>& pairs)
  {
    if (!hasPath(xml, path))
      return;

    XMLReader top(xml, path);
    const int n_pairs = top.count("elem");
    pairs.resize(n_pairs);
    for (int i = 0; i < n_pairs; ++i)
    {
      std::ostringstream xpath;
      xpath << "elem[" << (i + 1) << "]";
      read(top, xpath.str(), pairs[i]);
    }
  }

  void write(XMLWriter& xml,
             const std::string& path,
             const MeasurementResult& measurement)
  {
    push(xml, path);
    write(xml, "stream_id", measurement.stream_id);
    write(xml, "update_no", measurement.update_no);
    write(xml, "timeslice_operator", measurement.timeslice_operator);
    write(xml, "correlator", measurement.correlator);
    write(xml, "num_sources", measurement.num_sources);
    pop(xml);
  }

  void read(XMLReader& xml,
            const std::string& path,
            MeasurementResult& measurement)
  {
    XMLReader top(xml, path);
    read(top, "stream_id", measurement.stream_id);
    read(top, "update_no", measurement.update_no);
    read(top, "timeslice_operator", measurement.timeslice_operator);
    read(top, "correlator", measurement.correlator);
    read(top, "num_sources", measurement.num_sources);
  }

  void writeMeasurementList(XMLWriter& xml,
                            const std::string& path,
                            const std::vector<MeasurementResult>& measurements)
  {
    push(xml, path);
    for (std::size_t i = 0; i < measurements.size(); ++i)
      write(xml, "elem", measurements[i]);
    pop(xml);
  }

  void readMeasurementList(XMLReader& xml,
                           const std::string& path,
                           std::vector<MeasurementResult>& measurements)
  {
    if (!hasPath(xml, path))
      return;

    XMLReader top(xml, path);
    const int n_measurements = top.count("elem");
    measurements.resize(n_measurements);
    for (int i = 0; i < n_measurements; ++i)
    {
      std::ostringstream xpath;
      xpath << "elem[" << (i + 1) << "]";
      read(top, xpath.str(), measurements[i]);
    }
  }

  void write(XMLWriter& xml,
             const std::string& path,
             const LocalCorrelatorSummary& summary)
  {
    push(xml, path);
    write(xml, "delta_t_child_set", summary.delta_t_child_set);
    write(xml, "num_sources", summary.num_sources);
    write(xml, "mean", summary.mean);
    write(xml, "stddev", summary.stddev);
    write(xml, "stderr", summary.stderr);
    pop(xml);
  }

  void read(XMLReader& xml,
            const std::string& path,
            LocalCorrelatorSummary& summary)
  {
    XMLReader top(xml, path);
    read(top, "delta_t_child_set", summary.delta_t_child_set);
    read(top, "num_sources", summary.num_sources);
    read(top, "mean", summary.mean);
    read(top, "stddev", summary.stddev);
    read(top, "stderr", summary.stderr);
  }

  multi1d<Double> makeZeroDoubleArray(int size)
  {
    multi1d<Double> values(size);
    values = zero;
    return values;
  }

  multi1d<int> makeZeroIntArray(int size)
  {
    multi1d<int> values(size);
    for (int i = 0; i < size; ++i)
      values[i] = 0;
    return values;
  }

  MeasurementResult makeTimesliceMeasurement(int stream_id,
                                             unsigned long update_no,
                                             const multi1d<Double>& timeslice_operator,
                                             int correlator_size)
  {
    MeasurementResult result;
    result.stream_id = stream_id;
    result.update_no = update_no;
    result.timeslice_operator = timeslice_operator;
    result.correlator = makeZeroDoubleArray(correlator_size);
    result.num_sources = makeZeroIntArray(correlator_size);
    return result;
  }

  multi1d<Double> computeTimesliceOperator(XMLReader& qactden_xml,
                                           const multi1d<int>& nrow,
                                           int decay_dir)
  {
    const int length = nrow[decay_dir];
    const int volume = latticeVolume(nrow);

    multi1d<Double> timeslice_operator(length);
    timeslice_operator = zero;

    XMLReader olattice_xml(qactden_xml, "naiveTopCharge/OLattice");
    const int n_sites = olattice_xml.count("elem");
    check(n_sites == volume, "unexpected OLattice site count in naiveTopCharge");

    for (int site = 0; site < n_sites; ++site)
    {
      std::ostringstream xpath;
      xpath << "elem[" << (site + 1) << "]";

      double q = 0.0;
      read(olattice_xml, xpath.str(), q);

      const multi1d<int> coord = crtesn(site, nrow);
      timeslice_operator[coord[decay_dir]] += Double(q);
    }

    return timeslice_operator;
  }

  std::vector<MeasurementResult>
  collectRawMeasurements(XMLReader& hmc_xml,
                         const multi1d<int>& nrow,
                         int decay_dir,
                         int stream_id,
                         int discard_updates)
  {
    std::vector<MeasurementResult> measurements;

    XMLReader updates_xml(hmc_xml, "/hmc/doHMC/MCUpdates");
    const int n_updates = updates_xml.count("./elem");

    for (int i = 1; i <= n_updates; ++i)
    {
      std::ostringstream update_xpath;
      update_xpath << "./elem[" << i << "]/Update";
      XMLReader update_xml(updates_xml, update_xpath.str());

      if (!hasPath(update_xml, "InlineObservables"))
        continue;

      unsigned long update_no = 0;
      read(update_xml, "update_no", update_no);

      XMLReader observables_xml(update_xml, "InlineObservables");
      const int n_observables = observables_xml.count("./elem");

      for (int obs = 1; obs <= n_observables; ++obs)
      {
        std::ostringstream observable_xpath;
        observable_xpath << "./elem[" << obs << "]";
        XMLReader observable_xml(observables_xml, observable_xpath.str());

        if (!hasPath(observable_xml, "QActDen"))
          continue;

        XMLReader qactden_xml(observable_xml, "QActDen");

        unsigned long measurement_update_no = 0;
        read(qactden_xml, "update_no", measurement_update_no);
        check(measurement_update_no == update_no,
              "QActDen update number disagrees with enclosing HMC update");

        if (int(measurement_update_no) <= discard_updates)
          continue;

        const multi1d<Double> timeslice_operator =
          computeTimesliceOperator(qactden_xml, nrow, decay_dir);

        MeasurementResult measurement;
        measurement.stream_id = stream_id;
        measurement.update_no = measurement_update_no;
        measurement.timeslice_operator = timeslice_operator;
        measurements.push_back(measurement);
      }
    }

    return measurements;
  }

  MeasurementResult reducePeriodicMeasurement(const MeasurementResult& raw,
                                              int max_dt)
  {
    const int length = raw.timeslice_operator.size();
    MeasurementResult result =
      makeTimesliceMeasurement(raw.stream_id,
                               raw.update_no,
                               raw.timeslice_operator,
                               max_dt + 1);

    for (int dt = 0; dt <= max_dt; ++dt)
    {
      double accum = 0.0;
      for (int t0 = 0; t0 < length; ++t0)
      {
        const int t1 = (t0 + dt) % length;
        accum += toDouble(raw.timeslice_operator[t0]) *
                 toDouble(raw.timeslice_operator[t1]);
      }
      result.correlator[dt] = Double(accum / double(length));
      result.num_sources[dt] = length;
    }

    return result;
  }

  MeasurementResult reduceInteriorNowrapMeasurement(const MeasurementResult& raw,
                                                    int max_dt,
                                                    const multi1d<int>& safe_slices)
  {
    const int length = raw.timeslice_operator.size();
    MeasurementResult result =
      makeTimesliceMeasurement(raw.stream_id,
                               raw.update_no,
                               raw.timeslice_operator,
                               max_dt + 1);

    std::vector<int> safe = toStdVector(safe_slices);
    std::vector<bool> is_safe(length, false);
    for (std::size_t i = 0; i < safe.size(); ++i)
      is_safe[safe[i]] = true;

    for (int dt = 0; dt <= max_dt; ++dt)
    {
      double accum = 0.0;
      int count = 0;
      for (std::size_t i = 0; i < safe.size(); ++i)
      {
        const int t0 = safe[i];
        const int t1 = t0 + dt;
        if (t1 >= length)
          continue;
        if (!is_safe[t1])
          continue;

        accum += toDouble(raw.timeslice_operator[t0]) *
                 toDouble(raw.timeslice_operator[t1]);
        ++count;
      }

      if (count > 0)
        result.correlator[dt] = Double(accum / double(count));
      else
        result.correlator[dt] = Double(0.0);

      result.num_sources[dt] = count;
    }

    return result;
  }

  void sortMeasurements(std::vector<MeasurementResult>& measurements)
  {
    std::sort(measurements.begin(), measurements.end(),
              [](const MeasurementResult& a, const MeasurementResult& b) {
                if (a.stream_id != b.stream_id)
                  return a.stream_id < b.stream_id;
                return a.update_no < b.update_no;
              });
  }

  multi1d<Double> computeMeanTimesliceOperator(const std::vector<MeasurementResult>& measurements)
  {
    check(!measurements.empty(), "cannot compute mean timeslice operator with no measurements");

    const int length = measurements[0].timeslice_operator.size();
    multi1d<Double> mean(length);
    mean = zero;

    for (std::size_t i = 0; i < measurements.size(); ++i)
    {
      check(measurements[i].timeslice_operator.size() == length,
            "timeslice operator length mismatch");
      for (int t = 0; t < length; ++t)
        mean[t] += measurements[i].timeslice_operator[t];
    }

    for (int t = 0; t < length; ++t)
      mean[t] /= Double(measurements.size());

    return mean;
  }

  multi1d<Double> computeMeanCorrelator(const std::vector<MeasurementResult>& measurements)
  {
    check(!measurements.empty(), "cannot compute mean correlator with no measurements");

    const int length = measurements[0].correlator.size();
    multi1d<Double> mean(length);
    mean = zero;

    for (std::size_t i = 0; i < measurements.size(); ++i)
    {
      check(measurements[i].correlator.size() == length,
            "correlator length mismatch");
      for (int dt = 0; dt < length; ++dt)
        mean[dt] += measurements[i].correlator[dt];
    }

    for (int dt = 0; dt < length; ++dt)
      mean[dt] /= Double(measurements.size());

    return mean;
  }

  multi1d<Double> computeSampleStddevCorrelator(const std::vector<MeasurementResult>& measurements,
                                                const multi1d<Double>& mean)
  {
    multi1d<Double> stddev(mean.size());
    stddev = zero;

    if (measurements.size() <= 1)
      return stddev;

    for (std::size_t i = 0; i < measurements.size(); ++i)
    {
      for (int dt = 0; dt < mean.size(); ++dt)
      {
        const double diff = toDouble(measurements[i].correlator[dt] - mean[dt]);
        stddev[dt] += Double(diff * diff);
      }
    }

    for (int dt = 0; dt < mean.size(); ++dt)
      stddev[dt] = Double(std::sqrt(toDouble(stddev[dt]) /
                                    double(measurements.size() - 1)));

    return stddev;
  }

  multi1d<Double> computeSampleStderr(const multi1d<Double>& stddev, int count)
  {
    multi1d<Double> stderr(stddev.size());
    stderr = zero;

    if (count <= 0)
      return stderr;

    for (int i = 0; i < stddev.size(); ++i)
      stderr[i] = Double(toDouble(stddev[i]) / std::sqrt(double(count)));

    return stderr;
  }

  LocalCorrelatorSummary buildLocalCorrelatorSummary(const std::vector<MeasurementResult>& measurements)
  {
    check(!measurements.empty(), "cannot build local correlator summary without measurements");

    LocalCorrelatorSummary summary;
    const int size = measurements[0].correlator.size();
    summary.delta_t_child_set.resize(size);
    for (int i = 0; i < size; ++i)
      summary.delta_t_child_set[i] = i;
    summary.num_sources = measurements[0].num_sources;
    summary.mean = computeMeanCorrelator(measurements);
    summary.stddev = computeSampleStddevCorrelator(measurements, summary.mean);
    summary.stderr = computeSampleStderr(summary.stddev, measurements.size());
    return summary;
  }

  void writePeriodicSummary(const CheckerInput& input,
                            const std::vector<MeasurementResult>& measurements)
  {
    if (input.summary_file.empty())
      return;

    XMLFileWriter xml_out(input.summary_file);
    push(xml_out, "qactden_0mp_corr_summary");
    write(xml_out, "Mode", MODE_PERIODIC_ALL_T);

    push(xml_out, "Geometry");
    write(xml_out, "nrow", input.nrow);
    write(xml_out, "decay_dir", input.decay_dir);
    pop(xml_out);

    push(xml_out, "Reducer");
    write(xml_out, "max_dt", input.max_dt);
    pop(xml_out);

    write(xml_out, "measurement_count", int(measurements.size()));
    write(xml_out, "mean_correlator", computeMeanCorrelator(measurements));
    writeMeasurementList(xml_out, "Measurements", measurements);

    pop(xml_out);
    xml_out.close();
  }

  void writeChildSummary(const CheckerInput& input, const ChildSummary& summary)
  {
    if (input.summary_file.empty())
      return;

    XMLFileWriter xml_out(input.summary_file);
    push(xml_out, "qactden_0mp_corr_summary");
    write(xml_out, "Mode", MODE_INTERIOR_NOWRAP_CHILD);

    push(xml_out, "OuterSample");
    write(xml_out, "id", summary.outer_sample_id);
    write(xml_out, "update_no", parseOuterSampleUpdateNo(summary.outer_sample_id));
    pop(xml_out);

    push(xml_out, "Child");
    write(xml_out, "child_id", summary.child_id);
    write(xml_out, "nrow", summary.child_nrow);
    write(xml_out, "local_to_global_t", summary.child_local_to_global_t);
    write(xml_out, "frozen_local_intervals", summary.frozen_local_intervals);
    write(xml_out, "support_guard", summary.support_guard);
    write(xml_out, "safe_slices", summary.safe_slices);
    pop(xml_out);

    push(xml_out, "Reduction");
    write(xml_out, "stream_count", summary.stream_count);
    write(xml_out, "stream_ids", summary.stream_ids);
    write(xml_out, "retained_measurements_per_stream",
          summary.retained_measurements_per_stream);
    write(xml_out, "retained_measurement_count",
          summary.retained_measurement_count);
    pop(xml_out);

    writeMeasurementList(xml_out, "Measurements", summary.measurements);

    push(xml_out, "ConditionalMean");
    write(xml_out, "timeslice_operator", summary.conditional_mean_timeslice_operator);
    write(xml_out, "child_local_correlator", summary.child_local_summary.mean);
    pop(xml_out);

    push(xml_out, "ChildLocalCorrelator");
    write(xml_out, "delta_t_child_set", summary.child_local_summary.delta_t_child_set);
    write(xml_out, "num_sources", summary.child_local_summary.num_sources);
    write(xml_out, "stddev", summary.child_local_summary.stddev);
    write(xml_out, "stderr", summary.child_local_summary.stderr);
    pop(xml_out);

    pop(xml_out);
    xml_out.close();
  }

  ChildSummary readChildSummary(const std::string& filename)
  {
    XMLReader xml_in(filename);
    XMLReader summary_xml(xml_in, "/qactden_0mp_corr_summary");

    std::string mode;
    read(summary_xml, "Mode", mode);
    check(mode == MODE_INTERIOR_NOWRAP_CHILD,
          "summary file is not an INTERIOR_NOWRAP_CHILD summary: " + filename);

    ChildSummary summary;

    XMLReader outer_xml(summary_xml, "OuterSample");
    read(outer_xml, "id", summary.outer_sample_id);

    XMLReader child_xml(summary_xml, "Child");
    read(child_xml, "child_id", summary.child_id);
    read(child_xml, "nrow", summary.child_nrow);
    read(child_xml, "local_to_global_t", summary.child_local_to_global_t);
    read(child_xml, "frozen_local_intervals", summary.frozen_local_intervals);
    read(child_xml, "support_guard", summary.support_guard);
    read(child_xml, "safe_slices", summary.safe_slices);

    XMLReader reduction_xml(summary_xml, "Reduction");
    read(reduction_xml, "stream_count", summary.stream_count);
    read(reduction_xml, "stream_ids", summary.stream_ids);
    read(reduction_xml, "retained_measurements_per_stream",
         summary.retained_measurements_per_stream);
    read(reduction_xml, "retained_measurement_count",
         summary.retained_measurement_count);

    readMeasurementList(summary_xml, "Measurements", summary.measurements);

    XMLReader mean_xml(summary_xml, "ConditionalMean");
    read(mean_xml, "timeslice_operator",
         summary.conditional_mean_timeslice_operator);
    summary.child_local_summary.mean.resize(0);
    read(mean_xml, "child_local_correlator", summary.child_local_summary.mean);

    XMLReader local_xml(summary_xml, "ChildLocalCorrelator");
    read(local_xml, "delta_t_child_set",
         summary.child_local_summary.delta_t_child_set);
    read(local_xml, "num_sources", summary.child_local_summary.num_sources);
    read(local_xml, "stddev", summary.child_local_summary.stddev);
    read(local_xml, "stderr", summary.child_local_summary.stderr);

    return summary;
  }

  void writeParentWindowSummary(const CheckerInput& input,
                                const ParentWindowSummary& summary)
  {
    if (input.summary_file.empty())
      return;

    XMLFileWriter xml_out(input.summary_file);
    push(xml_out, "qactden_0mp_corr_summary");
    write(xml_out, "Mode", MODE_PARENT_WINDOW);

    push(xml_out, "OuterSample");
    write(xml_out, "id", summary.outer_sample_id);
    write(xml_out, "update_no", summary.retained_update_no);
    pop(xml_out);

    push(xml_out, "Geometry");
    write(xml_out, "nrow", summary.nrow);
    write(xml_out, "decay_dir", summary.decay_dir);
    pop(xml_out);

    push(xml_out, "PairSet");
    write(xml_out, "delta_t_parent_set", summary.delta_t_parent_set);
    write(xml_out, "pair_counts", summary.pair_counts);
    writePairList(xml_out, "Pairs", summary.pairs);
    pop(xml_out);

    push(xml_out, "Measurement");
    write(xml_out, "timeslice_operator", summary.timeslice_operator);
    write(xml_out, "correlator", summary.correlator);
    pop(xml_out);

    pop(xml_out);
    xml_out.close();
  }

  ParentWindowSummary readParentWindowSummary(const std::string& filename)
  {
    XMLReader xml_in(filename);
    XMLReader summary_xml(xml_in, "/qactden_0mp_corr_summary");

    std::string mode;
    read(summary_xml, "Mode", mode);
    check(mode == MODE_PARENT_WINDOW,
          "summary file is not a PARENT_WINDOW summary: " + filename);

    ParentWindowSummary summary;

    XMLReader outer_xml(summary_xml, "OuterSample");
    read(outer_xml, "id", summary.outer_sample_id);
    read(outer_xml, "update_no", summary.retained_update_no);

    XMLReader geometry_xml(summary_xml, "Geometry");
    read(geometry_xml, "nrow", summary.nrow);
    read(geometry_xml, "decay_dir", summary.decay_dir);

    XMLReader pair_xml(summary_xml, "PairSet");
    read(pair_xml, "delta_t_parent_set", summary.delta_t_parent_set);
    read(pair_xml, "pair_counts", summary.pair_counts);
    readPairList(pair_xml, "Pairs", summary.pairs);

    XMLReader measurement_xml(summary_xml, "Measurement");
    read(measurement_xml, "timeslice_operator", summary.timeslice_operator);
    read(measurement_xml, "correlator", summary.correlator);

    return summary;
  }

  void writeOuterSampleSummary(const CheckerInput& input,
                               const OuterSampleSummary& summary)
  {
    if (input.summary_file.empty())
      return;

    XMLFileWriter xml_out(input.summary_file);
    push(xml_out, "qactden_0mp_corr_summary");
    write(xml_out, "Mode", MODE_TWO_LEVEL_CROSS_DOMAIN);

    push(xml_out, "OuterSample");
    write(xml_out, "id", summary.outer_sample_id);
    write(xml_out, "update_no", parseOuterSampleUpdateNo(summary.outer_sample_id));
    pop(xml_out);

    push(xml_out, "PairSet");
    write(xml_out, "delta_t_parent_set", summary.delta_t_parent_set);
    write(xml_out, "pair_counts", summary.pair_counts);
    writePairList(xml_out, "Pairs", summary.pairs);
    pop(xml_out);

    push(xml_out, "Observables");
    write(xml_out, "two_level_correlator", summary.two_level_correlator);
    write(xml_out, "parent_window_correlator", summary.parent_window_correlator);
    write(xml_out, "delta_correlator", summary.delta_correlator);
    pop(xml_out);

    write(xml_out, "Child0LocalCorrelator", summary.child0_local_summary);
    write(xml_out, "Child1LocalCorrelator", summary.child1_local_summary);

    pop(xml_out);
    xml_out.close();
  }

  OuterSampleSummary readOuterSampleSummary(const std::string& filename)
  {
    XMLReader xml_in(filename);
    XMLReader summary_xml(xml_in, "/qactden_0mp_corr_summary");

    std::string mode;
    read(summary_xml, "Mode", mode);
    check(mode == MODE_TWO_LEVEL_CROSS_DOMAIN,
          "summary file is not a TWO_LEVEL_CROSS_DOMAIN summary: " + filename);

    OuterSampleSummary summary;

    XMLReader outer_xml(summary_xml, "OuterSample");
    read(outer_xml, "id", summary.outer_sample_id);

    XMLReader pair_xml(summary_xml, "PairSet");
    read(pair_xml, "delta_t_parent_set", summary.delta_t_parent_set);
    read(pair_xml, "pair_counts", summary.pair_counts);
    readPairList(pair_xml, "Pairs", summary.pairs);

    XMLReader observables_xml(summary_xml, "Observables");
    read(observables_xml, "two_level_correlator", summary.two_level_correlator);
    read(observables_xml, "parent_window_correlator",
         summary.parent_window_correlator);
    read(observables_xml, "delta_correlator", summary.delta_correlator);

    read(summary_xml, "Child0LocalCorrelator", summary.child0_local_summary);
    read(summary_xml, "Child1LocalCorrelator", summary.child1_local_summary);

    return summary;
  }

  void writeOuterEnsembleSummary(const CheckerInput& input,
                                 const OuterEnsembleSummary& summary)
  {
    if (input.summary_file.empty())
      return;

    XMLFileWriter xml_out(input.summary_file);
    push(xml_out, "qactden_0mp_corr_summary");
    write(xml_out, "Mode", MODE_TWO_LEVEL_OUTER_ENSEMBLE);

    push(xml_out, "OuterEnsemble");
    write(xml_out, "N0", summary.n_outer_samples);
    write(xml_out, "compare_delta_t_parent_set", summary.delta_t_parent_set);
    pop(xml_out);

    push(xml_out, "Means");
    write(xml_out, "two_level_correlator", summary.mean_two_level_correlator);
    write(xml_out, "parent_window_correlator",
          summary.mean_parent_window_correlator);
    write(xml_out, "delta_correlator", summary.mean_delta_correlator);
    pop(xml_out);

    push(xml_out, "Stderr");
    write(xml_out, "two_level_correlator", summary.stderr_two_level_correlator);
    write(xml_out, "parent_window_correlator",
          summary.stderr_parent_window_correlator);
    write(xml_out, "delta_correlator", summary.stderr_delta_correlator);
    pop(xml_out);

    write(xml_out, "Child0LocalCorrelator", summary.child0_local_summary);
    write(xml_out, "Child1LocalCorrelator", summary.child1_local_summary);

    push(xml_out, "Acceptance");
    write(xml_out, "pass", summary.pass);
    pop(xml_out);

    pop(xml_out);
    xml_out.close();
  }

  multi1d<Double> readArrayForDeltaSet(const multi1d<int>& source_delta_set,
                                       const multi1d<Double>& source_values,
                                       const multi1d<int>& target_delta_set,
                                       const std::string& label)
  {
    multi1d<Double> out(target_delta_set.size());
    out = zero;

    for (int i = 0; i < target_delta_set.size(); ++i)
    {
      const int idx = indexOf(source_delta_set, target_delta_set[i]);
      check(idx >= 0, label + " is missing a requested delta_t_parent");
      out[i] = source_values[idx];
    }

    return out;
  }

  multi1d<int> readIntArrayForDeltaSet(const multi1d<int>& source_delta_set,
                                       const multi1d<int>& source_values,
                                       const multi1d<int>& target_delta_set,
                                       const std::string& label)
  {
    multi1d<int> out(target_delta_set.size());
    for (int i = 0; i < target_delta_set.size(); ++i)
    {
      const int idx = indexOf(source_delta_set, target_delta_set[i]);
      check(idx >= 0, label + " is missing a requested delta_t_parent");
      out[i] = source_values[idx];
    }
    return out;
  }

  multi1d<Double> computeMeanSeries(const std::vector<multi1d<Double> >& series)
  {
    check(!series.empty(), "cannot compute mean over empty series set");

    const int length = series[0].size();
    multi1d<Double> mean(length);
    mean = zero;

    for (std::size_t i = 0; i < series.size(); ++i)
    {
      check(series[i].size() == length, "series length mismatch");
      for (int j = 0; j < length; ++j)
        mean[j] += series[i][j];
    }

    for (int j = 0; j < length; ++j)
      mean[j] /= Double(series.size());

    return mean;
  }

  multi1d<Double> computeStddevSeries(const std::vector<multi1d<Double> >& series,
                                      const multi1d<Double>& mean)
  {
    multi1d<Double> stddev(mean.size());
    stddev = zero;

    if (series.size() <= 1)
      return stddev;

    for (std::size_t i = 0; i < series.size(); ++i)
    {
      for (int j = 0; j < mean.size(); ++j)
      {
        const double diff = toDouble(series[i][j] - mean[j]);
        stddev[j] += Double(diff * diff);
      }
    }

    for (int j = 0; j < mean.size(); ++j)
      stddev[j] = Double(std::sqrt(toDouble(stddev[j]) /
                                   double(series.size() - 1)));

    return stddev;
  }

  void runMeasurementChecks(const CheckerInput& input,
                            const std::vector<MeasurementResult>& measurements,
                            const multi1d<int>* expected_safe_slices)
  {
    if (input.expected_measurements >= 0)
    {
      check(int(measurements.size()) == input.expected_measurements,
            "unexpected number of QActDen measurements");
    }

    check(!measurements.empty(), "no QActDen measurements found");

    for (std::size_t i = 0; i < measurements.size(); ++i)
    {
      if (input.expected_timeslices >= 0)
      {
        check(measurements[i].timeslice_operator.size() == input.expected_timeslices,
              "unexpected timeslice operator length");
      }

      if (input.expected_num_sources.size() > 0)
      {
        checkExpectedCounts(input.expected_num_sources,
                            measurements[i].num_sources,
                            "num_sources");
      }
    }

    if (expected_safe_slices != 0 && input.expected_safe_slices.size() > 0)
    {
      check(sameIntList(*expected_safe_slices, input.expected_safe_slices),
            "safe_slices check failed");
    }
  }

  MeasurementResult findMeasurementByUpdateNo(const std::vector<MeasurementResult>& measurements,
                                              unsigned long update_no)
  {
    bool found = false;
    MeasurementResult match;

    for (std::size_t i = 0; i < measurements.size(); ++i)
    {
      if (measurements[i].update_no == update_no)
      {
        check(!found, "duplicate retained parent QACTDEN update number");
        match = measurements[i];
        found = true;
      }
    }

    check(found, "failed to locate retained QACTDEN measurement at requested update");
    return match;
  }

  ParentWindowSummary buildParentWindowSummary(const CheckerInput& input)
  {
    GaugeSubdomainSplitPlan plan = readSplitPlan(input.sidecar_file);
    check(sameIntList(plan.parent_nrow, input.nrow),
          "PARENT_WINDOW nrow does not match split sidecar parent_nrow");

    XMLReader hmc_xml(input.input_files[0]);
    checkHmcGeometry(hmc_xml, input.nrow);

    std::vector<MeasurementResult> raw =
      collectRawMeasurements(hmc_xml, input.nrow, input.decay_dir, 0, -1);
    sortMeasurements(raw);

    const MeasurementResult raw_measurement =
      findMeasurementByUpdateNo(raw, input.retained_update_no);

    ParentWindowSummary summary;
    summary.retained_update_no = input.retained_update_no;
    summary.outer_sample_id = input.outer_sample_id.empty()
      ? makeOuterSampleId(input.retained_update_no)
      : input.outer_sample_id;
    check(summary.outer_sample_id == makeOuterSampleId(input.retained_update_no),
          "PARENT_WINDOW outer_sample_id must be outer_<retained_update_no>");
    summary.nrow = input.nrow;
    summary.decay_dir = input.decay_dir;
    summary.timeslice_operator = raw_measurement.timeslice_operator;

    const std::vector<ParentTimePair> all_pairs =
      buildPairs(plan, input.support_guard);
    const multi1d<int> derived_delta_set = deriveDeltaSet(all_pairs);
    summary.delta_t_parent_set =
      (input.compare_delta_t_parent_set.size() > 0)
        ? input.compare_delta_t_parent_set
        : derived_delta_set;
    summary.pairs = filterPairs(all_pairs, summary.delta_t_parent_set);
    summary.pair_counts = countPairsByDelta(summary.pairs, summary.delta_t_parent_set);
    summary.correlator.resize(summary.delta_t_parent_set.size());
    summary.correlator = zero;

    for (int i = 0; i < summary.delta_t_parent_set.size(); ++i)
    {
      const int delta_t_parent = summary.delta_t_parent_set[i];
      const int pair_count = summary.pair_counts[i];
      if (pair_count == 0)
        continue;

      double accum = 0.0;
      for (std::size_t j = 0; j < summary.pairs.size(); ++j)
      {
        if (summary.pairs[j].delta_t_parent != delta_t_parent)
          continue;
        accum += toDouble(summary.timeslice_operator[summary.pairs[j].parent_t1]) *
                 toDouble(summary.timeslice_operator[summary.pairs[j].parent_t0]);
      }
      summary.correlator[i] = Double(accum / double(pair_count));
    }

    if (input.expected_pairs.size() > 0)
      checkExpectedCountsForLabels(input.expected_pairs,
                                   summary.delta_t_parent_set,
                                   summary.pair_counts,
                                   "expected_pairs");

    writeParentWindowSummary(input, summary);
    return summary;
  }

  ChildSummary runChildMode(const CheckerInput& input)
  {
    GaugeSubdomainSplitPlan plan = readSplitPlan(input.sidecar_file);
    const multi1d<int>& expected_child_nrow = childNrow(plan, input.child_id);
    const multi1d<int>& local_to_global_t = childLocalToGlobalT(plan, input.child_id);
    const multi1d<GaugeSubdomainSplitInterval>& frozen_intervals =
      childFrozenIntervals(plan, input.child_id);

    check(sameIntList(expected_child_nrow, input.nrow),
          "INTERIOR_NOWRAP_CHILD nrow does not match split sidecar child_nrow");

    const multi1d<int> safe_slices =
      computeSafeSlices(frozen_intervals,
                        expected_child_nrow[plan.param.t_dir],
                        input.support_guard);

    std::vector<MeasurementResult> reduced;
    std::vector<int> retained_per_stream(input.input_files.size(), 0);

    for (std::size_t stream = 0; stream < input.input_files.size(); ++stream)
    {
      XMLReader hmc_xml(input.input_files[stream]);
      checkHmcGeometry(hmc_xml, input.nrow);

      const std::vector<MeasurementResult> raw =
        collectRawMeasurements(hmc_xml,
                               input.nrow,
                               input.decay_dir,
                               int(stream),
                               input.discard_updates);

      retained_per_stream[stream] = raw.size();
      for (std::size_t i = 0; i < raw.size(); ++i)
        reduced.push_back(reduceInteriorNowrapMeasurement(raw[i], input.max_dt, safe_slices));
    }

    sortMeasurements(reduced);
    runMeasurementChecks(input, reduced, &safe_slices);

    if (retained_per_stream.size() > 1)
    {
      for (std::size_t i = 1; i < retained_per_stream.size(); ++i)
      {
        check(retained_per_stream[i] == retained_per_stream[0],
              "all child streams must retain the same number of measurements");
      }
    }

    ChildSummary summary;
    summary.outer_sample_id = input.outer_sample_id;
    summary.child_id = input.child_id;
    summary.child_nrow = expected_child_nrow;
    summary.child_local_to_global_t = local_to_global_t;
    summary.frozen_local_intervals = frozen_intervals;
    summary.support_guard = input.support_guard;
    summary.safe_slices = safe_slices;
    summary.stream_count = input.input_files.size();
    summary.stream_ids.resize(summary.stream_count);
    summary.retained_measurements_per_stream.resize(summary.stream_count);
    for (int i = 0; i < summary.stream_count; ++i)
    {
      summary.stream_ids[i] = i;
      summary.retained_measurements_per_stream[i] = retained_per_stream[i];
    }
    summary.retained_measurement_count = reduced.size();
    summary.measurements = reduced;
    summary.conditional_mean_timeslice_operator =
      computeMeanTimesliceOperator(reduced);
    summary.child_local_summary = buildLocalCorrelatorSummary(reduced);

    writeChildSummary(input, summary);
    return summary;
  }

  void writeMeasurementCsv(const CheckerInput& input,
                           const std::vector<MeasurementResult>& measurements,
                           const std::string& correlator_series)
  {
    if (input.csv_file.empty())
      return;

    std::ofstream out(input.csv_file.c_str());
    check(out.good(), "failed to open CSV output file");

    out << "stream_id,update_no,series,index,value,num_sources\n";
    out << std::setprecision(17);

    for (std::size_t i = 0; i < measurements.size(); ++i)
    {
      for (int t = 0; t < measurements[i].timeslice_operator.size(); ++t)
      {
        out << measurements[i].stream_id
            << "," << measurements[i].update_no
            << ",O_q,"
            << t
            << "," << toDouble(measurements[i].timeslice_operator[t])
            << ",\n";
      }

      for (int dt = 0; dt < measurements[i].correlator.size(); ++dt)
      {
        out << measurements[i].stream_id
            << "," << measurements[i].update_no
            << "," << correlator_series
            << "," << dt
            << "," << toDouble(measurements[i].correlator[dt])
            << "," << measurements[i].num_sources[dt]
            << "\n";
      }
    }
  }

  void runPeriodicMode(const CheckerInput& input)
  {
    XMLReader hmc_xml(input.input_files[0]);
    checkHmcGeometry(hmc_xml, input.nrow);

    const std::vector<MeasurementResult> raw =
      collectRawMeasurements(hmc_xml, input.nrow, input.decay_dir, 0, -1);

    std::vector<MeasurementResult> reduced;
    reduced.reserve(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i)
      reduced.push_back(reducePeriodicMeasurement(raw[i], input.max_dt));

    sortMeasurements(reduced);
    runMeasurementChecks(input, reduced, 0);
    writePeriodicSummary(input, reduced);
    writeMeasurementCsv(input, reduced, "C_cfg");

    QDPIO::cout << "t_qactden_0mp_corr: passed" << std::endl;
    QDPIO::cout << "  mode: " << input.mode << std::endl;
    QDPIO::cout << "  measurements: " << reduced.size() << std::endl;
    QDPIO::cout << "  hmc_xml: " << input.input_files[0] << std::endl;
  }

  void runChildSummaryMode(const CheckerInput& input)
  {
    const ChildSummary summary = runChildMode(input);
    writeMeasurementCsv(input, summary.measurements, "C_local");

    QDPIO::cout << "t_qactden_0mp_corr: passed" << std::endl;
    QDPIO::cout << "  mode: " << input.mode << std::endl;
    QDPIO::cout << "  outer_sample_id: " << summary.outer_sample_id << std::endl;
    QDPIO::cout << "  child_id: " << summary.child_id << std::endl;
    QDPIO::cout << "  retained_measurements: "
                << summary.retained_measurement_count << std::endl;
  }

  void writeParentWindowCsv(const CheckerInput& input,
                            const ParentWindowSummary& summary)
  {
    if (input.csv_file.empty())
      return;

    std::ofstream out(input.csv_file.c_str());
    check(out.good(), "failed to open CSV output file");
    out << "delta_t_parent,series,child0_local_t,child1_local_t,parent_t0,parent_t1,value,pair_count\n";
    out << std::setprecision(17);

    for (std::size_t i = 0; i < summary.pairs.size(); ++i)
    {
      out << summary.pairs[i].delta_t_parent
          << ",pair,"
          << summary.pairs[i].child0_local_t
          << "," << summary.pairs[i].child1_local_t
          << "," << summary.pairs[i].parent_t0
          << "," << summary.pairs[i].parent_t1
          << ",,\n";
    }

    for (int i = 0; i < summary.delta_t_parent_set.size(); ++i)
    {
      out << summary.delta_t_parent_set[i]
          << ",parent_window,,,,,"
          << toDouble(summary.correlator[i])
          << "," << summary.pair_counts[i]
          << "\n";
    }
  }

  void runParentWindowMode(const CheckerInput& input)
  {
    const ParentWindowSummary summary = buildParentWindowSummary(input);
    writeParentWindowCsv(input, summary);

    QDPIO::cout << "t_qactden_0mp_corr: passed" << std::endl;
    QDPIO::cout << "  mode: " << input.mode << std::endl;
    QDPIO::cout << "  outer_sample_id: " << summary.outer_sample_id << std::endl;
    QDPIO::cout << "  retained_update_no: " << summary.retained_update_no << std::endl;
  }

  void checkChildSummaryConsistency(const ChildSummary& child,
                                    const GaugeSubdomainSplitPlan& plan,
                                    int expected_child_id,
                                    const std::string& expected_outer_sample_id)
  {
    check(child.child_id == expected_child_id,
          "child summary child_id mismatch");
    check(child.outer_sample_id == expected_outer_sample_id,
          "child summary outer_sample_id mismatch");
    check(sameIntList(child.child_nrow, childNrow(plan, expected_child_id)),
          "child summary nrow mismatch with split sidecar");
    check(sameIntList(child.child_local_to_global_t,
                      childLocalToGlobalT(plan, expected_child_id)),
          "child summary local_to_global_t mismatch with split sidecar");
    check(sameIntervalList(child.frozen_local_intervals,
                           childFrozenIntervals(plan, expected_child_id)),
          "child summary frozen_local_intervals mismatch with split sidecar");
    check(sameIntList(child.safe_slices,
                      computeSafeSlices(childFrozenIntervals(plan, expected_child_id),
                                        childNrow(plan, expected_child_id)[plan.param.t_dir],
                                        child.support_guard)),
          "child summary safe_slices mismatch with split sidecar");
  }

  void checkPairListConsistency(const std::vector<ParentTimePair>& a,
                                const std::vector<ParentTimePair>& b,
                                const std::string& label)
  {
    check(a.size() == b.size(), label + " pair list size mismatch");
    for (std::size_t i = 0; i < a.size(); ++i)
    {
      check(a[i].delta_t_parent == b[i].delta_t_parent &&
            a[i].child0_local_t == b[i].child0_local_t &&
            a[i].child1_local_t == b[i].child1_local_t &&
            a[i].parent_t0 == b[i].parent_t0 &&
            a[i].parent_t1 == b[i].parent_t1,
            label + " pair list mismatch");
    }
  }

  OuterSampleSummary buildOuterSampleSummary(const CheckerInput& input)
  {
    GaugeSubdomainSplitPlan plan = readSplitPlan(input.sidecar_file);
    ChildSummary child0 = readChildSummary(input.child0_summary_file);
    ChildSummary child1 = readChildSummary(input.child1_summary_file);
    ParentWindowSummary parent = readParentWindowSummary(input.parent_summary_file);

    std::string outer_sample_id = input.outer_sample_id.empty()
      ? parent.outer_sample_id
      : input.outer_sample_id;

    check(parent.outer_sample_id == outer_sample_id,
          "parent summary outer_sample_id mismatch");
    checkChildSummaryConsistency(child0, plan, 0, outer_sample_id);
    checkChildSummaryConsistency(child1, plan, 1, outer_sample_id);
    check(child0.support_guard == child1.support_guard,
          "child summary support_guard mismatch");

    const std::vector<ParentTimePair> all_pairs =
      buildPairs(plan, child0.support_guard);
    multi1d<int> delta_t_parent_set =
      (input.compare_delta_t_parent_set.size() > 0)
        ? input.compare_delta_t_parent_set
        : parent.delta_t_parent_set;
    std::vector<ParentTimePair> pairs = filterPairs(all_pairs, delta_t_parent_set);
    const multi1d<int> pair_counts = countPairsByDelta(pairs, delta_t_parent_set);

    checkExpectedCountsForLabels(input.expected_pairs,
                                 delta_t_parent_set,
                                 pair_counts,
                                 "expected_pairs");
    checkPairListConsistency(pairs, filterPairs(parent.pairs, delta_t_parent_set),
                             "parent summary");

    OuterSampleSummary summary;
    summary.outer_sample_id = outer_sample_id;
    summary.delta_t_parent_set = delta_t_parent_set;
    summary.pair_counts = pair_counts;
    summary.pairs = pairs;
    summary.two_level_correlator.resize(delta_t_parent_set.size());
    summary.two_level_correlator = zero;
    summary.parent_window_correlator =
      readArrayForDeltaSet(parent.delta_t_parent_set,
                           parent.correlator,
                           delta_t_parent_set,
                           "parent window correlator");
    summary.delta_correlator.resize(delta_t_parent_set.size());
    summary.delta_correlator = zero;

    for (int i = 0; i < delta_t_parent_set.size(); ++i)
    {
      const int delta_t_parent = delta_t_parent_set[i];
      const int pair_count = pair_counts[i];
      if (pair_count == 0)
        continue;

      double accum = 0.0;
      for (std::size_t j = 0; j < pairs.size(); ++j)
      {
        if (pairs[j].delta_t_parent != delta_t_parent)
          continue;
        accum += toDouble(child1.conditional_mean_timeslice_operator[pairs[j].child1_local_t]) *
                 toDouble(child0.conditional_mean_timeslice_operator[pairs[j].child0_local_t]);
      }

      summary.two_level_correlator[i] = Double(accum / double(pair_count));
      summary.delta_correlator[i] =
        summary.two_level_correlator[i] - summary.parent_window_correlator[i];
    }

    summary.child0_local_summary = child0.child_local_summary;
    summary.child1_local_summary = child1.child_local_summary;

    writeOuterSampleSummary(input, summary);
    return summary;
  }

  void writeOuterSampleCsv(const CheckerInput& input,
                           const OuterSampleSummary& summary)
  {
    if (input.csv_file.empty())
      return;

    std::ofstream out(input.csv_file.c_str());
    check(out.good(), "failed to open CSV output file");
    out << "delta_t_parent,series,child0_local_t,child1_local_t,parent_t0,parent_t1,value,pair_count\n";
    out << std::setprecision(17);

    for (std::size_t i = 0; i < summary.pairs.size(); ++i)
    {
      out << summary.pairs[i].delta_t_parent
          << ",pair,"
          << summary.pairs[i].child0_local_t
          << "," << summary.pairs[i].child1_local_t
          << "," << summary.pairs[i].parent_t0
          << "," << summary.pairs[i].parent_t1
          << ",,\n";
    }

    for (int i = 0; i < summary.delta_t_parent_set.size(); ++i)
    {
      out << summary.delta_t_parent_set[i]
          << ",two_level,,,,,"
          << toDouble(summary.two_level_correlator[i])
          << "," << summary.pair_counts[i]
          << "\n";
      out << summary.delta_t_parent_set[i]
          << ",parent_window,,,,,"
          << toDouble(summary.parent_window_correlator[i])
          << "," << summary.pair_counts[i]
          << "\n";
      out << summary.delta_t_parent_set[i]
          << ",delta,,,,,"
          << toDouble(summary.delta_correlator[i])
          << "," << summary.pair_counts[i]
          << "\n";
    }
  }

  void runOuterSampleMode(const CheckerInput& input)
  {
    const OuterSampleSummary summary = buildOuterSampleSummary(input);
    writeOuterSampleCsv(input, summary);

    QDPIO::cout << "t_qactden_0mp_corr: passed" << std::endl;
    QDPIO::cout << "  mode: " << input.mode << std::endl;
    QDPIO::cout << "  outer_sample_id: " << summary.outer_sample_id << std::endl;
  }

  OuterEnsembleSummary buildOuterEnsembleSummary(const CheckerInput& input)
  {
    std::vector<OuterSampleSummary> summaries;
    summaries.reserve(input.outer_sample_summary_files.size());

    for (std::size_t i = 0; i < input.outer_sample_summary_files.size(); ++i)
      summaries.push_back(readOuterSampleSummary(input.outer_sample_summary_files[i]));

    std::sort(summaries.begin(), summaries.end(),
              [](const OuterSampleSummary& a, const OuterSampleSummary& b) {
                return parseOuterSampleUpdateNo(a.outer_sample_id) <
                       parseOuterSampleUpdateNo(b.outer_sample_id);
              });

    if (input.min_outer_samples >= 0)
    {
      check(int(summaries.size()) >= input.min_outer_samples,
            "outer sample count is below min_outer_samples");
    }

    check(!summaries.empty(), "no outer sample summaries were provided");

    const multi1d<int> delta_t_parent_set =
      (input.compare_delta_t_parent_set.size() > 0)
        ? input.compare_delta_t_parent_set
        : summaries[0].delta_t_parent_set;

    std::vector<multi1d<Double> > two_level_values;
    std::vector<multi1d<Double> > parent_values;
    std::vector<multi1d<Double> > delta_values;
    std::vector<multi1d<Double> > child0_local_values;
    std::vector<multi1d<Double> > child1_local_values;

    const multi1d<int> child0_delta_t_child_set =
      summaries[0].child0_local_summary.delta_t_child_set;
    const multi1d<int> child1_delta_t_child_set =
      summaries[0].child1_local_summary.delta_t_child_set;
    const multi1d<int> child0_num_sources =
      summaries[0].child0_local_summary.num_sources;
    const multi1d<int> child1_num_sources =
      summaries[0].child1_local_summary.num_sources;

    for (std::size_t i = 0; i < summaries.size(); ++i)
    {
      two_level_values.push_back(readArrayForDeltaSet(summaries[i].delta_t_parent_set,
                                                      summaries[i].two_level_correlator,
                                                      delta_t_parent_set,
                                                      "two-level correlator"));
      parent_values.push_back(readArrayForDeltaSet(summaries[i].delta_t_parent_set,
                                                   summaries[i].parent_window_correlator,
                                                   delta_t_parent_set,
                                                   "parent-window correlator"));
      delta_values.push_back(readArrayForDeltaSet(summaries[i].delta_t_parent_set,
                                                  summaries[i].delta_correlator,
                                                  delta_t_parent_set,
                                                  "delta correlator"));

      check(sameIntList(summaries[i].child0_local_summary.delta_t_child_set,
                        child0_delta_t_child_set),
            "child0 local delta_t_child_set mismatch across outer samples");
      check(sameIntList(summaries[i].child1_local_summary.delta_t_child_set,
                        child1_delta_t_child_set),
            "child1 local delta_t_child_set mismatch across outer samples");
      check(sameIntList(summaries[i].child0_local_summary.num_sources,
                        child0_num_sources),
            "child0 local num_sources mismatch across outer samples");
      check(sameIntList(summaries[i].child1_local_summary.num_sources,
                        child1_num_sources),
            "child1 local num_sources mismatch across outer samples");

      child0_local_values.push_back(summaries[i].child0_local_summary.mean);
      child1_local_values.push_back(summaries[i].child1_local_summary.mean);
    }

    OuterEnsembleSummary summary;
    summary.n_outer_samples = summaries.size();
    summary.delta_t_parent_set = delta_t_parent_set;
    summary.mean_two_level_correlator = computeMeanSeries(two_level_values);
    summary.stderr_two_level_correlator =
      computeSampleStderr(computeStddevSeries(two_level_values,
                                             summary.mean_two_level_correlator),
                          two_level_values.size());
    summary.mean_parent_window_correlator = computeMeanSeries(parent_values);
    summary.stderr_parent_window_correlator =
      computeSampleStderr(computeStddevSeries(parent_values,
                                             summary.mean_parent_window_correlator),
                          parent_values.size());
    summary.mean_delta_correlator = computeMeanSeries(delta_values);

    const multi1d<Double> delta_stddev =
      computeStddevSeries(delta_values, summary.mean_delta_correlator);
    summary.stderr_delta_correlator =
      computeSampleStderr(delta_stddev, delta_values.size());

    summary.child0_local_summary.delta_t_child_set = child0_delta_t_child_set;
    summary.child0_local_summary.num_sources = child0_num_sources;
    summary.child0_local_summary.mean = computeMeanSeries(child0_local_values);
    summary.child0_local_summary.stddev =
      computeStddevSeries(child0_local_values, summary.child0_local_summary.mean);
    summary.child0_local_summary.stderr =
      computeSampleStderr(summary.child0_local_summary.stddev,
                          child0_local_values.size());

    summary.child1_local_summary.delta_t_child_set = child1_delta_t_child_set;
    summary.child1_local_summary.num_sources = child1_num_sources;
    summary.child1_local_summary.mean = computeMeanSeries(child1_local_values);
    summary.child1_local_summary.stddev =
      computeStddevSeries(child1_local_values, summary.child1_local_summary.mean);
    summary.child1_local_summary.stderr =
      computeSampleStderr(summary.child1_local_summary.stddev,
                          child1_local_values.size());

    summary.pass = true;
    for (int i = 0; i < summary.delta_t_parent_set.size(); ++i)
    {
      const double mu = toDouble(summary.mean_delta_correlator[i]);
      const double se = toDouble(summary.stderr_delta_correlator[i]);
      if (se == 0.0)
      {
        if (mu != 0.0)
          summary.pass = false;
      }
      else if (std::fabs(mu) > 3.0 * se)
      {
        summary.pass = false;
      }
    }

    writeOuterEnsembleSummary(input, summary);
    return summary;
  }

  void writeOuterEnsembleCsv(const CheckerInput& input,
                             const OuterEnsembleSummary& summary)
  {
    if (input.csv_file.empty())
      return;

    std::ofstream out(input.csv_file.c_str());
    check(out.good(), "failed to open CSV output file");
    out << "series,index,value,stderr\n";
    out << std::setprecision(17);

    for (int i = 0; i < summary.delta_t_parent_set.size(); ++i)
    {
      out << "two_level," << summary.delta_t_parent_set[i]
          << "," << toDouble(summary.mean_two_level_correlator[i])
          << "," << toDouble(summary.stderr_two_level_correlator[i]) << "\n";
      out << "parent_window," << summary.delta_t_parent_set[i]
          << "," << toDouble(summary.mean_parent_window_correlator[i])
          << "," << toDouble(summary.stderr_parent_window_correlator[i]) << "\n";
      out << "delta," << summary.delta_t_parent_set[i]
          << "," << toDouble(summary.mean_delta_correlator[i])
          << "," << toDouble(summary.stderr_delta_correlator[i]) << "\n";
    }

    for (int i = 0; i < summary.child0_local_summary.delta_t_child_set.size(); ++i)
    {
      out << "child0_local," << summary.child0_local_summary.delta_t_child_set[i]
          << "," << toDouble(summary.child0_local_summary.mean[i])
          << "," << toDouble(summary.child0_local_summary.stderr[i]) << "\n";
    }

    for (int i = 0; i < summary.child1_local_summary.delta_t_child_set.size(); ++i)
    {
      out << "child1_local," << summary.child1_local_summary.delta_t_child_set[i]
          << "," << toDouble(summary.child1_local_summary.mean[i])
          << "," << toDouble(summary.child1_local_summary.stderr[i]) << "\n";
    }
  }

  void runOuterEnsembleMode(const CheckerInput& input)
  {
    const OuterEnsembleSummary summary = buildOuterEnsembleSummary(input);
    writeOuterEnsembleCsv(input, summary);
    check(summary.pass,
          "outer-ensemble two-level estimator is not statistically consistent "
          "with the parent-window comparator");

    QDPIO::cout << "t_qactden_0mp_corr: passed" << std::endl;
    QDPIO::cout << "  mode: " << input.mode << std::endl;
    QDPIO::cout << "  outer_samples: " << summary.n_outer_samples << std::endl;
    QDPIO::cout << "  pass: " << summary.pass << std::endl;
  }
}

int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  int status = 0;

  try
  {
    XMLReader xml_in(Chroma::getXMLInputFileName());
    CheckerInput input;
    readCheckerInput(xml_in, "/qactden_0mp_corr_check", input);
    validateInput(input);

    if (input.mode == MODE_PERIODIC_ALL_T)
      runPeriodicMode(input);
    else if (input.mode == MODE_INTERIOR_NOWRAP_CHILD)
      runChildSummaryMode(input);
    else if (input.mode == MODE_PARENT_WINDOW)
      runParentWindowMode(input);
    else if (input.mode == MODE_TWO_LEVEL_CROSS_DOMAIN)
      runOuterSampleMode(input);
    else if (input.mode == MODE_TWO_LEVEL_OUTER_ENSEMBLE)
      runOuterEnsembleMode(input);
    else
      fail("unsupported reducer mode in main");

    if (!input.summary_file.empty())
      QDPIO::cout << "  summary_xml: " << input.summary_file << std::endl;
    if (!input.csv_file.empty())
      QDPIO::cout << "  csv: " << input.csv_file << std::endl;
  }
  catch (const std::string& e)
  {
    QDPIO::cerr << e << std::endl;
    status = 1;
  }
  catch (std::exception& e)
  {
    QDPIO::cerr << "t_qactden_0mp_corr: standard exception: "
                << e.what() << std::endl;
    status = 1;
  }
  catch (...)
  {
    QDPIO::cerr << "t_qactden_0mp_corr: unknown exception" << std::endl;
    status = 1;
  }

  Chroma::finalize();
  return status;
}
