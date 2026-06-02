#include "chroma.h"
#include "util/gauge/gauge_subdomain_split.h"
#include "qdp_util.h"

#include <sstream>

using namespace Chroma;

namespace
{
  struct HostGaugeField
  {
    multi1d<int> nrow;
    multi1d< multi1d<ColorMatrix> > link;
  };

  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_gauge_subdomain_split: " << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  void activateLayout(const multi1d<int>& nrow)
  {
    bool same = (Layout::lattSize().size() == nrow.size());
    if (same)
    {
      for (int i = 0; i < nrow.size(); ++i)
      {
        if (Layout::lattSize()[i] != nrow[i])
        {
          same = false;
          break;
        }
      }
    }

    if (same)
      return;

    Layout::destroy();
    Layout::setLattSize(nrow);
    Layout::create();
  }

  int latticeVolume(const multi1d<int>& nrow)
  {
    int vol = 1;
    for (int i = 0; i < nrow.size(); ++i)
      vol *= nrow[i];
    return vol;
  }

  bool exactMatch(const ColorMatrix& lhs, const ColorMatrix& rhs)
  {
    return toDouble(norm2(lhs - rhs)) == 0.0;
  }

  HostGaugeField extractHost(const multi1d<LatticeColorMatrix>& u, const multi1d<int>& nrow)
  {
    HostGaugeField host;
    host.nrow = nrow;
    host.link.resize(Nd);

    const int vol = latticeVolume(nrow);
    for (int mu = 0; mu < Nd; ++mu)
      host.link[mu].resize(vol);

    for (int site = 0; site < vol; ++site)
    {
      multi1d<int> coord = crtesn(site, nrow);
      for (int mu = 0; mu < Nd; ++mu)
        host.link[mu][site] = peekSite(u[mu], coord);
    }

    return host;
  }

  multi1d<int> buildGlobalToLocal(const multi1d<int>& local_to_global_t, int Lt)
  {
    multi1d<int> global_to_local(Lt);
    global_to_local = -1;

    for (int local_t = 0; local_t < local_to_global_t.size(); ++local_t)
    {
      const int global_t = local_to_global_t[local_t];
      check(global_to_local[global_t] == -1,
            "local_to_global_t should not repeat global times");
      global_to_local[global_t] = local_t;
    }

    return global_to_local;
  }

  void checkTimeOrder(const multi1d<int>& got,
                      const int* expected,
                      int expected_size,
                      const std::string& label)
  {
    check(got.size() == expected_size,
          label + ": unexpected local_to_global_t length");

    for (int i = 0; i < expected_size; ++i)
    {
      if (got[i] != expected[i])
      {
        std::ostringstream os;
        os << label << ": expected local_to_global_t[" << i << "]="
           << expected[i] << " but found " << got[i];
        fail(os.str());
      }
    }
  }

  void checkFrozenIntervals(const GaugeSubdomainSplitChild& child,
                            int first_start,
                            int first_end,
                            int second_start,
                            int second_end,
                            const std::string& label)
  {
    check(child.frozen_local_intervals.size() == 2,
          label + ": expected exactly two frozen intervals");

    check(child.frozen_local_intervals[0].t_start == first_start &&
          child.frozen_local_intervals[0].t_end == first_end,
          label + ": unexpected first frozen interval");

    check(child.frozen_local_intervals[1].t_start == second_start &&
          child.frozen_local_intervals[1].t_end == second_end,
          label + ": unexpected second frozen interval");
  }

  void modifyActiveSlices(GaugeSubdomainSplitChild& child, int t_dir)
  {
    activateLayout(child.nrow);

    const int first_frozen_end = child.frozen_local_intervals[0].t_end;
    const int second_frozen_start = child.frozen_local_intervals[1].t_start;
    const int active_start = first_frozen_end + 1;
    const int active_end = second_frozen_start - 1;

    check(active_start <= active_end,
          "expected a nonempty active interval in child lattice");

    LatticeInteger t = Layout::latticeCoordinate(t_dir);
    LatticeBoolean active_mask = (t >= active_start) && (t <= active_end);

    for (int mu = 0; mu < child.u.size(); ++mu)
    {
      LatticeColorMatrix replacement;
      gaussian(replacement);
      copymask(child.u[mu], active_mask, replacement);
    }
  }

  void checkExactFieldEquality(const multi1d<LatticeColorMatrix>& lhs,
                               const multi1d<LatticeColorMatrix>& rhs,
                               const std::string& label)
  {
    check(lhs.size() == rhs.size(), label + ": direction count mismatch");

    Double total_diff = Double(0);
    for (int mu = 0; mu < lhs.size(); ++mu)
      total_diff += norm2(lhs[mu] - rhs[mu]);

    if (toDouble(total_diff) != 0.0)
    {
      std::ostringstream os;
      os << label << ": fields differ with total norm2 " << toDouble(total_diff);
      fail(os.str());
    }
  }

  void verifyOwnershipAfterStitch(const GaugeSubdomainSplitResult& split,
                                  const multi1d<LatticeColorMatrix>& stitched,
                                  const std::string& label)
  {
    activateLayout(split.parent_nrow);
    HostGaugeField parent_host = extractHost(stitched, split.parent_nrow);

    activateLayout(split.child0.nrow);
    HostGaugeField child0_host = extractHost(split.child0.u, split.child0.nrow);

    activateLayout(split.child1.nrow);
    HostGaugeField child1_host = extractHost(split.child1.u, split.child1.nrow);

    const int Lt = split.parent_nrow[split.param.t_dir];
    const multi1d<int> child0_global_to_local =
      buildGlobalToLocal(split.child0.local_to_global_t, Lt);
    const multi1d<int> child1_global_to_local =
      buildGlobalToLocal(split.child1.local_to_global_t, Lt);

    const int parent_vol = latticeVolume(split.parent_nrow);
    for (int site = 0; site < parent_vol; ++site)
    {
      multi1d<int> parent_coord = crtesn(site, split.parent_nrow);
      const int global_t = parent_coord[split.param.t_dir];
      const int child0_local_t = child0_global_to_local[global_t];
      const int child1_local_t = child1_global_to_local[global_t];

      check(child0_local_t >= 0 || child1_local_t >= 0,
            label + ": parent slice missing from both children");

      multi1d<int> child0_coord = parent_coord;
      multi1d<int> child1_coord = parent_coord;

      if (child0_local_t >= 0)
        child0_coord[split.param.t_dir] = child0_local_t;
      if (child1_local_t >= 0)
        child1_coord[split.param.t_dir] = child1_local_t;

      const int child0_site =
        (child0_local_t >= 0) ? local_site(child0_coord, split.child0.nrow) : -1;
      const int child1_site =
        (child1_local_t >= 0) ? local_site(child1_coord, split.child1.nrow) : -1;

      for (int mu = 0; mu < Nd; ++mu)
      {
        const ColorMatrix& stitched_link = parent_host.link[mu][site];

        if (child0_site >= 0 && child1_site >= 0)
        {
          check(exactMatch(child0_host.link[mu][child0_site],
                           child1_host.link[mu][child1_site]),
                label + ": frozen boundary copies should match exactly");

          check(exactMatch(stitched_link, child0_host.link[mu][child0_site]),
                label + ": stitched frozen link did not come from the shared boundary");
        }
        else if (child0_site >= 0)
        {
          check(exactMatch(stitched_link, child0_host.link[mu][child0_site]),
                label + ": active child0 link was not stitched correctly");
        }
        else
        {
          check(exactMatch(stitched_link, child1_host.link[mu][child1_site]),
                label + ": active child1 link was not stitched correctly");
        }
      }
    }

    activateLayout(split.parent_nrow);
  }

  void expectSplitThrows(const multi1d<LatticeColorMatrix>& parent_u,
                         const GaugeSubdomainSplitParams& param,
                         const std::string& needle,
                         const std::string& label)
  {
    try
    {
      (void)splitGaugeSubdomains(parent_u, param);
      fail(label + ": expected splitGaugeSubdomains() to throw");
    }
    catch (const std::string& e)
    {
      if (e.find(needle) == std::string::npos)
      {
        std::ostringstream os;
        os << label << ": expected error containing \"" << needle
           << "\" but received \"" << e << "\"";
        fail(os.str());
      }
    }
  }

  void expectStitchThrows(const GaugeSubdomainSplitResult& split,
                          const std::string& needle,
                          const std::string& label)
  {
    try
    {
      (void)stitchGaugeSubdomains(split);
      fail(label + ": expected stitchGaugeSubdomains() to throw");
    }
    catch (const std::string& e)
    {
      if (e.find(needle) == std::string::npos)
      {
        std::ostringstream os;
        os << label << ": expected error containing \"" << needle
           << "\" but received \"" << e << "\"";
        fail(os.str());
      }
    }
  }
}

int main(int argc, char* argv[])
{
  Chroma::initialize(&argc, &argv);

  try
  {
    {
      const int nrow_arr[] = {2, 2, 2, 8};
      multi1d<int> nrow(Nd);
      nrow = nrow_arr;
      Layout::setLattSize(nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> parent_u(Nd);
      for (int mu = 0; mu < Nd; ++mu)
        gaussian(parent_u[mu]);

      multi1d<LatticeColorMatrix> parent_before = parent_u;

      GaugeSubdomainSplitParams param;
      param.t_dir = Nd - 1;
      param.cut0 = 1;
      param.cut1 = 5;
      param.frozen_width = 1;

      GaugeSubdomainSplitResult split = splitGaugeSubdomains(parent_u, param);

      const int child0_expected[] = {1, 2, 3, 4, 5};
      const int child1_expected[] = {5, 6, 7, 0, 1};

      check(split.param.cut0 == 1 && split.param.cut1 == 5,
            "NonWrappingSplit: cuts should be stored in increasing order");
      checkTimeOrder(split.child0.local_to_global_t, child0_expected, 5, "NonWrappingSplit child0");
      checkTimeOrder(split.child1.local_to_global_t, child1_expected, 5, "NonWrappingSplit child1");
      checkFrozenIntervals(split.child0, 0, 0, 4, 4, "NonWrappingSplit child0");
      checkFrozenIntervals(split.child1, 0, 0, 4, 4, "NonWrappingSplit child1");
      check(split.child0.nrow[split.param.t_dir] == 5,
            "NonWrappingSplit: child0 extent mismatch");
      check(split.child1.nrow[split.param.t_dir] == 5,
            "NonWrappingSplit: child1 extent mismatch");

      activateLayout(split.parent_nrow);
      multi1d<LatticeColorMatrix> stitched = stitchGaugeSubdomains(split);
      checkExactFieldEquality(stitched, parent_before, "ImmediateRoundTrip");

      modifyActiveSlices(split.child0, split.param.t_dir);
      modifyActiveSlices(split.child1, split.param.t_dir);

      activateLayout(split.parent_nrow);
      stitched = stitchGaugeSubdomains(split);
      verifyOwnershipAfterStitch(split, stitched, "InteriorOnlyEvolutionSurrogate");

      activateLayout(split.child0.nrow);
      ColorMatrix z = zero;
      multi1d<int> coord(Nd);
      coord = 0;
      pokeSite(split.child0.u[0], z, coord);
      expectStitchThrows(split, "duplicated frozen-boundary link mismatch",
                         "FrozenMismatchDetection");
    }

    {
      const int nrow_arr[] = {2, 2, 2, 16};
      multi1d<int> nrow(Nd);
      nrow = nrow_arr;
      activateLayout(nrow);

      multi1d<LatticeColorMatrix> parent_u(Nd);
      for (int mu = 0; mu < Nd; ++mu)
        gaussian(parent_u[mu]);

      GaugeSubdomainSplitParams param;
      param.t_dir = Nd - 1;
      param.cut0 = 4;
      param.cut1 = 14;
      param.frozen_width = 3;

      GaugeSubdomainSplitResult split = splitGaugeSubdomains(parent_u, param);

      const int child0_expected[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0};
      const int child1_expected[] = {14, 15, 0, 1, 2, 3, 4, 5, 6};

      checkTimeOrder(split.child0.local_to_global_t, child0_expected, 13, "WrappingFrozenBlock child0");
      checkTimeOrder(split.child1.local_to_global_t, child1_expected, 9, "WrappingFrozenBlock child1");
      checkFrozenIntervals(split.child0, 0, 2, 10, 12, "WrappingFrozenBlock child0");
      checkFrozenIntervals(split.child1, 0, 2, 6, 8, "WrappingFrozenBlock child1");
    }

    {
      const int nrow_arr[] = {2, 2, 2, 8};
      multi1d<int> nrow(Nd);
      nrow = nrow_arr;
      activateLayout(nrow);

      multi1d<LatticeColorMatrix> parent_u(Nd);
      for (int mu = 0; mu < Nd; ++mu)
        gaussian(parent_u[mu]);

      GaugeSubdomainSplitParams param;
      param.t_dir = Nd;
      param.cut0 = 1;
      param.cut1 = 5;
      param.frozen_width = 1;
      expectSplitThrows(parent_u, param, "t_dir=", "InvalidTDir");

      param.t_dir = Nd - 1;
      param.cut0 = 2;
      param.cut1 = 2;
      expectSplitThrows(parent_u, param, "must be different", "EqualCuts");

      param.cut0 = 1;
      param.cut1 = 5;
      param.frozen_width = 0;
      expectSplitThrows(parent_u, param, "must be positive", "NonPositiveWidth");

      param.cut0 = 1;
      param.cut1 = 2;
      param.frozen_width = 2;
      expectSplitThrows(parent_u, param, "overlap", "OverlappingFrozenBlocks");

      param.cut0 = 1;
      param.cut1 = 3;
      param.frozen_width = 2;
      expectSplitThrows(parent_u, param, "do not leave two nonempty independent subdomains",
                        "TouchingFrozenBlocks");
    }
  }
  catch (const std::string& e)
  {
    fail(e);
  }

  Chroma::finalize();
  return 0;
}
