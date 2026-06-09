/*! \file
 *  \brief Gauge subdomain split/stitch utilities
 */

#include "util/gauge/gauge_subdomain_split.h"
#include "qdp_util.h"

#include <sstream>
#include <vector>

namespace Chroma
{

  namespace
  {
    void failGaugeSubdomainSplit(const std::string& message)
    {
      throw std::string("GaugeSubdomainSplit: " + message);
    }

    bool sameShape(const multi1d<int>& a, const multi1d<int>& b)
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

    int latticeVolume(const multi1d<int>& nrow)
    {
      if (nrow.size() != Nd)
        failGaugeSubdomainSplit("lattice extent must have length Nd");

      int vol = 1;
      for (int mu = 0; mu < Nd; ++mu)
      {
        if (nrow[mu] <= 0)
        {
          std::ostringstream os;
          os << "lattice extent nrow[" << mu << "]=" << nrow[mu]
             << " must be positive";
          failGaugeSubdomainSplit(os.str());
        }

        vol *= nrow[mu];
      }

      return vol;
    }

    int advanceTime(int t, int n, int Lt)
    {
      int v = (t + n) % Lt;
      if (v < 0)
        v += Lt;
      return v;
    }

    void validateParams(const multi1d<int>& parent_nrow,
                        const GaugeSubdomainSplitParams& param,
                        int& Lt,
                        int& cut_left,
                        int& cut_right,
                        int& lenA,
                        int& lenB)
    {
      latticeVolume(parent_nrow);

      if (param.t_dir < 0 || param.t_dir >= Nd)
      {
        std::ostringstream os;
        os << "t_dir=" << param.t_dir << " is outside [0," << Nd - 1 << "]";
        failGaugeSubdomainSplit(os.str());
      }

      Lt = parent_nrow[param.t_dir];

      if (param.cut0 < 0 || param.cut0 >= Lt)
      {
        std::ostringstream os;
        os << "cut0=" << param.cut0 << " is outside [0," << Lt - 1 << "]";
        failGaugeSubdomainSplit(os.str());
      }

      if (param.cut1 < 0 || param.cut1 >= Lt)
      {
        std::ostringstream os;
        os << "cut1=" << param.cut1 << " is outside [0," << Lt - 1 << "]";
        failGaugeSubdomainSplit(os.str());
      }

      if (param.cut0 == param.cut1)
        failGaugeSubdomainSplit("cut0 and cut1 must differ");

      if (param.frozen_width <= 0)
      {
        std::ostringstream os;
        os << "frozen_width=" << param.frozen_width << " must be positive";
        failGaugeSubdomainSplit(os.str());
      }

      cut_left = std::min(param.cut0, param.cut1);
      cut_right = std::max(param.cut0, param.cut1);

      std::vector<int> frozen_owner(Lt, 0);
      for (int k = 0; k < param.frozen_width; ++k)
      {
        const int t0 = advanceTime(cut_left, k, Lt);
        const int t1 = advanceTime(cut_right, k, Lt);

        if (frozen_owner[t0] != 0)
          failGaugeSubdomainSplit("left frozen block overlaps itself");
        frozen_owner[t0] = 1;

        if (frozen_owner[t1] != 0)
          failGaugeSubdomainSplit("frozen boundary blocks overlap");
        frozen_owner[t1] = 2;
      }

      lenA = advanceTime(cut_right, -cut_left - param.frozen_width, Lt);
      lenB = advanceTime(cut_left, -cut_right - param.frozen_width, Lt);

      if (lenA <= 0 || lenB <= 0)
        failGaugeSubdomainSplit("frozen blocks must leave two nonempty active subdomains");
    }

    multi1d<int> buildOrdering(int first_cut,
                               int second_cut,
                               int frozen_width,
                               int active_len,
                               int Lt)
    {
      multi1d<int> ordering(2 * frozen_width + active_len);
      int p = 0;

      for (int k = 0; k < frozen_width; ++k)
        ordering[p++] = advanceTime(first_cut, k, Lt);

      for (int k = 0; k < active_len; ++k)
        ordering[p++] = advanceTime(first_cut + frozen_width, k, Lt);

      for (int k = 0; k < frozen_width; ++k)
        ordering[p++] = advanceTime(second_cut, k, Lt);

      return ordering;
    }

    multi1d<GaugeSubdomainSplitInterval>
    buildFrozenLocalIntervals(const multi1d<int>& child_nrow, int t_dir, int frozen_width)
    {
      multi1d<GaugeSubdomainSplitInterval> intervals(2);

      intervals[0].t_start = 0;
      intervals[0].t_end = frozen_width - 1;

      intervals[1].t_start = child_nrow[t_dir] - frozen_width;
      intervals[1].t_end = child_nrow[t_dir] - 1;

      return intervals;
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

    void validateIntervals(const multi1d<GaugeSubdomainSplitInterval>& intervals,
                           int child_extent,
                           const std::string& label)
    {
      if (intervals.size() != 2)
      {
        std::ostringstream os;
        os << label << " must contain exactly two frozen local intervals";
        failGaugeSubdomainSplit(os.str());
      }

      for (int i = 0; i < intervals.size(); ++i)
      {
        if (intervals[i].t_start < 0 ||
            intervals[i].t_end < intervals[i].t_start ||
            intervals[i].t_end >= child_extent)
        {
          std::ostringstream os;
          os << label << " interval " << i << " is invalid for extent "
             << child_extent;
          failGaugeSubdomainSplit(os.str());
        }
      }
    }

    void validatePlan(const GaugeSubdomainSplitPlan& plan)
    {
      GaugeSubdomainSplitPlan derived =
        makeGaugeSubdomainSplitPlan(plan.parent_nrow, plan.param);

      if (!sameShape(plan.child0_nrow, derived.child0_nrow) ||
          !sameShape(plan.child1_nrow, derived.child1_nrow) ||
          !sameShape(plan.child0_local_to_global_t, derived.child0_local_to_global_t) ||
          !sameShape(plan.child1_local_to_global_t, derived.child1_local_to_global_t))
      {
        failGaugeSubdomainSplit("split plan metadata is inconsistent with the stored parameters");
      }

      if (plan.child0_nrow.size() != Nd || plan.child1_nrow.size() != Nd)
        failGaugeSubdomainSplit("split plan child lattice extents must have length Nd");

      if (plan.child0_frozen_local_intervals.size() != derived.child0_frozen_local_intervals.size() ||
          plan.child1_frozen_local_intervals.size() != derived.child1_frozen_local_intervals.size())
      {
        failGaugeSubdomainSplit("split plan frozen local interval count is inconsistent");
      }

      for (int i = 0; i < plan.child0_frozen_local_intervals.size(); ++i)
      {
        if (plan.child0_frozen_local_intervals[i].t_start != derived.child0_frozen_local_intervals[i].t_start ||
            plan.child0_frozen_local_intervals[i].t_end != derived.child0_frozen_local_intervals[i].t_end)
        {
          failGaugeSubdomainSplit("child0 frozen local intervals are inconsistent with the split plan");
        }
      }

      for (int i = 0; i < plan.child1_frozen_local_intervals.size(); ++i)
      {
        if (plan.child1_frozen_local_intervals[i].t_start != derived.child1_frozen_local_intervals[i].t_start ||
            plan.child1_frozen_local_intervals[i].t_end != derived.child1_frozen_local_intervals[i].t_end)
        {
          failGaugeSubdomainSplit("child1 frozen local intervals are inconsistent with the split plan");
        }
      }

      validateIntervals(plan.child0_frozen_local_intervals,
                        plan.child0_nrow[plan.param.t_dir],
                        "child0_frozen_local_intervals");
      validateIntervals(plan.child1_frozen_local_intervals,
                        plan.child1_nrow[plan.param.t_dir],
                        "child1_frozen_local_intervals");
    }

    int linkIndex(int mu, int site, int vol)
    {
      return mu * vol + site;
    }
  }

  GaugeSubdomainSplitParams::GaugeSubdomainSplitParams()
    : t_dir(Nd - 1), cut0(0), cut1(0), frozen_width(1)
  {
  }

  GaugeSubdomainSplitParams::GaugeSubdomainSplitParams(XMLReader& xml,
                                                       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "t_dir", t_dir);
    read(paramtop, "cut0", cut0);
    read(paramtop, "cut1", cut1);
    read(paramtop, "frozen_width", frozen_width);
  }

  GaugeSubdomainSplitInterval::GaugeSubdomainSplitInterval()
    : t_start(0), t_end(-1)
  {
  }

  GaugeSubdomainSplitPlan::GaugeSubdomainSplitPlan()
  {
  }

  GaugeSubdomainSplitPlan::GaugeSubdomainSplitPlan(XMLReader& xml,
                                                   const std::string& path)
  {
    XMLReader plan_xml(xml, path);

    read(plan_xml, "Param", param);

    XMLReader parent_xml(plan_xml, "Parent");
    read(parent_xml, "nrow", parent_nrow);

    XMLReader child0_xml(plan_xml, "Child0");
    read(child0_xml, "nrow", child0_nrow);
    read(child0_xml, "local_to_global_t", child0_local_to_global_t);
    read(child0_xml, "frozen_local_intervals", child0_frozen_local_intervals);

    XMLReader child1_xml(plan_xml, "Child1");
    read(child1_xml, "nrow", child1_nrow);
    read(child1_xml, "local_to_global_t", child1_local_to_global_t);
    read(child1_xml, "frozen_local_intervals", child1_frozen_local_intervals);

    validatePlan(*this);
  }

  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitParams& p)
  {
    GaugeSubdomainSplitParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const GaugeSubdomainSplitParams& p)
  {
    push(xml, path);
    write(xml, "t_dir", p.t_dir);
    write(xml, "cut0", p.cut0);
    write(xml, "cut1", p.cut1);
    write(xml, "frozen_width", p.frozen_width);
    pop(xml);
  }

  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitInterval& p)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "t_start", p.t_start);
    read(paramtop, "t_end", p.t_end);
  }

  void write(XMLWriter& xml, const std::string& path, const GaugeSubdomainSplitInterval& p)
  {
    push(xml, path);
    write(xml, "t_start", p.t_start);
    write(xml, "t_end", p.t_end);
    pop(xml);
  }

  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitPlan& p)
  {
    GaugeSubdomainSplitPlan tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const GaugeSubdomainSplitPlan& p)
  {
    validatePlan(p);

    push(xml, path);
    write(xml, "Param", p.param);

    push(xml, "Parent");
    write(xml, "nrow", p.parent_nrow);
    pop(xml);

    push(xml, "Child0");
    write(xml, "nrow", p.child0_nrow);
    write(xml, "local_to_global_t", p.child0_local_to_global_t);
    write(xml, "frozen_local_intervals", p.child0_frozen_local_intervals);
    pop(xml);

    push(xml, "Child1");
    write(xml, "nrow", p.child1_nrow);
    write(xml, "local_to_global_t", p.child1_local_to_global_t);
    write(xml, "frozen_local_intervals", p.child1_frozen_local_intervals);
    pop(xml);

    pop(xml);
  }

  GaugeSubdomainSplitPlan
  makeGaugeSubdomainSplitPlan(const multi1d<int>& parent_nrow,
                              const GaugeSubdomainSplitParams& param)
  {
    GaugeSubdomainSplitPlan plan;
    plan.param = param;
    plan.parent_nrow = parent_nrow;

    int Lt = 0;
    int cut_left = 0;
    int cut_right = 0;
    int lenA = 0;
    int lenB = 0;
    validateParams(parent_nrow, param, Lt, cut_left, cut_right, lenA, lenB);

    plan.child0_nrow = parent_nrow;
    plan.child1_nrow = parent_nrow;
    plan.child0_nrow[param.t_dir] = 2 * param.frozen_width + lenA;
    plan.child1_nrow[param.t_dir] = 2 * param.frozen_width + lenB;

    plan.child0_local_to_global_t =
      buildOrdering(cut_left, cut_right, param.frozen_width, lenA, Lt);
    plan.child1_local_to_global_t =
      buildOrdering(cut_right, cut_left, param.frozen_width, lenB, Lt);

    plan.child0_frozen_local_intervals =
      buildFrozenLocalIntervals(plan.child0_nrow, param.t_dir, param.frozen_width);
    plan.child1_frozen_local_intervals =
      buildFrozenLocalIntervals(plan.child1_nrow, param.t_dir, param.frozen_width);

    return plan;
  }

  void snapshotGaugeField(const multi1d<LatticeColorMatrix>& u,
                          multi1d<ColorMatrix>& site_links)
  {
    if (u.size() != Nd)
      failGaugeSubdomainSplit("snapshotGaugeField expects Nd gauge links");

    const multi1d<int>& nrow = Layout::lattSize();
    const int vol = latticeVolume(nrow);
    site_links.resize(Nd * vol);

    for (int site = 0; site < vol; ++site)
    {
      multi1d<int> coord = crtesn(site, nrow);
      for (int mu = 0; mu < Nd; ++mu)
        site_links[linkIndex(mu, site, vol)] = peekSite(u[mu], coord);
    }
  }

  void materializeGaugeField(const multi1d<int>& nrow,
                             const multi1d<ColorMatrix>& site_links,
                             multi1d<LatticeColorMatrix>& u)
  {
    if (!sameShape(Layout::lattSize(), nrow))
      failGaugeSubdomainSplit("current layout does not match the requested gauge-field materialization extent");

    const int vol = latticeVolume(nrow);
    if (site_links.size() != Nd * vol)
      failGaugeSubdomainSplit("materializeGaugeField received a buffer with the wrong size");

    u.resize(Nd);
    for (int mu = 0; mu < Nd; ++mu)
      u[mu] = zero;

    for (int site = 0; site < vol; ++site)
    {
      multi1d<int> coord = crtesn(site, nrow);
      for (int mu = 0; mu < Nd; ++mu)
        pokeSite(u[mu], site_links[linkIndex(mu, site, vol)], coord);
    }
  }

  void extractGaugeSubdomain(const multi1d<ColorMatrix>& parent_links,
                             const GaugeSubdomainSplitPlan& plan,
                             int child_id,
                             multi1d<ColorMatrix>& child_links)
  {
    validatePlan(plan);

    if (child_id != 0 && child_id != 1)
      failGaugeSubdomainSplit("extractGaugeSubdomain child_id must be 0 or 1");

    const multi1d<int>& child_nrow = (child_id == 0) ? plan.child0_nrow : plan.child1_nrow;
    const multi1d<int>& local_to_global_t =
      (child_id == 0) ? plan.child0_local_to_global_t : plan.child1_local_to_global_t;

    const int parent_vol = latticeVolume(plan.parent_nrow);
    const int child_vol = latticeVolume(child_nrow);

    if (parent_links.size() != Nd * parent_vol)
      failGaugeSubdomainSplit("extractGaugeSubdomain parent buffer has the wrong size");

    if (local_to_global_t.size() != child_nrow[plan.param.t_dir])
      failGaugeSubdomainSplit("extractGaugeSubdomain local_to_global_t length does not match child extent");

    child_links.resize(Nd * child_vol);

    for (int child_site = 0; child_site < child_vol; ++child_site)
    {
      multi1d<int> child_coord = crtesn(child_site, child_nrow);
      multi1d<int> parent_coord = child_coord;
      parent_coord[plan.param.t_dir] = local_to_global_t[child_coord[plan.param.t_dir]];

      const int parent_site = local_site(parent_coord, plan.parent_nrow);
      for (int mu = 0; mu < Nd; ++mu)
      {
        child_links[linkIndex(mu, child_site, child_vol)] =
          parent_links[linkIndex(mu, parent_site, parent_vol)];
      }
    }
  }

  void stitchGaugeSubdomains(const GaugeSubdomainSplitPlan& plan,
                             const multi1d<ColorMatrix>& child0_links,
                             const multi1d<ColorMatrix>& child1_links,
                             multi1d<ColorMatrix>& parent_links)
  {
    validatePlan(plan);

    const int parent_vol = latticeVolume(plan.parent_nrow);
    const int child0_vol = latticeVolume(plan.child0_nrow);
    const int child1_vol = latticeVolume(plan.child1_nrow);

    if (child0_links.size() != Nd * child0_vol)
      failGaugeSubdomainSplit("child0 buffer has the wrong size for stitchGaugeSubdomains");
    if (child1_links.size() != Nd * child1_vol)
      failGaugeSubdomainSplit("child1 buffer has the wrong size for stitchGaugeSubdomains");

    parent_links.resize(Nd * parent_vol);
    std::vector<int> parent_written(parent_vol, 0);

    for (int child_site = 0; child_site < child0_vol; ++child_site)
    {
      multi1d<int> child_coord = crtesn(child_site, plan.child0_nrow);
      multi1d<int> parent_coord = child_coord;
      parent_coord[plan.param.t_dir] =
        plan.child0_local_to_global_t[child_coord[plan.param.t_dir]];
      const int parent_site = local_site(parent_coord, plan.parent_nrow);

      for (int mu = 0; mu < Nd; ++mu)
      {
        parent_links[linkIndex(mu, parent_site, parent_vol)] =
          child0_links[linkIndex(mu, child_site, child0_vol)];
      }

      parent_written[parent_site] = 1;
    }

    for (int child_site = 0; child_site < child1_vol; ++child_site)
    {
      multi1d<int> child_coord = crtesn(child_site, plan.child1_nrow);
      multi1d<int> parent_coord = child_coord;
      parent_coord[plan.param.t_dir] =
        plan.child1_local_to_global_t[child_coord[plan.param.t_dir]];
      const int parent_site = local_site(parent_coord, plan.parent_nrow);
      const bool frozen =
        isFrozenLocalTime(child_coord[plan.param.t_dir], plan.child1_frozen_local_intervals);

      if (frozen)
      {
        if (!parent_written[parent_site])
          failGaugeSubdomainSplit("child1 frozen slice has no matching child0 source during stitch");

        for (int mu = 0; mu < Nd; ++mu)
        {
          const ColorMatrix& lhs = parent_links[linkIndex(mu, parent_site, parent_vol)];
          const ColorMatrix& rhs = child1_links[linkIndex(mu, child_site, child1_vol)];
          if (toDouble(norm2(lhs - rhs)) != 0.0)
          {
            std::ostringstream os;
            os << "duplicated frozen boundary mismatch at global t="
               << parent_coord[plan.param.t_dir] << " direction " << mu;
            failGaugeSubdomainSplit(os.str());
          }
        }
      }
      else
      {
        if (parent_written[parent_site])
          failGaugeSubdomainSplit("child active regions overlap unexpectedly during stitch");

        for (int mu = 0; mu < Nd; ++mu)
        {
          parent_links[linkIndex(mu, parent_site, parent_vol)] =
            child1_links[linkIndex(mu, child_site, child1_vol)];
        }

        parent_written[parent_site] = 1;
      }
    }

    for (int parent_site = 0; parent_site < parent_vol; ++parent_site)
    {
      if (!parent_written[parent_site])
      {
        std::ostringstream os;
        os << "stitched parent site " << parent_site
           << " was not covered by either child";
        failGaugeSubdomainSplit(os.str());
      }
    }
  }

}
