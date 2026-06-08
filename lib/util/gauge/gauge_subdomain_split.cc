/*! \file
 * \brief Split and stitch utilities for temporal gauge subdomains
 */

#include "util/gauge/gauge_subdomain_split.h"
#include "qdp_util.h"

#include <algorithm>
#include <sstream>

namespace Chroma
{

  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitParams& p)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "t_dir", p.t_dir);
    read(paramtop, "cut0", p.cut0);
    read(paramtop, "cut1", p.cut1);
    read(paramtop, "frozen_width", p.frozen_width);
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

  namespace
  {
    struct HostGaugeField
    {
      multi1d<int> nrow;
      multi1d< multi1d<ColorMatrix> > link;
    };

    struct SplitPlan
    {
      GaugeSubdomainSplitParams param;
      multi1d<int> parent_nrow;
      int Lt;
      int lenA;
      int lenB;
      multi1d<int> F0;
      multi1d<int> F1;
      multi1d<int> A;
      multi1d<int> B;
      multi1d<int> child0_order;
      multi1d<int> child1_order;
      multi1d<int> child0_nrow;
      multi1d<int> child1_nrow;
    };

    class LayoutGuard
    {
    public:
      LayoutGuard() : saved_nrow(Layout::lattSize()) {}

      ~LayoutGuard()
      {
        activate(saved_nrow);
      }

      static void activate(const multi1d<int>& nrow)
      {
        if (Layout::lattSize().size() == nrow.size())
        {
          bool same = true;
          for (int i = 0; i < nrow.size(); ++i)
          {
            if (Layout::lattSize()[i] != nrow[i])
            {
              same = false;
              break;
            }
          }

          if (same)
            return;
        }

        Layout::destroy();
        Layout::setLattSize(nrow);
        Layout::create();
      }

      void setRestore(const multi1d<int>& nrow)
      {
        saved_nrow = nrow;
      }

    private:
      multi1d<int> saved_nrow;
    };

    void throwGaugeSubdomainSplit(const std::string& message)
    {
      throw std::string("GaugeSubdomainSplit: " + message);
    }

    int latticeVolume(const multi1d<int>& nrow)
    {
      int vol = 1;
      for (int i = 0; i < nrow.size(); ++i)
        vol *= nrow[i];
      return vol;
    }

    int advanceTime(int t, int amount, int Lt)
    {
      return (t + amount) % Lt;
    }

    bool sameIntArray(const multi1d<int>& lhs, const multi1d<int>& rhs)
    {
      if (lhs.size() != rhs.size())
        return false;

      for (int i = 0; i < lhs.size(); ++i)
      {
        if (lhs[i] != rhs[i])
          return false;
      }

      return true;
    }

    bool exactMatch(const ColorMatrix& lhs, const ColorMatrix& rhs)
    {
      return toDouble(norm2(lhs - rhs)) == 0.0;
    }

    multi1d<int> buildBlockSlices(int start, int width, int Lt, const std::string& label)
    {
      multi1d<int> block(width);
      multi1d<bool> seen(Lt);
      seen = false;

      for (int i = 0; i < width; ++i)
      {
        const int t = advanceTime(start, i, Lt);
        if (seen[t])
        {
          std::ostringstream os;
          os << label << " revisits time slice " << t
             << "; frozen_width=" << width << " is invalid for Lt=" << Lt;
          throwGaugeSubdomainSplit(os.str());
        }

        seen[t] = true;
        block[i] = t;
      }

      return block;
    }

    multi1d<int> buildForwardInterval(int start, int length, int Lt)
    {
      multi1d<int> interval(length);
      for (int i = 0; i < length; ++i)
        interval[i] = advanceTime(start, i, Lt);
      return interval;
    }

    multi1d<int> concatTimeOrders(const multi1d<int>& lhs,
                                  const multi1d<int>& middle,
                                  const multi1d<int>& rhs)
    {
      multi1d<int> out(lhs.size() + middle.size() + rhs.size());
      int offset = 0;

      for (int i = 0; i < lhs.size(); ++i)
        out[offset++] = lhs[i];

      for (int i = 0; i < middle.size(); ++i)
        out[offset++] = middle[i];

      for (int i = 0; i < rhs.size(); ++i)
        out[offset++] = rhs[i];

      return out;
    }

    SplitPlan buildSplitPlan(const multi1d<int>& parent_nrow,
                             const GaugeSubdomainSplitParams& input_param)
    {
      SplitPlan plan;
      plan.parent_nrow = parent_nrow;
      plan.param = input_param;

      if (input_param.t_dir < 0 || input_param.t_dir >= Nd)
      {
        std::ostringstream os;
        os << "t_dir=" << input_param.t_dir << " is outside [0," << Nd - 1 << "]";
        throwGaugeSubdomainSplit(os.str());
      }

      plan.Lt = parent_nrow[input_param.t_dir];

      if (input_param.cut0 < 0 || input_param.cut0 >= plan.Lt)
      {
        std::ostringstream os;
        os << "cut0=" << input_param.cut0 << " is outside [0," << plan.Lt - 1 << "]";
        throwGaugeSubdomainSplit(os.str());
      }

      if (input_param.cut1 < 0 || input_param.cut1 >= plan.Lt)
      {
        std::ostringstream os;
        os << "cut1=" << input_param.cut1 << " is outside [0," << plan.Lt - 1 << "]";
        throwGaugeSubdomainSplit(os.str());
      }

      if (input_param.cut0 == input_param.cut1)
        throwGaugeSubdomainSplit("cut0 and cut1 must be different");

      if (input_param.frozen_width <= 0)
      {
        std::ostringstream os;
        os << "frozen_width=" << input_param.frozen_width << " must be positive";
        throwGaugeSubdomainSplit(os.str());
      }

      const int cut_left = std::min(input_param.cut0, input_param.cut1);
      const int cut_right = std::max(input_param.cut0, input_param.cut1);

      plan.param.cut0 = cut_left;
      plan.param.cut1 = cut_right;

      plan.F0 = buildBlockSlices(cut_left, input_param.frozen_width, plan.Lt, "F0");
      plan.F1 = buildBlockSlices(cut_right, input_param.frozen_width, plan.Lt, "F1");

      multi1d<int> owner(plan.Lt);
      owner = 0;

      for (int i = 0; i < plan.F0.size(); ++i)
        owner[plan.F0[i]] = 1;

      for (int i = 0; i < plan.F1.size(); ++i)
      {
        if (owner[plan.F1[i]] != 0)
        {
          std::ostringstream os;
          os << "frozen boundary blocks overlap at parent time slice " << plan.F1[i];
          throwGaugeSubdomainSplit(os.str());
        }

        owner[plan.F1[i]] = 2;
      }

      plan.lenA = (cut_right - (cut_left + input_param.frozen_width) + plan.Lt) % plan.Lt;
      plan.lenB = (cut_left - (cut_right + input_param.frozen_width) + plan.Lt) % plan.Lt;

      if (plan.lenA == 0 || plan.lenB == 0)
      {
        std::ostringstream os;
        os << "cuts " << cut_left << " and " << cut_right
           << " with frozen_width=" << input_param.frozen_width
           << " do not leave two nonempty independent subdomains";
        throwGaugeSubdomainSplit(os.str());
      }

      plan.A = buildForwardInterval(cut_left + input_param.frozen_width, plan.lenA, plan.Lt);
      plan.B = buildForwardInterval(cut_right + input_param.frozen_width, plan.lenB, plan.Lt);
      plan.child0_order = concatTimeOrders(plan.F0, plan.A, plan.F1);
      plan.child1_order = concatTimeOrders(plan.F1, plan.B, plan.F0);

      plan.child0_nrow = parent_nrow;
      plan.child1_nrow = parent_nrow;
      plan.child0_nrow[input_param.t_dir] = plan.child0_order.size();
      plan.child1_nrow[input_param.t_dir] = plan.child1_order.size();

      return plan;
    }

    HostGaugeField extractGaugeField(const multi1d<LatticeColorMatrix>& u,
                                     const multi1d<int>& nrow)
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

    multi1d<LatticeColorMatrix> materializeGaugeField(const HostGaugeField& host)
    {
      multi1d<LatticeColorMatrix> u(Nd);
      for (int mu = 0; mu < Nd; ++mu)
        u[mu] = zero;

      const int vol = latticeVolume(host.nrow);
      for (int site = 0; site < vol; ++site)
      {
        multi1d<int> coord = crtesn(site, host.nrow);
        for (int mu = 0; mu < Nd; ++mu)
          pokeSite(u[mu], host.link[mu][site], coord);
      }

      return u;
    }

    HostGaugeField buildChildHostField(const HostGaugeField& parent,
                                       const multi1d<int>& child_nrow,
                                       const multi1d<int>& child_order,
                                       int t_dir)
    {
      HostGaugeField child;
      child.nrow = child_nrow;
      child.link.resize(Nd);

      const int child_vol = latticeVolume(child_nrow);
      for (int mu = 0; mu < Nd; ++mu)
        child.link[mu].resize(child_vol);

      for (int site = 0; site < child_vol; ++site)
      {
        multi1d<int> child_coord = crtesn(site, child_nrow);
        multi1d<int> parent_coord = child_coord;
        parent_coord[t_dir] = child_order[child_coord[t_dir]];

        const int parent_site = local_site(parent_coord, parent.nrow);
        for (int mu = 0; mu < Nd; ++mu)
          child.link[mu][site] = parent.link[mu][parent_site];
      }

      return child;
    }

    multi1d<int> buildGlobalToLocalMap(const multi1d<int>& local_to_global_t, int Lt)
    {
      multi1d<int> global_to_local(Lt);
      global_to_local = -1;

      for (int local_t = 0; local_t < local_to_global_t.size(); ++local_t)
      {
        const int global_t = local_to_global_t[local_t];
        if (global_to_local[global_t] != -1)
        {
          std::ostringstream os;
          os << "duplicate global time slice " << global_t << " in child local_to_global_t";
          throwGaugeSubdomainSplit(os.str());
        }

        global_to_local[global_t] = local_t;
      }

      return global_to_local;
    }
  }

  GaugeSubdomainSplitChild::GaugeSubdomainSplitChild()
  {
  }

  GaugeSubdomainSplitChild::GaugeSubdomainSplitChild(const GaugeSubdomainSplitChild& rhs) :
    nrow(rhs.nrow),
    local_to_global_t(rhs.local_to_global_t),
    frozen_local_intervals(rhs.frozen_local_intervals)
  {
    if (rhs.u.size() == 0 || rhs.nrow.size() == 0)
      return;

    LayoutGuard guard;
    LayoutGuard::activate(rhs.nrow);

    u.resize(rhs.u.size());
    for (int mu = 0; mu < rhs.u.size(); ++mu)
      u[mu] = rhs.u[mu];
  }

  GaugeSubdomainSplitChild&
  GaugeSubdomainSplitChild::operator=(const GaugeSubdomainSplitChild& rhs)
  {
    if (this == &rhs)
      return *this;

    nrow = rhs.nrow;
    local_to_global_t = rhs.local_to_global_t;
    frozen_local_intervals = rhs.frozen_local_intervals;

    if (rhs.u.size() == 0 || rhs.nrow.size() == 0)
    {
      u.resize(0);
      return *this;
    }

    LayoutGuard guard;
    LayoutGuard::activate(rhs.nrow);

    u.resize(rhs.u.size());
    for (int mu = 0; mu < rhs.u.size(); ++mu)
      u[mu] = rhs.u[mu];

    return *this;
  }

  GaugeSubdomainSplitResult
  splitGaugeSubdomains(const multi1d<LatticeColorMatrix>& parent_u,
                       const GaugeSubdomainSplitParams& input_param)
  {
    START_CODE();

    LayoutGuard guard;

    if (parent_u.size() != Nd)
      throwGaugeSubdomainSplit("parent_u must contain Nd link directions");

    const multi1d<int> parent_nrow = Layout::lattSize();
    const SplitPlan plan = buildSplitPlan(parent_nrow, input_param);
    const HostGaugeField parent_host = extractGaugeField(parent_u, parent_nrow);

    GaugeSubdomainSplitResult result;
    result.param = plan.param;
    result.parent_nrow = parent_nrow;

    result.child0.nrow = plan.child0_nrow;
    result.child0.local_to_global_t = plan.child0_order;
    result.child0.frozen_local_intervals.resize(2);
    result.child0.frozen_local_intervals[0] =
      GaugeSubdomainSplitInterval(0, input_param.frozen_width - 1);
    result.child0.frozen_local_intervals[1] =
      GaugeSubdomainSplitInterval(plan.child0_nrow[input_param.t_dir] - input_param.frozen_width,
                                  plan.child0_nrow[input_param.t_dir] - 1);

    result.child1.nrow = plan.child1_nrow;
    result.child1.local_to_global_t = plan.child1_order;
    result.child1.frozen_local_intervals.resize(2);
    result.child1.frozen_local_intervals[0] =
      GaugeSubdomainSplitInterval(0, input_param.frozen_width - 1);
    result.child1.frozen_local_intervals[1] =
      GaugeSubdomainSplitInterval(plan.child1_nrow[input_param.t_dir] - input_param.frozen_width,
                                  plan.child1_nrow[input_param.t_dir] - 1);

    HostGaugeField child0_host =
      buildChildHostField(parent_host, plan.child0_nrow, plan.child0_order, plan.param.t_dir);
    LayoutGuard::activate(plan.child0_nrow);
    result.child0.u = materializeGaugeField(child0_host);

    HostGaugeField child1_host =
      buildChildHostField(parent_host, plan.child1_nrow, plan.child1_order, plan.param.t_dir);
    LayoutGuard::activate(plan.child1_nrow);
    result.child1.u = materializeGaugeField(child1_host);

    LayoutGuard::activate(parent_nrow);

    END_CODE();
    return result;
  }

  multi1d<LatticeColorMatrix>
  stitchGaugeSubdomains(const GaugeSubdomainSplitResult& split)
  {
    START_CODE();

    LayoutGuard guard;

    const SplitPlan plan = buildSplitPlan(split.parent_nrow, split.param);

    if (!sameIntArray(split.parent_nrow, plan.parent_nrow))
      throwGaugeSubdomainSplit("stored parent_nrow does not match the split plan");

    if (!sameIntArray(split.child0.nrow, plan.child0_nrow))
      throwGaugeSubdomainSplit("child0 nrow metadata does not match the split plan");

    if (!sameIntArray(split.child1.nrow, plan.child1_nrow))
      throwGaugeSubdomainSplit("child1 nrow metadata does not match the split plan");

    if (!sameIntArray(split.child0.local_to_global_t, plan.child0_order))
      throwGaugeSubdomainSplit("child0 local_to_global_t does not match the split plan");

    if (!sameIntArray(split.child1.local_to_global_t, plan.child1_order))
      throwGaugeSubdomainSplit("child1 local_to_global_t does not match the split plan");

    if (split.child0.u.size() != Nd)
      throwGaugeSubdomainSplit("child0 gauge field must contain Nd link directions");

    if (split.child1.u.size() != Nd)
      throwGaugeSubdomainSplit("child1 gauge field must contain Nd link directions");

    if (split.child0.frozen_local_intervals.size() != 2 ||
        split.child0.frozen_local_intervals[0].t_start != 0 ||
        split.child0.frozen_local_intervals[0].t_end != split.param.frozen_width - 1 ||
        split.child0.frozen_local_intervals[1].t_start != plan.child0_nrow[split.param.t_dir] - split.param.frozen_width ||
        split.child0.frozen_local_intervals[1].t_end != plan.child0_nrow[split.param.t_dir] - 1)
      throwGaugeSubdomainSplit("child0 frozen_local_intervals do not match the split plan");

    if (split.child1.frozen_local_intervals.size() != 2 ||
        split.child1.frozen_local_intervals[0].t_start != 0 ||
        split.child1.frozen_local_intervals[0].t_end != split.param.frozen_width - 1 ||
        split.child1.frozen_local_intervals[1].t_start != plan.child1_nrow[split.param.t_dir] - split.param.frozen_width ||
        split.child1.frozen_local_intervals[1].t_end != plan.child1_nrow[split.param.t_dir] - 1)
      throwGaugeSubdomainSplit("child1 frozen_local_intervals do not match the split plan");

    LayoutGuard::activate(plan.child0_nrow);
    HostGaugeField child0_host = extractGaugeField(split.child0.u, plan.child0_nrow);

    LayoutGuard::activate(plan.child1_nrow);
    HostGaugeField child1_host = extractGaugeField(split.child1.u, plan.child1_nrow);

    const multi1d<int> child0_global_to_local = buildGlobalToLocalMap(plan.child0_order, plan.Lt);
    const multi1d<int> child1_global_to_local = buildGlobalToLocalMap(plan.child1_order, plan.Lt);

    HostGaugeField parent_host;
    parent_host.nrow = split.parent_nrow;
    parent_host.link.resize(Nd);
    const int parent_vol = latticeVolume(split.parent_nrow);
    for (int mu = 0; mu < Nd; ++mu)
      parent_host.link[mu].resize(parent_vol);

    for (int site = 0; site < parent_vol; ++site)
    {
      multi1d<int> parent_coord = crtesn(site, split.parent_nrow);
      const int global_t = parent_coord[split.param.t_dir];
      const int child0_local_t = child0_global_to_local[global_t];
      const int child1_local_t = child1_global_to_local[global_t];

      if (child0_local_t < 0 && child1_local_t < 0)
      {
        std::ostringstream os;
        os << "parent time slice " << global_t
           << " is missing from both children during stitch";
        throwGaugeSubdomainSplit(os.str());
      }

      multi1d<int> child0_coord = parent_coord;
      multi1d<int> child1_coord = parent_coord;
      if (child0_local_t >= 0)
        child0_coord[split.param.t_dir] = child0_local_t;
      if (child1_local_t >= 0)
        child1_coord[split.param.t_dir] = child1_local_t;

      const int child0_site =
        (child0_local_t >= 0) ? local_site(child0_coord, plan.child0_nrow) : -1;
      const int child1_site =
        (child1_local_t >= 0) ? local_site(child1_coord, plan.child1_nrow) : -1;

      for (int mu = 0; mu < Nd; ++mu)
      {
        if (child0_site >= 0 && child1_site >= 0)
        {
          const ColorMatrix& lhs = child0_host.link[mu][child0_site];
          const ColorMatrix& rhs = child1_host.link[mu][child1_site];

          if (!exactMatch(lhs, rhs))
          {
            std::ostringstream os;
            os << "duplicated frozen-boundary link mismatch at global t=" << global_t
               << ", mu=" << mu;
            throwGaugeSubdomainSplit(os.str());
          }

          parent_host.link[mu][site] = lhs;
        }
        else if (child0_site >= 0)
        {
          parent_host.link[mu][site] = child0_host.link[mu][child0_site];
        }
        else
        {
          parent_host.link[mu][site] = child1_host.link[mu][child1_site];
        }
      }
    }

    LayoutGuard::activate(split.parent_nrow);
    multi1d<LatticeColorMatrix> parent_u = materializeGaugeField(parent_host);
    guard.setRestore(split.parent_nrow);

    END_CODE();
    return parent_u;
  }

}
