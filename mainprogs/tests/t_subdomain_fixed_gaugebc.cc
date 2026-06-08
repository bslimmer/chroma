/*! \file
 *  \brief Focused checks for the subdomain fixed gauge BC
 */

#include "chroma.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugebcs/subdomain_fixed_gaugebc.h"
#include "actions/gauge/gaugestates/simple_gaugestate.h"
#include "actions/gauge/gaugestates/periodic_gaugestate.h"
#include "util/gauge/gauge_subdomain_split.h"

#include <cmath>
#include <sstream>

using namespace Chroma;

namespace
{
  class LayoutGuard
  {
  public:
    LayoutGuard() : saved_nrow(Layout::lattSize()) {}

    ~LayoutGuard()
    {
      if (saved_nrow.size() > 0)
      {
        Layout::setLattSize(saved_nrow);
        Layout::create();
      }
    }

    static void activate(const multi1d<int>& nrow)
    {
      Layout::setLattSize(nrow);
      Layout::create();
    }

  private:
    multi1d<int> saved_nrow;
  };

  void fail(const std::string& message)
  {
    QDPIO::cerr << "FAIL: " << message << std::endl;
    QDP_abort(1);
  }

  LatticeBoolean frozenMask(const int t_dir,
                            const multi1d<GaugeSubdomainSplitInterval>& intervals)
  {
    LatticeInteger t = Layout::latticeCoordinate(t_dir);
    LatticeBoolean mask = false;
    for (int i = 0; i < intervals.size(); ++i)
      mask |= (t >= intervals[i].t_start) && (t <= intervals[i].t_end);
    return mask;
  }

  LatticeBoolean edgeMask(const int t_dir)
  {
    LatticeInteger t = Layout::latticeCoordinate(t_dir);
    return (t == (Layout::lattSize()[t_dir] - 1));
  }

  Double maskedNorm2(const LatticeColorMatrix& field, const LatticeBoolean& mask)
  {
    LatticeColorMatrix tmp = zero;
    copymask(tmp, mask, field);
    return sum(localNorm2(tmp));
  }

  Double plaquetteAction(const Handle<CreateGaugeState<multi1d<LatticeColorMatrix>,
                                                       multi1d<LatticeColorMatrix> > >& cgs,
                         const multi1d<LatticeColorMatrix>& u)
  {
    AnisoParam_t aniso;
    PlaqGaugeAct act(cgs, Real(1), aniso);
    Handle<GaugeState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state(
      act.createState(u));
    return act.S(state);
  }

  Double rectAction(const Handle<CreateGaugeState<multi1d<LatticeColorMatrix>,
                                                  multi1d<LatticeColorMatrix> > >& cgs,
                    const multi1d<LatticeColorMatrix>& u)
  {
    RectGaugeAct act(cgs, Real(1));
    Handle<GaugeState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state(
      act.createState(u));
    return act.S(state);
  }
}

int main(int argc, char* argv[])
{
  Chroma::initialize(&argc, &argv);
  START_CODE();

  {
    multi1d<int> parent_nrow(Nd);
    parent_nrow[0] = 4;
    parent_nrow[1] = 4;
    parent_nrow[2] = 4;
    parent_nrow[3] = 8;
    Layout::setLattSize(parent_nrow);
    Layout::create();

    multi1d<LatticeColorMatrix> parent_u(Nd);
    for (int mu = 0; mu < Nd; ++mu)
      parent_u[mu] = 1;

    GaugeSubdomainSplitParams split_param;
    split_param.t_dir = 3;
    split_param.cut0 = 1;
    split_param.cut1 = 5;
    split_param.frozen_width = 1;

    GaugeSubdomainSplitResult split = splitGaugeSubdomains(parent_u, split_param);

    {
      LayoutGuard guard;
      LayoutGuard::activate(split.child0.nrow);

      Handle<GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > bc(
        new SubdomainFixedGaugeBC(split.child0.u,
                                  split.child0.frozen_local_intervals,
                                  split_param.t_dir,
                                  2));
      Handle<CreateGaugeState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >
        subdomain_cgs(new CreateSimpleGaugeState<multi1d<LatticeColorMatrix>,
                                                 multi1d<LatticeColorMatrix> >(bc));
      Handle<CreateGaugeState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >
        periodic_cgs(new CreatePeriodicGaugeState<multi1d<LatticeColorMatrix>,
                                                  multi1d<LatticeColorMatrix> >());

      multi1d<LatticeColorMatrix> restored = split.child0.u;
      gaussian(restored[0]);
      bc->modify(restored);

      multi1d<LatticeColorMatrix> clamped = split.child0.u;
      gaussian(clamped[0]);
      gaussian(clamped[split_param.t_dir]);
      dynamic_cast<SubdomainFixedGaugeBC&>(*bc).clampToFrozen(clamped);

      const LatticeBoolean frozen =
        frozenMask(split_param.t_dir, split.child0.frozen_local_intervals);
      const LatticeBoolean edge = edgeMask(split_param.t_dir);

      if (toDouble(maskedNorm2(restored[0] - split.child0.u[0], frozen)) != 0.0)
        fail("modify() did not restore frozen links exactly");

      if (toDouble(maskedNorm2(clamped[0] - split.child0.u[0], frozen)) != 0.0)
        fail("clampToFrozen() did not restore frozen links exactly");

      if (toDouble(maskedNorm2(restored[split_param.t_dir], edge)) != 0.0)
        fail("modify() did not zero copied edge temporal links in the gauge state");

      if (toDouble(maskedNorm2(clamped[split_param.t_dir] - split.child0.u[split_param.t_dir],
                               edge)) != 0.0)
        fail("clampToFrozen() changed stored edge temporal links");

      multi1d<LatticeColorMatrix> changed_edge = split.child0.u;
      LatticeColorMatrix edge_tweak;
      gaussian(edge_tweak);
      copymask(changed_edge[split_param.t_dir], edge, edge_tweak);

      const Double periodic_plaq_before = plaquetteAction(periodic_cgs, split.child0.u);
      const Double periodic_plaq_after = plaquetteAction(periodic_cgs, changed_edge);
      if (std::fabs(toDouble(periodic_plaq_after - periodic_plaq_before)) < 1.0e-12)
        fail("periodic plaquette action did not see changed edge temporal links");

      const Double subdomain_plaq_before = plaquetteAction(subdomain_cgs, split.child0.u);
      const Double subdomain_plaq_after = plaquetteAction(subdomain_cgs, changed_edge);
      if (std::fabs(toDouble(subdomain_plaq_after - subdomain_plaq_before)) > 1.0e-12)
        fail("subdomain plaquette action still depends on copied edge temporal links");

      const Double periodic_rect_before = rectAction(periodic_cgs, split.child0.u);
      const Double periodic_rect_after = rectAction(periodic_cgs, changed_edge);
      if (std::fabs(toDouble(periodic_rect_after - periodic_rect_before)) < 1.0e-12)
        fail("periodic rectangle action did not see changed edge temporal links");

      const Double subdomain_rect_before = rectAction(subdomain_cgs, split.child0.u);
      const Double subdomain_rect_after = rectAction(subdomain_cgs, changed_edge);
      if (std::fabs(toDouble(subdomain_rect_after - subdomain_rect_before)) > 1.0e-12)
        fail("subdomain rectangle action still depends on copied edge temporal links");

      multi1d<LatticeColorMatrix> changed_bulk = split.child0.u;
      LatticeBoolean bulk_mask = (Layout::latticeCoordinate(split_param.t_dir) ==
                                  (Layout::lattSize()[split_param.t_dir] - 2));
      LatticeColorMatrix bulk_tweak;
      gaussian(bulk_tweak);
      copymask(changed_bulk[split_param.t_dir], bulk_mask, bulk_tweak);

      const Double subdomain_bulk_before = plaquetteAction(subdomain_cgs, split.child0.u);
      const Double subdomain_bulk_after = plaquetteAction(subdomain_cgs, changed_bulk);
      if (std::fabs(toDouble(subdomain_bulk_after - subdomain_bulk_before)) < 1.0e-12)
        fail("subdomain plaquette action ignored a non-edge temporal link change");
    }
  }

  QDPIO::cout << "PASS: t_subdomain_fixed_gaugebc" << std::endl;

  END_CODE();
  Chroma::finalize();
  return 0;
}
