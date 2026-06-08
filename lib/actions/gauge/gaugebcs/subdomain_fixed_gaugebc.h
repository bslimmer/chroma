// -*- C++ -*-
/*! \file
 *  \brief Runtime child-boundary gauge BC for subdomain HMC
 */

#ifndef __subdomain_fixed_gaugebc_h__
#define __subdomain_fixed_gaugebc_h__

#include "gaugebc.h"
#include "handle.h"
#include "util/gauge/gauge_subdomain_split.h"

namespace Chroma
{

  class SubdomainFixedGaugeBC;

  /*! @ingroup gaugebcs */
  namespace SubdomainFixedGaugeBCEnv
  {
    extern const std::string name;
    bool registerAll();

    void setRuntimePrototype(const Handle<SubdomainFixedGaugeBC>& bc);
    void clearRuntimePrototype();
  }

  //! Child fixed-boundary BC for subdomain HMC
  /*! @ingroup gaugebcs */
  class SubdomainFixedGaugeBC :
    public GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    SubdomainFixedGaugeBC(const multi1d<LatticeColorMatrix>& child_u,
                          const multi1d<GaugeSubdomainSplitInterval>& frozen_intervals,
                          int t_dir,
                          int loop_extent);

    SubdomainFixedGaugeBC(const SubdomainFixedGaugeBC& rhs);

    ~SubdomainFixedGaugeBC() {}

    void modify(multi1d<LatticeColorMatrix>& u) const;
    void clampToFrozen(multi1d<LatticeColorMatrix>& u) const;
    void zero(multi1d<LatticeColorMatrix>& ds_u) const;
    bool nontrivialP() const;

    int getDir() const { return t_dir; }
    int getLoopExtent() const { return loop_extent; }
    const multi1d<GaugeSubdomainSplitInterval>& getFrozenIntervals() const { return frozen_intervals; }

  private:
    SubdomainFixedGaugeBC();
    void operator=(const SubdomainFixedGaugeBC&);

    void buildMasks();

  private:
    int t_dir;
    int loop_extent;
    multi1d<GaugeSubdomainSplitInterval> frozen_intervals;
    multi1d<LatticeColorMatrix> fixed_links;
    multi1d<LatticeBoolean> frozen_mask;
    LatticeBoolean edge_temporal_mask;
  };

}

#endif
