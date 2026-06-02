// -*- C++ -*-
/*! \file
 * \brief Split and stitch utilities for temporal gauge subdomains
 */

#ifndef __gauge_subdomain_split_h__
#define __gauge_subdomain_split_h__

#include "chromabase.h"

namespace Chroma
{

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitParams
  {
    GaugeSubdomainSplitParams() : t_dir(Nd - 1), cut0(0), cut1(0), frozen_width(1) {}

    int t_dir;
    int cut0;
    int cut1;
    int frozen_width;
  };

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitInterval
  {
    GaugeSubdomainSplitInterval() : t_start(0), t_end(-1) {}
    GaugeSubdomainSplitInterval(int start_, int end_) : t_start(start_), t_end(end_) {}

    int t_start;
    int t_end;
  };

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitChild
  {
    GaugeSubdomainSplitChild();
    GaugeSubdomainSplitChild(const GaugeSubdomainSplitChild& rhs);
    GaugeSubdomainSplitChild& operator=(const GaugeSubdomainSplitChild& rhs);

    // The gauge field is allocated under the lattice extent stored in nrow.
    // Callers must activate that matching Layout before operating on u.
    multi1d<LatticeColorMatrix> u;
    multi1d<int> nrow;
    multi1d<int> local_to_global_t;
    multi1d<GaugeSubdomainSplitInterval> frozen_local_intervals;
  };

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitResult
  {
    GaugeSubdomainSplitParams param;
    multi1d<int> parent_nrow;
    GaugeSubdomainSplitChild child0;
    GaugeSubdomainSplitChild child1;
  };

  /*! @ingroup gauge */
  GaugeSubdomainSplitResult
  splitGaugeSubdomains(const multi1d<LatticeColorMatrix>& parent_u,
                       const GaugeSubdomainSplitParams& param);

  /*! @ingroup gauge */
  multi1d<LatticeColorMatrix>
  stitchGaugeSubdomains(const GaugeSubdomainSplitResult& split);

}

#endif
