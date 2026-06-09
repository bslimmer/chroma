// -*- C++ -*-
/*! \file
 * \brief Autodiscover a gauge-monomial BC source for refreshed-momentum masking
 */

#ifndef HMC_GAUGE_MONOMIAL_MOMENTUM_BC_H
#define HMC_GAUGE_MONOMIAL_MOMENTUM_BC_H

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/monomial/gauge_monomial.h"

namespace Chroma
{
  //! Discover a compatible nontrivial gauge monomial BC source from Hamiltonian monomial ids
  Handle<GaugeMonomial>
  discoverMomentumMaskingGaugeMonomial(const multi1d<std::string>& monomial_ids);
}

#endif
