// -*- C++ -*-
/*! \file
 * \brief Autodiscover a gauge-monomial BC source for refreshed-momentum masking
 */

#include "chromabase.h"

#include "update/molecdyn/hmc/gauge_monomial_momentum_bc.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "meas/inline/io/named_objmap.h"

#include <sstream>

namespace Chroma
{
  namespace
  {
    typedef multi1d<LatticeColorMatrix> LCM;
    typedef Monomial<LCM, LCM> BaseMonomial;

    bool compatibleMomentumMasks(const GaugeMonomial& lhs,
                                 const GaugeMonomial& rhs)
    {
      multi1d<LatticeColorMatrix> probe(Nd);
      multi1d<LatticeColorMatrix> lhs_probe(Nd);
      multi1d<LatticeColorMatrix> rhs_probe(Nd);

      for (int mu = 0; mu < Nd; ++mu)
      {
        probe[mu] = Real(mu + 1);
        lhs_probe[mu] = probe[mu];
        rhs_probe[mu] = probe[mu];
      }

      lhs.zeroGaugeLikeField(lhs_probe);
      rhs.zeroGaugeLikeField(rhs_probe);

      for (int mu = 0; mu < Nd; ++mu)
      {
        if (toBool(norm2(lhs_probe[mu] - rhs_probe[mu]) > Double(0)))
          return false;
      }

      return true;
    }
  }

  Handle<GaugeMonomial>
  discoverMomentumMaskingGaugeMonomial(const multi1d<std::string>& monomial_ids)
  {
    Handle<GaugeMonomial> candidate;
    std::string candidate_id;

    for (int i = 0; i < monomial_ids.size(); ++i)
    {
      const std::string& monomial_id = monomial_ids[i];
      Handle<BaseMonomial> monomial_handle =
        TheNamedObjMap::Instance().getData< Handle<BaseMonomial> >(monomial_id);

      GaugeMonomial* gauge_monomial_ptr =
        dynamic_cast<GaugeMonomial*>(monomial_handle.operator->());

      if (gauge_monomial_ptr == 0)
        continue;

      Handle<GaugeMonomial> gauge_monomial = monomial_handle.cast<GaugeMonomial>();
      if (!gauge_monomial->hasNontrivialGaugeBC())
        continue;

      if (candidate.operator->() == 0)
      {
        candidate = gauge_monomial;
        candidate_id = monomial_id;
        continue;
      }

      if (!compatibleMomentumMasks(*candidate, *gauge_monomial))
      {
        std::ostringstream error;
        error << "discoverMomentumMaskingGaugeMonomial: incompatible nontrivial "
              << "gauge BC masks for monomial ids '" << candidate_id
              << "' and '" << monomial_id << "'";
        throw error.str();
      }
    }

    return candidate;
  }
}
