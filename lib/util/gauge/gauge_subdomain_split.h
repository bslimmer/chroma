// -*- C++ -*-
/*! \file
 *  \brief Gauge subdomain split/stitch utilities
 */

#ifndef __gauge_subdomain_split_h__
#define __gauge_subdomain_split_h__

#include "chromabase.h"

namespace Chroma
{

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitParams
  {
    GaugeSubdomainSplitParams();
    GaugeSubdomainSplitParams(XMLReader& xml, const std::string& path);

    int t_dir;
    int cut0;
    int cut1;
    int frozen_width;
  };

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitInterval
  {
    GaugeSubdomainSplitInterval();

    int t_start;
    int t_end;
  };

  /*! @ingroup gauge */
  struct GaugeSubdomainSplitPlan
  {
    GaugeSubdomainSplitPlan();
    GaugeSubdomainSplitPlan(XMLReader& xml, const std::string& path);

    GaugeSubdomainSplitParams param;
    multi1d<int> parent_nrow;
    multi1d<int> child0_nrow;
    multi1d<int> child1_nrow;
    multi1d<int> child0_local_to_global_t;
    multi1d<int> child1_local_to_global_t;
    multi1d<GaugeSubdomainSplitInterval> child0_frozen_local_intervals;
    multi1d<GaugeSubdomainSplitInterval> child1_frozen_local_intervals;
  };

  /*! @ingroup gauge */
  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitParams& p);

  /*! @ingroup gauge */
  void write(XMLWriter& xml, const std::string& path, const GaugeSubdomainSplitParams& p);

  /*! @ingroup gauge */
  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitInterval& p);

  /*! @ingroup gauge */
  void write(XMLWriter& xml, const std::string& path, const GaugeSubdomainSplitInterval& p);

  /*! @ingroup gauge */
  void read(XMLReader& xml, const std::string& path, GaugeSubdomainSplitPlan& p);

  /*! @ingroup gauge */
  void write(XMLWriter& xml, const std::string& path, const GaugeSubdomainSplitPlan& p);

  /*! @ingroup gauge */
  GaugeSubdomainSplitPlan
  makeGaugeSubdomainSplitPlan(const multi1d<int>& parent_nrow,
                              const GaugeSubdomainSplitParams& param);

  /*! @ingroup gauge */
  void
  snapshotGaugeField(const multi1d<LatticeColorMatrix>& u,
                     multi1d<ColorMatrix>& site_links);

  /*! @ingroup gauge */
  void
  materializeGaugeField(const multi1d<int>& nrow,
                        const multi1d<ColorMatrix>& site_links,
                        multi1d<LatticeColorMatrix>& u);

  /*! @ingroup gauge */
  void
  extractGaugeSubdomain(const multi1d<ColorMatrix>& parent_links,
                        const GaugeSubdomainSplitPlan& plan,
                        int child_id,
                        multi1d<ColorMatrix>& child_links);

  /*! @ingroup gauge */
  void
  stitchGaugeSubdomains(const GaugeSubdomainSplitPlan& plan,
                        const multi1d<ColorMatrix>& child0_links,
                        const multi1d<ColorMatrix>& child1_links,
                        multi1d<ColorMatrix>& parent_links);

}

#endif
