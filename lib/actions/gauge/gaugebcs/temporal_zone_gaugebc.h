// -*- C++ -*-
/*! \file
 *  \brief Temporal zone gauge boundary conditions
 */

#ifndef __temporal_zone_gaugebc_h__
#define __temporal_zone_gaugebc_h__

#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugebcs */
  namespace TemporalZoneGaugeBCEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  /*! @ingroup gaugebcs */
  struct TemporalZoneInterval
  {
    TemporalZoneInterval() : t_start(0), t_end(-1) {}

    int t_start;
    int t_end;
  };

  /*! @ingroup gaugebcs */
  struct TemporalZoneGaugeBCParams
  {
    TemporalZoneGaugeBCParams();
    TemporalZoneGaugeBCParams(XMLReader& xml, const std::string& path);

    int t_dir;
    multi1d<TemporalZoneInterval> zero_intervals;
  };

  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, TemporalZoneInterval& p);

  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const TemporalZoneInterval& p);

  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, TemporalZoneGaugeBCParams& p);

  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const TemporalZoneGaugeBCParams& p);

  //! Concrete class for suppressing gauge-like fields on temporal zones
  /*! @ingroup gaugebcs */
  class TemporalZoneGaugeBC :
    public GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Only full constructor
    TemporalZoneGaugeBC(const TemporalZoneGaugeBCParams& p);

    //! Destructor is automatic
    ~TemporalZoneGaugeBC() {}

    //! This BC only suppresses force/update fields
    void modify(multi1d<LatticeColorMatrix>& u) const;

    //! Zero the gauge-like field on selected time zones
    void zero(multi1d<LatticeColorMatrix>& ds_u) const;

    //! Says if there are zeroed links within the lattice
    bool nontrivialP() const;

  private:
    TemporalZoneGaugeBC();
    void operator=(const TemporalZoneGaugeBC&);

  private:
    TemporalZoneGaugeBCParams param;
    multi1d<LatticeBoolean> mask;
  };

}

#endif
