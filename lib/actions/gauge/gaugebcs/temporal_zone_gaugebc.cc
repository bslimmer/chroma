/*! \file
 *  \brief Temporal zone gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/temporal_zone_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

#include <sstream>

namespace Chroma
{

  namespace TemporalZoneGaugeBCEnv
  {
    const std::string name = "TEMPORAL_ZONE_GAUGEBC";
  }

  namespace
  {
    void abortTemporalZoneGaugeBC(const std::string& message)
    {
      QDPIO::cerr << TemporalZoneGaugeBCEnv::name << ": " << message << std::endl;
      QDP_abort(1);
    }

    std::string intervalAsString(const TemporalZoneInterval& interval)
    {
      std::ostringstream os;
      os << "[" << interval.t_start << "," << interval.t_end << "]";
      return os.str();
    }
  }

  namespace TemporalZoneGaugeBCEnv
  {
    //! Callback function to register with the factory
    GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >*
    createGaugeBC(XMLReader& xml, const std::string& path)
    {
      return new TemporalZoneGaugeBC(TemporalZoneGaugeBCParams(xml, path));
    }

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
        success &= TheGaugeBCFactory::Instance().registerObject(name, createGaugeBC);
        registered = true;
      }
      return success;
    }
  }

  TemporalZoneGaugeBCParams::TemporalZoneGaugeBCParams()
  {
    t_dir = Nd - 1;
  }

  TemporalZoneGaugeBCParams::TemporalZoneGaugeBCParams(XMLReader& xml,
                                                       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    try
    {
      t_dir = Nd - 1;
      if (paramtop.count("t_dir") != 0)
        read(paramtop, "t_dir", t_dir);

      if (paramtop.count("zero_intervals") == 0)
        abortTemporalZoneGaugeBC("missing required zero_intervals block");

      read(paramtop, "zero_intervals", zero_intervals);

      if (zero_intervals.size() == 0)
        abortTemporalZoneGaugeBC("zero_intervals must contain at least one interval");
    }
    catch (const std::string& e)
    {
      abortTemporalZoneGaugeBC("error reading XML: " + e);
    }
  }

  void read(XMLReader& xml, const std::string& path, TemporalZoneInterval& p)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "t_start", p.t_start);
    read(paramtop, "t_end", p.t_end);
  }

  void write(XMLWriter& xml, const std::string& path, const TemporalZoneInterval& p)
  {
    push(xml, path);
    write(xml, "t_start", p.t_start);
    write(xml, "t_end", p.t_end);
    pop(xml);
  }

  void read(XMLReader& xml, const std::string& path, TemporalZoneGaugeBCParams& p)
  {
    TemporalZoneGaugeBCParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const TemporalZoneGaugeBCParams& p)
  {
    push(xml, path);
    write(xml, "t_dir", p.t_dir);
    write(xml, "zero_intervals", p.zero_intervals);
    pop(xml);
  }

  TemporalZoneGaugeBC::TemporalZoneGaugeBC(const TemporalZoneGaugeBCParams& p) :
    param(p)
  {
    if (param.t_dir < 0 || param.t_dir >= Nd)
    {
      std::ostringstream os;
      os << "t_dir=" << param.t_dir << " is outside [0," << Nd - 1 << "]";
      abortTemporalZoneGaugeBC(os.str());
    }

    if (param.zero_intervals.size() == 0)
      abortTemporalZoneGaugeBC("zero_intervals must contain at least one interval");

    const int t_extent = Layout::lattSize()[param.t_dir];

    for (int i = 0; i < param.zero_intervals.size(); ++i)
    {
      const TemporalZoneInterval& interval = param.zero_intervals[i];

      if (interval.t_start < 0 || interval.t_start > interval.t_end || interval.t_end >= t_extent)
      {
        std::ostringstream os;
        os << "interval " << i << " " << intervalAsString(interval)
           << " is invalid for t_dir=" << param.t_dir
           << " with extent " << t_extent;
        abortTemporalZoneGaugeBC(os.str());
      }
    }

    mask.resize(Nd);

    LatticeInteger t = Layout::latticeCoordinate(param.t_dir);
    LatticeBoolean site_mask = false;

    for (int i = 0; i < param.zero_intervals.size(); ++i)
    {
      const TemporalZoneInterval& interval = param.zero_intervals[i];
      site_mask |= (t >= interval.t_start) && (t <= interval.t_end);
    }

    for (int mu = 0; mu < Nd; ++mu)
      mask[mu] = site_mask;
  }

  void TemporalZoneGaugeBC::modify(multi1d<LatticeColorMatrix>& u) const
  {
    (void)u;
  }

  void TemporalZoneGaugeBC::zero(multi1d<LatticeColorMatrix>& ds_u) const
  {
    START_CODE();

    LatticeColorMatrix z = QDP::zero;

    for (int mu = 0; mu < ds_u.size(); ++mu)
      copymask(ds_u[mu], mask[mu], z);

    END_CODE();
  }

  bool TemporalZoneGaugeBC::nontrivialP() const
  {
    return param.zero_intervals.size() > 0;
  }

}
