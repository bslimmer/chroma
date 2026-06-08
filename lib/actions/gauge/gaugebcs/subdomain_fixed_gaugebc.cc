/*! \file
 *  \brief Runtime child-boundary gauge BC for subdomain HMC
 */

#include "actions/gauge/gaugebcs/subdomain_fixed_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

#include <sstream>

namespace Chroma
{

  namespace SubdomainFixedGaugeBCEnv
  {
    const std::string name = "SUBDOMAIN_FIXED_GAUGEBC";
  }

  namespace
  {
    Handle<SubdomainFixedGaugeBC> runtime_prototype;
    bool runtime_prototype_valid = false;

    void abortSubdomainFixedGaugeBC(const std::string& message)
    {
      QDPIO::cerr << SubdomainFixedGaugeBCEnv::name << ": " << message << std::endl;
      QDP_abort(1);
    }

    int readLoopExtent(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);
      int loop_extent = 1;
      if (paramtop.count("loop_extent") != 0)
        read(paramtop, "loop_extent", loop_extent);
      return loop_extent;
    }
  }

  namespace SubdomainFixedGaugeBCEnv
  {
    GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >*
    createGaugeBC(XMLReader& xml, const std::string& path)
    {
      if (!runtime_prototype_valid)
        abortSubdomainFixedGaugeBC("runtime prototype is not set");

      const int requested_loop_extent = readLoopExtent(xml, path);
      if (requested_loop_extent != runtime_prototype->getLoopExtent())
      {
        std::ostringstream os;
        os << "loop_extent mismatch between XML (" << requested_loop_extent
           << ") and runtime prototype (" << runtime_prototype->getLoopExtent() << ")";
        abortSubdomainFixedGaugeBC(os.str());
      }

      return new SubdomainFixedGaugeBC(*runtime_prototype);
    }

    static bool registered = false;

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

    void setRuntimePrototype(const Handle<SubdomainFixedGaugeBC>& bc)
    {
      runtime_prototype = bc;
      runtime_prototype_valid = true;
    }

    void clearRuntimePrototype()
    {
      runtime_prototype = Handle<SubdomainFixedGaugeBC>();
      runtime_prototype_valid = false;
    }
  }

  SubdomainFixedGaugeBC::SubdomainFixedGaugeBC(const multi1d<LatticeColorMatrix>& child_u,
                                               const multi1d<GaugeSubdomainSplitInterval>& frozen_intervals_,
                                               int t_dir_,
                                               int loop_extent_) :
    t_dir(t_dir_),
    loop_extent(loop_extent_),
    frozen_intervals(frozen_intervals_),
    fixed_links(child_u)
  {
    buildMasks();
  }

  SubdomainFixedGaugeBC::SubdomainFixedGaugeBC(const SubdomainFixedGaugeBC& rhs) :
    t_dir(rhs.t_dir),
    loop_extent(rhs.loop_extent),
    frozen_intervals(rhs.frozen_intervals),
    fixed_links(rhs.fixed_links),
    frozen_mask(rhs.frozen_mask),
    edge_temporal_mask(rhs.edge_temporal_mask)
  {
  }

  void SubdomainFixedGaugeBC::buildMasks()
  {
    if (t_dir < 0 || t_dir >= Nd)
    {
      std::ostringstream os;
      os << "t_dir=" << t_dir << " is outside [0," << Nd - 1 << "]";
      abortSubdomainFixedGaugeBC(os.str());
    }

    if (loop_extent < 1)
    {
      std::ostringstream os;
      os << "loop_extent must be positive, got " << loop_extent;
      abortSubdomainFixedGaugeBC(os.str());
    }

    if (fixed_links.size() != Nd)
      abortSubdomainFixedGaugeBC("child gauge field must contain Nd link directions");

    const int t_extent = Layout::lattSize()[t_dir];
    if (t_extent < 2)
      abortSubdomainFixedGaugeBC("child temporal extent must be at least 2");

    if (frozen_intervals.size() == 0)
      abortSubdomainFixedGaugeBC("frozen_intervals must not be empty");

    LatticeInteger t = Layout::latticeCoordinate(t_dir);
    LatticeBoolean site_mask = false;

    for (int i = 0; i < frozen_intervals.size(); ++i)
    {
      const GaugeSubdomainSplitInterval& interval = frozen_intervals[i];
      if (interval.t_start < 0 || interval.t_start > interval.t_end || interval.t_end >= t_extent)
      {
        std::ostringstream os;
        os << "invalid frozen interval [" << interval.t_start << "," << interval.t_end
           << "] for child temporal extent " << t_extent;
        abortSubdomainFixedGaugeBC(os.str());
      }

      site_mask |= (t >= interval.t_start) && (t <= interval.t_end);
    }

    frozen_mask.resize(Nd);
    for (int mu = 0; mu < Nd; ++mu)
      frozen_mask[mu] = site_mask;

    edge_temporal_mask = (t == (t_extent - 1));
  }

  void SubdomainFixedGaugeBC::modify(multi1d<LatticeColorMatrix>& u) const
  {
    START_CODE();

    if (u.size() != Nd)
      abortSubdomainFixedGaugeBC("modify() received a gauge field with the wrong number of directions");

    clampToFrozen(u);

    LatticeColorMatrix z = QDP::zero;
    copymask(u[t_dir], edge_temporal_mask, z);

    END_CODE();
  }

  void SubdomainFixedGaugeBC::clampToFrozen(multi1d<LatticeColorMatrix>& u) const
  {
    START_CODE();

    if (u.size() != Nd)
      abortSubdomainFixedGaugeBC("clampToFrozen() received a gauge field with the wrong number of directions");

    for (int mu = 0; mu < Nd; ++mu)
      copymask(u[mu], frozen_mask[mu], fixed_links[mu]);

    END_CODE();
  }

  void SubdomainFixedGaugeBC::zero(multi1d<LatticeColorMatrix>& ds_u) const
  {
    START_CODE();

    if (ds_u.size() != Nd)
      abortSubdomainFixedGaugeBC("zero() received a gauge-like field with the wrong number of directions");

    LatticeColorMatrix z = QDP::zero;
    for (int mu = 0; mu < Nd; ++mu)
      copymask(ds_u[mu], frozen_mask[mu], z);

    END_CODE();
  }

  bool SubdomainFixedGaugeBC::nontrivialP() const
  {
    return true;
  }

}
