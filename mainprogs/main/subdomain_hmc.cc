/*! \file
 *  \brief Main code for quenched subdomain HMC
 */

#include "chroma.h"
#include "actions/gauge/gaugebcs/subdomain_fixed_gaugebc.h"
#include "io/monomial_io.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/io/default_gauge_field.h"
#include "meas/inline/io/named_objmap.h"
#include "update/molecdyn/hmc/lcm_hmc.h"
#include "util/gauge/gauge_subdomain_split.h"
#include "util/gauge/szinqio_gauge_init.h"

#include <fstream>
#include <sstream>

using namespace Chroma;

namespace Chroma
{

  class LayoutGuard
  {
  public:
    LayoutGuard() : saved_nrow(Layout::lattSize()) {}

    ~LayoutGuard()
    {
      if (saved_nrow.size() > 0)
        activate(saved_nrow);
    }

    static void activate(const multi1d<int>& nrow)
    {
      Layout::setLattSize(nrow);
      Layout::create();
    }

  private:
    multi1d<int> saved_nrow;
  };

  class SubdomainLatColMatHMCTrj : public LatColMatHMCTrj
  {
  public:
    SubdomainLatColMatHMCTrj(
      Handle<AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_MC,
      Handle<AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& MD_int,
      const GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& bc_) :
      LatColMatHMCTrj(H_MC, MD_int),
      bc(bc_)
    {
    }

  protected:
    void refreshP(AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
    {
      LatColMatHMCTrj::refreshP(s);
      bc.zero(s.getP());
    }

  private:
    const GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& bc;
  };

  std::string renameOuterTag(const std::string& xml_string,
                             const std::string& old_tag,
                             const std::string& new_tag);

  bool fileExists(const std::string& file_name)
  {
    std::ifstream f(file_name.c_str());
    return f.good();
  }

  QDP_serialparallel_t readParallelIO(XMLReader& paramtop)
  {
    bool parioP = Layout::isIOGridDefined() && (Layout::numIONodeGrid() > 1);

    if (paramtop.count("./parallel_io") > 0)
    {
      read(paramtop, "./parallel_io", parioP);
    }
    else if (paramtop.count("./ParallelIO") > 0)
    {
      read(paramtop, "./ParallelIO", parioP);
    }

    if (parioP)
    {
      QDPIO::cout << "Setting parallel write mode for saving configurations" << std::endl;
      return QDPIO_PARALLEL;
    }

    QDPIO::cout << "Setting serial write mode for saving configurations" << std::endl;
    return QDPIO_SERIAL;
  }

  struct ChainControl
  {
    QDP::Seed rng_seed;
    unsigned long start_update_num;
    unsigned long n_warm_up_updates;
    unsigned long n_production_updates;
    unsigned int n_updates_this_run;
    unsigned int save_interval;
    std::string save_prefix;
    QDP_volfmt_t save_volfmt;
    QDP_serialparallel_t save_pario;
    bool repro_checkP;
    int repro_check_frequency;
    bool rev_checkP;
    int rev_check_frequency;
    bool monitorForcesP;
  };

  void read(XMLReader& xml, const std::string& path, ChainControl& p)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "./RNG", p.rng_seed);
    read(paramtop, "./StartUpdateNum", p.start_update_num);
    read(paramtop, "./NWarmUpUpdates", p.n_warm_up_updates);
    read(paramtop, "./NProductionUpdates", p.n_production_updates);
    read(paramtop, "./NUpdatesThisRun", p.n_updates_this_run);
    read(paramtop, "./SaveInterval", p.save_interval);
    read(paramtop, "./SavePrefix", p.save_prefix);
    read(paramtop, "./SaveVolfmt", p.save_volfmt);
    p.save_pario = readParallelIO(paramtop);

    p.repro_checkP = true;
    p.repro_check_frequency = 10;
    if (paramtop.count("./ReproCheckP") == 1)
      read(paramtop, "./ReproCheckP", p.repro_checkP);
    if (p.repro_checkP && paramtop.count("./ReproCheckFrequency") == 1)
      read(paramtop, "./ReproCheckFrequency", p.repro_check_frequency);

    p.rev_checkP = true;
    p.rev_check_frequency = 10;
    if (paramtop.count("./ReverseCheckP") == 1)
      read(paramtop, "./ReverseCheckP", p.rev_checkP);
    if (p.rev_checkP && paramtop.count("./ReverseCheckFrequency") == 1)
      read(paramtop, "./ReverseCheckFrequency", p.rev_check_frequency);

    if (paramtop.count("./MonitorForces") == 1)
      read(paramtop, "./MonitorForces", p.monitorForcesP);
    else
      p.monitorForcesP = true;
  }

  void write(XMLWriter& xml, const std::string& path, const ChainControl& p)
  {
    push(xml, path);
    write(xml, "RNG", p.rng_seed);
    write(xml, "StartUpdateNum", p.start_update_num);
    write(xml, "NWarmUpUpdates", p.n_warm_up_updates);
    write(xml, "NProductionUpdates", p.n_production_updates);
    write(xml, "NUpdatesThisRun", p.n_updates_this_run);
    write(xml, "SaveInterval", p.save_interval);
    write(xml, "SavePrefix", p.save_prefix);
    write(xml, "SaveVolfmt", p.save_volfmt);
    {
      bool pario = (p.save_pario == QDPIO_PARALLEL);
      write(xml, "ParallelIO", pario);
    }
    write(xml, "ReproCheckP", p.repro_checkP);
    if (p.repro_checkP)
      write(xml, "ReproCheckFrequency", p.repro_check_frequency);
    write(xml, "ReverseCheckP", p.rev_checkP);
    if (p.rev_checkP)
      write(xml, "ReverseCheckFrequency", p.rev_check_frequency);
    write(xml, "MonitorForces", p.monitorForcesP);
    pop(xml);
  }

  struct HMCDefinition
  {
    std::string monomials_xml;
    std::string hamiltonian_xml;
    std::string integrator_xml;
  };

  void readHMCDefinition(XMLReader& xml, const std::string& path, HMCDefinition& p)
  {
    XMLReader top(xml, path);

    {
      XMLReader reader(top, "./Monomials");
      std::ostringstream os;
      reader.print(os);
      p.monomials_xml = os.str();
    }

    {
      XMLReader reader(top, "./Hamiltonian");
      std::ostringstream os;
      reader.print(os);
      p.hamiltonian_xml = os.str();
    }

    {
      XMLReader reader(top, "./MDIntegrator");
      std::ostringstream os;
      reader.print(os);
      p.integrator_xml = os.str();
    }
  }

  struct HMCOverrides
  {
    bool has_monomials;
    bool has_hamiltonian;
    bool has_integrator;
    std::string monomials_xml;
    std::string hamiltonian_xml;
    std::string integrator_xml;
  };

  void read(XMLReader& xml, const std::string& path, HMCOverrides& p)
  {
    XMLReader top(xml, path);

    p.has_monomials = false;
    p.has_hamiltonian = false;
    p.has_integrator = false;

    if (top.count("./Monomials") > 0)
    {
      XMLReader reader(top, "./Monomials");
      std::ostringstream os;
      reader.print(os);
      p.monomials_xml = os.str();
      p.has_monomials = true;
    }

    if (top.count("./Hamiltonian") > 0)
    {
      XMLReader reader(top, "./Hamiltonian");
      std::ostringstream os;
      reader.print(os);
      p.hamiltonian_xml = os.str();
      p.has_hamiltonian = true;
    }

    if (top.count("./MDIntegrator") > 0)
    {
      XMLReader reader(top, "./MDIntegrator");
      std::ostringstream os;
      reader.print(os);
      p.integrator_xml = os.str();
      p.has_integrator = true;
    }
  }

  struct SharedHMC
  {
    int loop_extent;
    HMCDefinition definition;
  };

  void read(XMLReader& xml, const std::string& path, SharedHMC& p)
  {
    XMLReader top(xml, path);
    std::string bc_name;
    read(top, "./Boundary/GaugeBC", bc_name);
    if (bc_name != SubdomainFixedGaugeBCEnv::name)
    {
      QDPIO::cerr << "SubdomainHMC: SharedHMC/Boundary/GaugeBC must be "
                  << SubdomainFixedGaugeBCEnv::name << std::endl;
      QDP_abort(1);
    }

    p.loop_extent = 1;
    if (top.count("./Boundary/loop_extent") > 0)
      read(top, "./Boundary/loop_extent", p.loop_extent);

    readHMCDefinition(top, ".", p.definition);
  }

  struct ChildConfig
  {
    ChainControl chain;
    bool has_overrides;
    HMCOverrides overrides;
    std::string inline_measurement_xml;
  };

  void readInlineMeasurementXML(XMLReader& top,
                                const std::string& path,
                                std::string& xml_out_string)
  {
    if (top.count(path) == 0)
    {
      XMLBufferWriter dummy;
      push(dummy, "InlineMeasurements");
      pop(dummy);
      xml_out_string = dummy.printCurrentContext();
      return;
    }

    XMLReader reader(top, path);
    std::ostringstream os;
    reader.print(os);
    xml_out_string = os.str();
  }

  void read(XMLReader& xml, const std::string& path, ChildConfig& p)
  {
    XMLReader top(xml, path);
    read(top, "./ChainControl", p.chain);
    p.has_overrides = (top.count("./HMCOverrides") > 0);
    if (p.has_overrides)
      read(top, "./HMCOverrides", p.overrides);
    readInlineMeasurementXML(top, "./InlineMeasurements", p.inline_measurement_xml);
  }

  enum StartModeType
  {
    START_FROM_PARENT_CFG,
    START_FROM_CHILD_RESTART
  };

  struct StartConfig
  {
    StartModeType mode;
    GroupXML_t cfg;
    GaugeSubdomainSplitParams split;
    std::string child_restart_manifest;
  };

  void read(XMLReader& xml, const std::string& path, StartConfig& p)
  {
    XMLReader top(xml, path);
    std::string mode;
    read(top, "./Mode", mode);

    if (mode == "FROM_PARENT_CFG")
    {
      p.mode = START_FROM_PARENT_CFG;
      p.cfg = readXMLGroup(top, "Cfg", "cfg_type");
      read(top, "./Split", p.split);
    }
    else if (mode == "FROM_CHILD_RESTART")
    {
      p.mode = START_FROM_CHILD_RESTART;
      read(top, "./ChildRestartManifest", p.child_restart_manifest);
    }
    else
    {
      QDPIO::cerr << "SubdomainHMC: unknown Start/Mode = " << mode << std::endl;
      QDP_abort(1);
    }
  }

  struct OutputConfig
  {
    bool stitch_at_end;
    bool write_stitched_parent;
    bool keep_child_restart_artifacts;
    std::string child_restart_manifest_file;
    std::string stitched_parent_file;
    QDP_volfmt_t stitched_parent_volfmt;
    QDP_serialparallel_t stitched_parent_pario;
    std::string post_stitch_inline_measurement_xml;
  };

  void read(XMLReader& xml, const std::string& path, OutputConfig& p)
  {
    XMLReader top(xml, path);

    p.stitch_at_end = false;
    if (top.count("./StitchAtEnd") > 0)
      read(top, "./StitchAtEnd", p.stitch_at_end);

    p.write_stitched_parent = false;
    if (top.count("./WriteStitchedParent") > 0)
      read(top, "./WriteStitchedParent", p.write_stitched_parent);

    p.keep_child_restart_artifacts = true;
    if (top.count("./KeepChildRestartArtifacts") > 0)
      read(top, "./KeepChildRestartArtifacts", p.keep_child_restart_artifacts);

    read(top, "./ChildRestartManifestFile", p.child_restart_manifest_file);

    p.stitched_parent_file = "";
    p.stitched_parent_volfmt = QDPIO_SINGLEFILE;
    p.stitched_parent_pario = readParallelIO(top);
    if (top.count("./StitchedParentFile") > 0)
      read(top, "./StitchedParentFile", p.stitched_parent_file);
    if (top.count("./StitchedParentSaveVolfmt") > 0)
      read(top, "./StitchedParentSaveVolfmt", p.stitched_parent_volfmt);
    if (top.count("./StitchedParentParallelIO") > 0)
    {
      bool pario = false;
      read(top, "./StitchedParentParallelIO", pario);
      p.stitched_parent_pario = pario ? QDPIO_PARALLEL : QDPIO_SERIAL;
    }

    readInlineMeasurementXML(top,
                             "./PostStitchInlineMeasurements",
                             p.post_stitch_inline_measurement_xml);
  }

  struct SubdomainHMCParams
  {
    multi1d<int> nrow;
    StartConfig start;
    SharedHMC shared_hmc;
    ChildConfig child0;
    ChildConfig child1;
    OutputConfig outputs;
  };

  void read(XMLReader& xml, const std::string& path, SubdomainHMCParams& p)
  {
    XMLReader top(xml, path);
    p.nrow.resize(0);
    if (top.count("./nrow") > 0)
      read(top, "./nrow", p.nrow);
    read(top, "./Start", p.start);
    read(top, "./SharedHMC", p.shared_hmc);
    read(top, "./Child0", p.child0);
    read(top, "./Child1", p.child1);
    read(top, "./Outputs", p.outputs);
  }

  struct ChildRestartEntry
  {
    multi1d<int> nrow;
    multi1d<int> local_to_global_t;
    multi1d<GaugeSubdomainSplitInterval> frozen_local_intervals;
    QDP::Seed rng_seed;
    unsigned long start_update_num;
    GroupXML_t cfg;
    GroupXML_t mom;
  };

  void write(XMLWriter& xml, const std::string& path, const ChildRestartEntry& p)
  {
    push(xml, path);
    write(xml, "nrow", p.nrow);
    write(xml, "local_to_global_t", p.local_to_global_t);
    write(xml, "frozen_local_intervals", p.frozen_local_intervals);
    write(xml, "RNG", p.rng_seed);
    write(xml, "StartUpdateNum", p.start_update_num);
    xml << p.cfg.xml;
    xml << renameOuterTag(p.mom.xml, "Cfg", "Momenta");
    pop(xml);
  }

  void read(XMLReader& xml, const std::string& path, ChildRestartEntry& p)
  {
    XMLReader top(xml, path);
    read(top, "./nrow", p.nrow);
    read(top, "./local_to_global_t", p.local_to_global_t);
    read(top, "./frozen_local_intervals", p.frozen_local_intervals);
    read(top, "./RNG", p.rng_seed);
    read(top, "./StartUpdateNum", p.start_update_num);
    p.cfg = readXMLGroup(top, "Cfg", "cfg_type");
    p.mom = readXMLGroup(top, "Momenta", "cfg_type");
  }

  struct ChildRestartManifest
  {
    int version;
    GaugeSubdomainSplitParams split;
    multi1d<int> parent_nrow;
    ChildRestartEntry child0;
    ChildRestartEntry child1;
  };

  void write(XMLWriter& xml, const std::string& path, const ChildRestartManifest& p)
  {
    push(xml, path);
    write(xml, "Version", p.version);
    write(xml, "Split", p.split);
    write(xml, "parent_nrow", p.parent_nrow);
    write(xml, "Child0", p.child0);
    write(xml, "Child1", p.child1);
    pop(xml);
  }

  void read(XMLReader& xml, const std::string& path, ChildRestartManifest& p)
  {
    XMLReader top(xml, path);
    read(top, "./Version", p.version);
    read(top, "./Split", p.split);
    read(top, "./parent_nrow", p.parent_nrow);
    read(top, "./Child0", p.child0);
    read(top, "./Child1", p.child1);
  }

  struct ChildChainRuntime
  {
    GaugeSubdomainSplitChild split_child;
    multi1d<LatticeColorMatrix> p;
    QDP::Seed rng_seed;
    unsigned long start_update_num;
  };

  struct RunState
  {
    GaugeSubdomainSplitParams split;
    multi1d<int> parent_nrow;
    ChildChainRuntime child0;
    ChildChainRuntime child1;
  };

  HMCDefinition effectiveHMCDefinition(const SharedHMC& shared, const ChildConfig& child)
  {
    HMCDefinition out = shared.definition;
    if (child.has_overrides)
    {
      if (child.overrides.has_monomials)
        out.monomials_xml = child.overrides.monomials_xml;
      if (child.overrides.has_hamiltonian)
        out.hamiltonian_xml = child.overrides.hamiltonian_xml;
      if (child.overrides.has_integrator)
        out.integrator_xml = child.overrides.integrator_xml;
    }
    return out;
  }

  std::string subdomainGaugeBCXML(const int loop_extent)
  {
    std::ostringstream os;
    os << "<GaugeBC>"
       << "<Name>" << SubdomainFixedGaugeBCEnv::name << "</Name>"
       << "<loop_extent>" << loop_extent << "</loop_extent>"
       << "</GaugeBC>";
    return os.str();
  }

  std::string injectSubdomainGaugeBC(const std::string& monomials_xml, const int loop_extent)
  {
    std::string output;
    std::string remaining = monomials_xml;
    const std::string open_tag = "<GaugeAction>";
    const std::string close_tag = "</GaugeAction>";
    const std::string bc_open = "<GaugeBC>";
    const std::string bc_close = "</GaugeBC>";
    const std::string bc_xml = subdomainGaugeBCXML(loop_extent);

    std::string::size_type cursor = 0;
    while (true)
    {
      const std::string::size_type open = remaining.find(open_tag, cursor);
      if (open == std::string::npos)
      {
        output += remaining.substr(cursor);
        break;
      }

      const std::string::size_type close = remaining.find(close_tag, open);
      if (close == std::string::npos)
      {
        QDPIO::cerr << "SubdomainHMC: malformed Monomials XML while injecting GaugeBC" << std::endl;
        QDP_abort(1);
      }

      output += remaining.substr(cursor, open - cursor);

      std::string gauge_action_block =
        remaining.substr(open, close + close_tag.size() - open);

      const std::string::size_type bc_begin = gauge_action_block.find(bc_open);
      if (bc_begin != std::string::npos)
      {
        const std::string::size_type bc_end = gauge_action_block.find(bc_close, bc_begin);
        if (bc_end == std::string::npos)
        {
          QDPIO::cerr << "SubdomainHMC: malformed GaugeBC XML in Monomials block" << std::endl;
          QDP_abort(1);
        }
        gauge_action_block.erase(bc_begin, bc_end + bc_close.size() - bc_begin);
      }

      gauge_action_block.insert(gauge_action_block.size() - close_tag.size(), bc_xml);
      output += gauge_action_block;
      cursor = close + close_tag.size();
    }

    return output;
  }

  multi1d<std::string> collectMonomialIds(const std::string& monomials_xml)
  {
    std::istringstream is(monomials_xml);
    XMLReader reader(is);
    XMLReader top(reader, "/Monomials");
    const int n_items = top.count("./elem");
    multi1d<std::string> ids(n_items);
    for (int i = 0; i < n_items; ++i)
    {
      std::ostringstream path;
      path << "./elem[" << i + 1 << "]";
      XMLReader elem(top, path.str());
      std::string monomial_name;
      read(elem, "./Name", monomial_name);
      if (monomial_name != GaugeMonomialEnv::name)
      {
        QDPIO::cerr << "SubdomainHMC: only GAUGE_MONOMIAL is supported in the first version" << std::endl;
        QDP_abort(1);
      }
      read(elem, "./NamedObject/monomial_id", ids[i]);
    }
    return ids;
  }

  void eraseMonomials(const multi1d<std::string>& monomial_ids)
  {
    for (int i = 0; i < monomial_ids.size(); ++i)
    {
      if (TheNamedObjMap::Instance().check(monomial_ids[i]))
        TheNamedObjMap::Instance().erase(monomial_ids[i]);
    }
  }

  GroupXML_t restartGroupXML(const std::string& file_name, const QDP_serialparallel_t pario)
  {
    SZINQIOGaugeInitEnv::Params cfg;
    cfg.cfg_file = file_name;
    cfg.cfg_pario = pario;
    cfg.reunitP = false;
    return SZINQIOGaugeInitEnv::createXMLGroup(cfg);
  }

  std::string renameOuterTag(const std::string& xml_string,
                             const std::string& old_tag,
                             const std::string& new_tag)
  {
    std::string renamed = xml_string;
    const std::string old_open = "<" + old_tag + ">";
    const std::string old_close = "</" + old_tag + ">";
    const std::string new_open = "<" + new_tag + ">";
    const std::string new_close = "</" + new_tag + ">";

    const std::string::size_type open_pos = renamed.find(old_open);
    if (open_pos == std::string::npos)
    {
      QDPIO::cerr << "SubdomainHMC: could not find opening tag " << old_open
                  << " while writing restart metadata" << std::endl;
      QDP_abort(1);
    }
    renamed.replace(open_pos, old_open.size(), new_open);

    const std::string::size_type close_pos = renamed.rfind(old_close);
    if (close_pos == std::string::npos)
    {
      QDPIO::cerr << "SubdomainHMC: could not find closing tag " << old_close
                  << " while writing restart metadata" << std::endl;
      QDP_abort(1);
    }
    renamed.replace(close_pos, old_close.size(), new_close);

    return renamed;
  }

  void loadGaugeLikeField(const GroupXML_t& cfg,
                          multi1d<LatticeColorMatrix>& u,
                          const bool force_no_reunit)
  {
    XMLReader file_xml;
    XMLReader config_xml;
    GroupXML_t load_cfg = cfg;

    if (cfg.id == SZINQIOGaugeInitEnv::name || cfg.id == "SCIDAC")
    {
      std::istringstream xml_c(cfg.xml);
      XMLReader cfgtop(xml_c);
      SZINQIOGaugeInitEnv::Params params(cfgtop, cfg.path);
      if (force_no_reunit)
        params.reunitP = false;
      load_cfg = SZINQIOGaugeInitEnv::createXMLGroup(params);
    }

    std::istringstream xml_c(load_cfg.xml);
    XMLReader cfgtop(xml_c);
    Handle<GaugeInit> gaugeInit(
      TheGaugeInitFactory::Instance().createObject(load_cfg.id, cfgtop, load_cfg.path));
    (*gaugeInit)(file_xml, config_xml, u);
  }

  void saveGaugeLikeField(const std::string& file_name,
                          const multi1d<LatticeColorMatrix>& u,
                          const QDP_volfmt_t volfmt,
                          const QDP_serialparallel_t pario,
                          const std::string& root_tag)
  {
    XMLBufferWriter file_xml;
    push(file_xml, root_tag);
    proginfo(file_xml);
    pop(file_xml);

    XMLBufferWriter record_xml;
    push(record_xml, "Record");
    pop(record_xml);

    writeGauge(file_xml, record_xml, u, file_name, volfmt, pario);
  }

  ChildRestartManifest currentRestartManifest(const RunState& state,
                                              const ChildConfig& child0_cfg,
                                              const ChildConfig& child1_cfg)
  {
    ChildRestartManifest manifest;
    manifest.version = 1;
    manifest.split = state.split;
    manifest.parent_nrow = state.parent_nrow;

    manifest.child0.nrow = state.child0.split_child.nrow;
    manifest.child0.local_to_global_t = state.child0.split_child.local_to_global_t;
    manifest.child0.frozen_local_intervals = state.child0.split_child.frozen_local_intervals;
    manifest.child0.rng_seed = state.child0.rng_seed;
    manifest.child0.start_update_num = state.child0.start_update_num;
    manifest.child0.cfg = restartGroupXML(child0_cfg.chain.save_prefix + "_cfg.lime",
                                          child0_cfg.chain.save_pario);
    manifest.child0.mom = restartGroupXML(child0_cfg.chain.save_prefix + "_mom.lime",
                                          child0_cfg.chain.save_pario);

    manifest.child1.nrow = state.child1.split_child.nrow;
    manifest.child1.local_to_global_t = state.child1.split_child.local_to_global_t;
    manifest.child1.frozen_local_intervals = state.child1.split_child.frozen_local_intervals;
    manifest.child1.rng_seed = state.child1.rng_seed;
    manifest.child1.start_update_num = state.child1.start_update_num;
    manifest.child1.cfg = restartGroupXML(child1_cfg.chain.save_prefix + "_cfg.lime",
                                          child1_cfg.chain.save_pario);
    manifest.child1.mom = restartGroupXML(child1_cfg.chain.save_prefix + "_mom.lime",
                                          child1_cfg.chain.save_pario);

    return manifest;
  }

  void saveCombinedChildRestart(const RunState& state,
                                const ChildConfig& child0_cfg,
                                const ChildConfig& child1_cfg,
                                const OutputConfig& outputs)
  {
    LayoutGuard guard;

    LayoutGuard::activate(state.child0.split_child.nrow);
    saveGaugeLikeField(child0_cfg.chain.save_prefix + "_cfg.lime",
                       state.child0.split_child.u,
                       child0_cfg.chain.save_volfmt,
                       child0_cfg.chain.save_pario,
                       "SubdomainHMCChild0Cfg");
    saveGaugeLikeField(child0_cfg.chain.save_prefix + "_mom.lime",
                       state.child0.p,
                       child0_cfg.chain.save_volfmt,
                       child0_cfg.chain.save_pario,
                       "SubdomainHMCChild0Mom");

    LayoutGuard::activate(state.child1.split_child.nrow);
    saveGaugeLikeField(child1_cfg.chain.save_prefix + "_cfg.lime",
                       state.child1.split_child.u,
                       child1_cfg.chain.save_volfmt,
                       child1_cfg.chain.save_pario,
                       "SubdomainHMCChild1Cfg");
    saveGaugeLikeField(child1_cfg.chain.save_prefix + "_mom.lime",
                       state.child1.p,
                       child1_cfg.chain.save_volfmt,
                       child1_cfg.chain.save_pario,
                       "SubdomainHMCChild1Mom");

    XMLFileWriter manifest_xml(outputs.child_restart_manifest_file);
    write(manifest_xml,
          "SubdomainChildRestart",
          currentRestartManifest(state, child0_cfg, child1_cfg));
    manifest_xml.close();
  }

  void loadChildRestartManifest(const std::string& file_name, ChildRestartManifest& manifest)
  {
    XMLReader xml(file_name);
    read(xml, "/SubdomainChildRestart", manifest);
  }

  void restoreChildRuntime(const ChildRestartEntry& entry, ChildChainRuntime& runtime)
  {
    LayoutGuard::activate(entry.nrow);
    runtime.split_child.nrow = entry.nrow;
    runtime.split_child.local_to_global_t = entry.local_to_global_t;
    runtime.split_child.frozen_local_intervals = entry.frozen_local_intervals;
    loadGaugeLikeField(entry.cfg, runtime.split_child.u, false);
    loadGaugeLikeField(entry.mom, runtime.p, true);
    runtime.rng_seed = entry.rng_seed;
    runtime.start_update_num = entry.start_update_num;
  }

  void clearChildRuntime(ChildChainRuntime& runtime)
  {
    if (runtime.split_child.nrow.size() != 0)
      LayoutGuard::activate(runtime.split_child.nrow);

    runtime.split_child.u.resize(0);
    runtime.p.resize(0);
    runtime.split_child.nrow.resize(0);
    runtime.split_child.local_to_global_t.resize(0);
    runtime.split_child.frozen_local_intervals.resize(0);
    runtime.start_update_num = 0;
  }

  void loadInlineMeasurements(const std::string& xml_string,
                              multi1d<Handle<AbsInlineMeasurement> >& measurements)
  {
    std::istringstream is(xml_string);
    XMLReader xml(is);
    read(xml, "/InlineMeasurements", measurements);
  }

  bool checkReproducability(const multi1d<LatticeColorMatrix>& P_new,
                            const multi1d<LatticeColorMatrix>& Q_new,
                            const QDP::Seed& seed_new,
                            const multi1d<LatticeColorMatrix>& P_old,
                            const multi1d<LatticeColorMatrix>& Q_old,
                            const QDP::Seed& seed_old)
  {
#if ! defined(QDP_IS_QDPJIT2)
    int diffs_found = 0;
    if (P_new.size() != P_old.size() || Q_new.size() != Q_old.size())
      return false;

    const int bytes =
      2 * Nc * Nc * Layout::sitesOnNode() * sizeof(WordType<LatticeColorMatrix>::Type_t);

    for (int mu = 0; mu < P_new.size(); ++mu)
    {
      const unsigned char* p1 = (const unsigned char*)P_new[mu].getF();
      const unsigned char* p2 = (const unsigned char*)P_old[mu].getF();
      for (int b = 0; b < bytes; ++b)
      {
        const unsigned char diff = *p1 - *p2;
        if (diff != 0)
          diffs_found++;
        ++p1;
        ++p2;
      }
    }

    QDPInternal::globalSum(diffs_found);
    if (diffs_found != 0)
      return false;

    diffs_found = 0;
    for (int mu = 0; mu < Q_new.size(); ++mu)
    {
      const unsigned char* p1 = (const unsigned char*)Q_new[mu].getF();
      const unsigned char* p2 = (const unsigned char*)Q_old[mu].getF();
      for (int b = 0; b < bytes; ++b)
      {
        const unsigned char diff = *p1 - *p2;
        if (diff != 0)
          diffs_found++;
        ++p1;
        ++p2;
      }
    }

    QDPInternal::globalSum(diffs_found);
    if (diffs_found != 0)
      return false;

    if (!toBool(seed_new == seed_old))
      return false;
#else
    QDPIO::cout << "qdp-jit2: skipping momentum repro check" << std::endl;
#endif

    return true;
  }

  void enforceFrozenChildState(GaugeFieldState& gauge_state, const SubdomainFixedGaugeBC& bc)
  {
    bc.clampToFrozen(gauge_state.getQ());
    bc.zero(gauge_state.getP());
  }

  void validateFrozenChildState(const std::string& label,
                                const multi1d<LatticeColorMatrix>& q,
                                const SubdomainFixedGaugeBC& bc)
  {
    multi1d<LatticeColorMatrix> restored = q;
    bc.clampToFrozen(restored);

    Double diff = zero;
    for (int mu = 0; mu < Nd; ++mu)
      diff += norm2(q[mu] - restored[mu]);

    if (toDouble(diff) != 0.0)
    {
      QDPIO::cerr << "SubdomainHMC: frozen-link drift detected in " << label
                  << " after enforcing the child boundary" << std::endl;
      QDP_abort(1);
    }
  }

  void doChildMeasurements(const std::string& label,
                           const unsigned long update_no,
                           const multi1d<LatticeColorMatrix>& q,
                           const SubdomainFixedGaugeBC& bc,
                           multi1d<Handle<AbsInlineMeasurement> >& default_measurements,
                           multi1d<Handle<AbsInlineMeasurement> >& user_measurements,
                           XMLWriter& xml_out)
  {
    XMLBufferWriter gauge_xml;
    push(gauge_xml, "SubdomainHMCChild");
    write(gauge_xml, "label", label);
    write(gauge_xml, "update_no", update_no);
    pop(gauge_xml);

    multi1d<LatticeColorMatrix> measurement_q = q;
    bc.modify(measurement_q);

    InlineDefaultGaugeField::reset();
    InlineDefaultGaugeField::set(measurement_q, gauge_xml);

    push(xml_out, "InlineObservables");

    for (int m = 0; m < default_measurements.size(); ++m)
    {
      push(xml_out, "elem");
      (*(default_measurements[m]))(update_no, xml_out);
      pop(xml_out);
    }

    for (int m = 0; m < user_measurements.size(); ++m)
    {
      AbsInlineMeasurement& meas = *(user_measurements[m]);
      if (update_no % meas.getFrequency() == 0)
      {
        push(xml_out, "elem");
        meas(update_no, xml_out);
        pop(xml_out);
      }
    }

    pop(xml_out);
    InlineDefaultGaugeField::reset();
  }

  void runChildChain(const std::string& label,
                     ChildChainRuntime& runtime,
                     const ChildConfig& child_cfg,
                     const ChildConfig& child0_cfg,
                     const ChildConfig& child1_cfg,
                     const SharedHMC& shared_hmc,
                     const HMCDefinition& hmc_definition,
                     RunState& state,
                     const OutputConfig& outputs)
  {
    LayoutGuard::activate(runtime.split_child.nrow);

    Handle<SubdomainFixedGaugeBC> bc(new SubdomainFixedGaugeBC(runtime.split_child.u,
                                                               runtime.split_child.frozen_local_intervals,
                                                               state.split.t_dir,
                                                               shared_hmc.loop_extent));
    SubdomainFixedGaugeBCEnv::setRuntimePrototype(bc);

    const std::string monomials_xml =
      injectSubdomainGaugeBC(hmc_definition.monomials_xml, shared_hmc.loop_extent);
    const multi1d<std::string> monomial_ids = collectMonomialIds(monomials_xml);

    Handle<AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > H_MC;
    Handle<AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > integrator;
    multi1d<Handle<AbsInlineMeasurement> > measurements;

    try
    {
      std::istringstream monomial_is(monomials_xml);
      XMLReader monomial_reader(monomial_is);
      readNamedMonomialArray(monomial_reader, "/Monomials");

      std::istringstream ham_is(hmc_definition.hamiltonian_xml);
      XMLReader ham_xml(ham_is);
      ExactHamiltonianParams ham_params(ham_xml, "/Hamiltonian");
      H_MC = new ExactHamiltonian(ham_params);

      std::istringstream int_is(hmc_definition.integrator_xml);
      XMLReader int_xml(int_is);
      LCMToplevelIntegratorParams int_par(int_xml, "/MDIntegrator");
      integrator = new LCMToplevelIntegrator(int_par);

      loadInlineMeasurements(child_cfg.inline_measurement_xml, measurements);
    }
    catch (const std::string& e)
    {
      eraseMonomials(monomial_ids);
      SubdomainFixedGaugeBCEnv::clearRuntimePrototype();
      QDPIO::cerr << "SubdomainHMC: failed to construct child chain objects: " << e << std::endl;
      QDP_abort(1);
    }

    SubdomainLatColMatHMCTrj hmc_trj(H_MC, integrator, *bc);

    multi1d<Handle<AbsInlineMeasurement> > default_measurements(1);
    InlinePlaquetteEnv::Params plaq_params;
    plaq_params.frequency = 1;
    default_measurements[0] = new InlinePlaquetteEnv::InlineMeas(plaq_params);

    XMLWriter& xml_out = TheXMLOutputWriter::Instance();
    XMLWriter& xml_log = TheXMLLogWriter::Instance();

    push(xml_out, label);
    push(xml_log, label);

    setForceMonitoring(child_cfg.chain.monitorForcesP);

    QDP::RNG::setrn(runtime.rng_seed);
    GaugeFieldState gauge_state(runtime.p, runtime.split_child.u);

    unsigned long cur_update = runtime.start_update_num;
    const unsigned long total_updates =
      child_cfg.chain.n_warm_up_updates + child_cfg.chain.n_production_updates;
    const unsigned long remaining =
      (total_updates > cur_update) ? (total_updates - cur_update) : 0;
    const unsigned long to_do =
      std::min<unsigned long>(remaining, child_cfg.chain.n_updates_this_run);

    push(xml_out, "MCUpdates");
    push(xml_log, "MCUpdates");

    for (unsigned long i = 0; i < to_do; ++i)
    {
      push(xml_out, "elem");
      push(xml_log, "elem");
      push(xml_out, "Update");
      push(xml_log, "Update");

      ++cur_update;
      const bool warm_up_p = cur_update <= child_cfg.chain.n_warm_up_updates;
      const bool do_reverse =
        child_cfg.chain.rev_checkP &&
        (cur_update % child_cfg.chain.rev_check_frequency == 0);

      write(xml_out, "update_no", cur_update);
      write(xml_log, "update_no", cur_update);
      write(xml_out, "WarmUpP", warm_up_p);
      write(xml_log, "WarmUpP", warm_up_p);

      QDP::StopWatch swatch;

      if (child_cfg.chain.repro_checkP &&
          (cur_update % child_cfg.chain.repro_check_frequency == 0))
      {
        GaugeFieldState repro_bkup_start(gauge_state.getP(), gauge_state.getQ());
        QDP::Seed rng_seed_bkup_start;
        QDP::RNG::savern(rng_seed_bkup_start);

        swatch.start();
        hmc_trj(gauge_state, warm_up_p, do_reverse);
        enforceFrozenChildState(gauge_state, *bc);
        swatch.stop();
        write(xml_out, "seconds_for_trajectory", swatch.getTimeInSeconds());
        write(xml_log, "seconds_for_trajectory", swatch.getTimeInSeconds());

        GaugeFieldState repro_bkup_end(gauge_state.getP(), gauge_state.getQ());
        QDP::Seed rng_seed_bkup_end;
        QDP::RNG::savern(rng_seed_bkup_end);

        gauge_state.getP() = repro_bkup_start.getP();
        gauge_state.getQ() = repro_bkup_start.getQ();
        QDP::RNG::setrn(rng_seed_bkup_start);

        swatch.reset();
        swatch.start();
        hmc_trj(gauge_state, warm_up_p, false);
        enforceFrozenChildState(gauge_state, *bc);
        swatch.stop();
        write(xml_out, "seconds_for_repro_trajectory", swatch.getTimeInSeconds());
        write(xml_log, "seconds_for_repro_trajectory", swatch.getTimeInSeconds());

        QDP::Seed rng_seed_end2;
        QDP::RNG::savern(rng_seed_end2);

        const bool pass =
          checkReproducability(gauge_state.getP(),
                               gauge_state.getQ(),
                               rng_seed_end2,
                               repro_bkup_end.getP(),
                               repro_bkup_end.getQ(),
                               rng_seed_bkup_end);
        write(xml_out, "ReproCheck", pass);
        write(xml_log, "ReproCheck", pass);
        if (!pass)
        {
          eraseMonomials(monomial_ids);
          SubdomainFixedGaugeBCEnv::clearRuntimePrototype();
          QDP_abort(1);
        }
      }
      else
      {
        swatch.start();
        hmc_trj(gauge_state, warm_up_p, do_reverse);
        enforceFrozenChildState(gauge_state, *bc);
        swatch.stop();
        write(xml_out, "seconds_for_trajectory", swatch.getTimeInSeconds());
        write(xml_log, "seconds_for_trajectory", swatch.getTimeInSeconds());
      }

      doChildMeasurements(label,
                          cur_update,
                          gauge_state.getQ(),
                          *bc,
                          default_measurements,
                          measurements,
                          xml_out);

      runtime.p = gauge_state.getP();
      runtime.split_child.u = gauge_state.getQ();
      validateFrozenChildState(label, runtime.split_child.u, *bc);
      runtime.start_update_num = cur_update;
      QDP::RNG::savern(runtime.rng_seed);

      if (cur_update % child_cfg.chain.save_interval == 0)
        saveCombinedChildRestart(state, child0_cfg, child1_cfg, outputs);

      pop(xml_log);
      pop(xml_out);
      pop(xml_log);
      pop(xml_out);
    }

    pop(xml_log);
    pop(xml_out);
    pop(xml_log);
    pop(xml_out);

    eraseMonomials(monomial_ids);
    SubdomainFixedGaugeBCEnv::clearRuntimePrototype();
  }

  void runParentMeasurements(const std::string& xml_string,
                             const unsigned long update_no,
                             const multi1d<LatticeColorMatrix>& q)
  {
    multi1d<Handle<AbsInlineMeasurement> > measurements;
    multi1d<LatticeColorMatrix> measurement_q = q;
    const std::string parent_gauge_id = "subdomain_hmc_parent_gauge";
    bool parent_gauge_field_set = false;

    try
    {
      std::string measurement_xml = xml_string;
      const std::string default_gauge_id = InlineDefaultGaugeField::getId();
      std::string::size_type pos = 0;
      while ((pos = measurement_xml.find(default_gauge_id, pos)) != std::string::npos)
      {
        measurement_xml.replace(pos, default_gauge_id.size(), parent_gauge_id);
        pos += parent_gauge_id.size();
      }

      loadInlineMeasurements(measurement_xml, measurements);

      XMLBufferWriter gauge_xml;
      push(gauge_xml, "SubdomainHMCParent");
      write(gauge_xml, "update_no", update_no);
      pop(gauge_xml);

      if (TheNamedObjMap::Instance().check(parent_gauge_id))
        TheNamedObjMap::Instance().erase(parent_gauge_id);

      XMLBufferWriter file_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);

      TheNamedObjMap::Instance().create<multi1d<LatticeColorMatrix> >(parent_gauge_id);
      TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(parent_gauge_id) = measurement_q;
      TheNamedObjMap::Instance().get(parent_gauge_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(parent_gauge_id).setRecordXML(gauge_xml);
      parent_gauge_field_set = true;

      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "PostStitchInlineMeasurements");
      for (int m = 0; m < measurements.size(); ++m)
      {
        AbsInlineMeasurement& meas = *(measurements[m]);
        if (update_no % meas.getFrequency() == 0)
        {
          push(xml_out, "elem");
          meas(update_no, xml_out);
          pop(xml_out);
        }
      }
      pop(xml_out);

      TheNamedObjMap::Instance().erase(parent_gauge_id);
      parent_gauge_field_set = false;
    }
    catch (const std::string& e)
    {
      if (parent_gauge_field_set && TheNamedObjMap::Instance().check(parent_gauge_id))
        TheNamedObjMap::Instance().erase(parent_gauge_id);

      QDPIO::cerr << "SubdomainHMC: post-stitch measurement failure: " << e << std::endl;
      QDP_abort(1);
    }
  }

  bool linkageHack(void)
  {
    bool ok = true;
    ok &= GaugeMonomialEnv::registerAll();
    ok &= LCMMDComponentIntegratorAggregateEnv::registerAll();
    ok &= ChronoPredictorAggregrateEnv::registerAll();
    ok &= InlineAggregateEnv::registerAll();
    ok &= GaugeInitEnv::registerAll();
    return ok;
  }
}

int main(int argc, char* argv[])
{
  Chroma::initialize(&argc, &argv);

  START_CODE();

  {
    QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

    XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
    XMLFileWriter& xml_log = Chroma::getXMLLogInstance();

    push(xml_out, "subdomain_hmc");
    push(xml_log, "subdomain_hmc");

    SubdomainHMCParams params;
    try
    {
      XMLReader xml_in(Chroma::getXMLInputFileName());
      read(xml_in, "/Params/SubdomainHMC", params);
      write(xml_out, "Input", xml_in);
      write(xml_log, "Input", xml_in);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << "subdomain_hmc: failed to read input: " << e << std::endl;
      QDP_abort(1);
    }

    RunState state;
    XMLReader config_xml;

    if (params.start.mode == START_FROM_PARENT_CFG)
    {
      if (params.nrow.size() == 0)
      {
        QDPIO::cerr << "subdomain_hmc: nrow is required for FROM_PARENT_CFG startup" << std::endl;
        QDP_abort(1);
      }

      Layout::setLattSize(params.nrow);
      Layout::create();
      state.parent_nrow = params.nrow;
      state.split = params.start.split;

      multi1d<LatticeColorMatrix> parent_u(Nd);
      try
      {
        XMLReader file_xml;
        std::istringstream xml_c(params.start.cfg.xml);
        XMLReader cfgtop(xml_c);
        Handle<GaugeInit> gaugeInit(
          TheGaugeInitFactory::Instance().createObject(params.start.cfg.id, cfgtop, params.start.cfg.path));
        (*gaugeInit)(file_xml, config_xml, parent_u);
        write(xml_out, "Config_info", config_xml);
        write(xml_log, "Config_info", config_xml);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "subdomain_hmc: failed to load parent gauge field: " << e << std::endl;
        QDP_abort(1);
      }

      GaugeSubdomainSplitResult split = splitGaugeSubdomains(parent_u, state.split);
      state.child0.split_child = split.child0;
      state.child1.split_child = split.child1;
      state.child0.p.resize(Nd);
      state.child1.p.resize(Nd);
      state.child0.p = zero;
      state.child1.p = zero;
      state.child0.rng_seed = params.child0.chain.rng_seed;
      state.child1.rng_seed = params.child1.chain.rng_seed;
      state.child0.start_update_num = params.child0.chain.start_update_num;
      state.child1.start_update_num = params.child1.chain.start_update_num;
    }
    else
    {
      ChildRestartManifest manifest;
      loadChildRestartManifest(params.start.child_restart_manifest, manifest);
      if (manifest.version != 1)
      {
        QDPIO::cerr << "subdomain_hmc: unsupported child restart manifest version "
                    << manifest.version << std::endl;
        QDP_abort(1);
      }

      state.parent_nrow = manifest.parent_nrow;
      state.split = manifest.split;

      Layout::setLattSize(state.parent_nrow);
      Layout::create();

      restoreChildRuntime(manifest.child0, state.child0);
      restoreChildRuntime(manifest.child1, state.child1);
      Layout::setLattSize(state.parent_nrow);
      Layout::create();
    }

    const HMCDefinition child0_hmc = effectiveHMCDefinition(params.shared_hmc, params.child0);
    const HMCDefinition child1_hmc = effectiveHMCDefinition(params.shared_hmc, params.child1);

    runChildChain("Child0",
                  state.child0,
                  params.child0,
                  params.child0,
                  params.child1,
                  params.shared_hmc,
                  child0_hmc,
                  state,
                  params.outputs);

    runChildChain("Child1",
                  state.child1,
                  params.child1,
                  params.child0,
                  params.child1,
                  params.shared_hmc,
                  child1_hmc,
                  state,
                  params.outputs);

    saveCombinedChildRestart(state, params.child0, params.child1, params.outputs);

    if (params.outputs.stitch_at_end)
    {
      try
      {
        GaugeSubdomainSplitResult split;
        split.param = state.split;
        split.parent_nrow = state.parent_nrow;
        split.child0 = state.child0.split_child;
        split.child1 = state.child1.split_child;

        multi1d<LatticeColorMatrix> parent_u = stitchGaugeSubdomains(split);
        Layout::setLattSize(state.parent_nrow);
        Layout::create();

        const unsigned long update_no =
          std::max(state.child0.start_update_num, state.child1.start_update_num);

        clearChildRuntime(state.child0);
        clearChildRuntime(state.child1);
        Layout::setLattSize(state.parent_nrow);
        Layout::create();

        if (params.outputs.post_stitch_inline_measurement_xml.find("<elem>") != std::string::npos)
          runParentMeasurements(params.outputs.post_stitch_inline_measurement_xml, update_no, parent_u);

        if (params.outputs.write_stitched_parent)
        {
          saveGaugeLikeField(params.outputs.stitched_parent_file,
                             parent_u,
                             params.outputs.stitched_parent_volfmt,
                             params.outputs.stitched_parent_pario,
                             "SubdomainHMCStitchedParent");
        }
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << "subdomain_hmc: stitch/post-stitch failure: " << e << std::endl;
        QDP_abort(1);
      }
    }

    clearChildRuntime(state.child0);
    clearChildRuntime(state.child1);
    Layout::setLattSize(state.parent_nrow);
    Layout::create();

    pop(xml_log);
    pop(xml_out);
  }

  END_CODE();

  Chroma::finalize();
  return 0;
}
