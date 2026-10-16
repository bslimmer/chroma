/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "fermact.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "qdp.h"
#include "qdp_disk_map_slice.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_map_obj_memory.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/factorized_distillation.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_distillation.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/mgproton.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/superb_contractions.h"
#include "util/ferm/transf.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/gauge/gauge_subdomain_split.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <sstream>
#include <vector>

// Adding for Dirac operator handling in factorized measurement
#include "actions/ferm/fermstates/fermstates.h"
#include "handle.h"
#include "linearop.h"

#ifdef BUILD_SB

namespace Chroma {

//----------------------------------------------------------------------------
namespace InlinePropAndMatElemDistillationSuperbEnv {
//! Propagator input
void read(XMLReader &xml, const std::string &path,
          Params::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "gauge_id", input.gauge_id);
  read(inputtop, "colorvec_files", input.colorvec_files);
  read(inputtop, "prop_op_file", input.prop_op_file);

  input.frozen_boundary_sidecar_file.clear();
  if (inputtop.count("frozen_boundary_sidecar_file") == 1) {
    read(inputtop, "frozen_boundary_sidecar_file",
         input.frozen_boundary_sidecar_file);
  }

  input.frozen_boundary_child_id = -1;
  if (inputtop.count("frozen_boundary_child_id") == 1) {
    read(inputtop, "frozen_boundary_child_id", input.frozen_boundary_child_id);
  }

  input.factorized_intermediate_file.clear();
  if (inputtop.count("factorized_intermediate_file") == 1) {
    read(inputtop, "factorized_intermediate_file",
         input.factorized_intermediate_file);
  }

  input.factorized_cross_file.clear();
  if (inputtop.count("factorized_cross_file") == 1) {
    read(inputtop, "factorized_cross_file", input.factorized_cross_file);
  }

  input.num_intermediate_projectors = -1;
  if (inputtop.count("num_intermediate_projectors") == 1) {
    read(inputtop, "num_intermediate_projectors",
         input.num_intermediate_projectors);
  }

  input.intermediate_projector_phase.resize(Nd - 1);
  for (int i = 0; i < Nd - 1; ++i)
    input.intermediate_projector_phase[i] = 0;
  if (inputtop.count("intermediate_projector_phase") == 1) {
    read(inputtop, "intermediate_projector_phase",
         input.intermediate_projector_phase);
    if (input.intermediate_projector_phase.size() != Nd - 1) {
      throw std::runtime_error(
          std::string("intermediate_projector_phase should have ") +
          std::to_string(Nd - 1) + " components");
    }
  }

  input.overlapping_factorization = false;
  if (inputtop.count("overlapping_factorization") == 1) {
    read(inputtop, "overlapping_factorization",
         input.overlapping_factorization);
  }
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const Params::NamedObject_t &input) {
  push(xml, path);

  write(xml, "gauge_id", input.gauge_id);
  write(xml, "colorvec_files", input.colorvec_files);
  write(xml, "prop_op_file", input.prop_op_file);
  if (!input.frozen_boundary_sidecar_file.empty()) {
    write(xml, "frozen_boundary_sidecar_file",
          input.frozen_boundary_sidecar_file);
  }
  if (input.frozen_boundary_child_id >= 0) {
    write(xml, "frozen_boundary_child_id", input.frozen_boundary_child_id);
  }
  if (!input.factorized_intermediate_file.empty()) {
    write(xml, "factorized_intermediate_file",
          input.factorized_intermediate_file);
  }
  if (!input.factorized_cross_file.empty()) {
    write(xml, "factorized_cross_file", input.factorized_cross_file);
  }
  if (input.num_intermediate_projectors >= 0) {
    write(xml, "num_intermediate_projectors",
          input.num_intermediate_projectors);
  }

  write(xml, "intermediate_projector_phase",
        input.intermediate_projector_phase);
  if (input.overlapping_factorization) {
    write(xml, "overlapping_factorization", input.overlapping_factorization);
  }

  pop(xml);
}

//! Propagator input
void read(XMLReader &xml, const std::string &path,
          Params::Param_t::Phasing_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "sink", input.snk);
  if (input.snk.size() != Nd - 1) {
    throw std::runtime_error(std::string("phase sink should have ") +
                             std::to_string(Nd - 1) + " components");
  }
  read(inputtop, "source", input.src);
  if (input.src.size() != Nd - 1) {
    throw std::runtime_error(std::string("phase source should have ") +
                             std::to_string(Nd - 1) + " components");
  }
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const Params::Param_t::Phasing_t &input) {
  push(xml, path);

  write(xml, "sink", input.snk);
  write(xml, "source", input.src);

  pop(xml);
}

//! Propagator input
void read(XMLReader &xml, const std::string &path,
          Params::Param_t::Contract_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "num_vecs", input.num_vecs);
  read(inputtop, "t_sources", input.t_sources);
  read(inputtop, "decay_dir", input.decay_dir);
  read(inputtop, "Nt_forward", input.Nt_forward);
  read(inputtop, "Nt_backward", input.Nt_backward);
  read(inputtop, "mass_label", input.mass_label);

  input.max_rhs = 8;
  if (inputtop.count("max_rhs") == 1) {
    read(inputtop, "max_rhs", input.max_rhs);
  }

  if (inputtop.count("phases") == 1) {
    read(inputtop, "phases", input.phases);
  } else if (inputtop.count("phase") == 1) {
    read(inputtop, "phase", input.phase);
    if (input.phase.size() != Nd - 1) {
      throw std::runtime_error(std::string("phase tag should have ") +
                               std::to_string(Nd - 1) + " components");
    }
  } else {
    input.phase.resize(Nd - 1);
    for (int i = 0; i < Nd - 1; ++i)
      input.phase[i] = 0;
  }

  input.use_device_for_contractions = true;
  if (inputtop.count("use_device_for_contractions") == 1) {
    read(inputtop, "use_device_for_contractions",
         input.use_device_for_contractions);
  }

  input.use_superb_format = true;
  if (inputtop.count("use_superb_format") == 1) {
    read(inputtop, "use_superb_format", input.use_superb_format);
  }

  input.output_file_is_local = false;
  if (inputtop.count("output_file_is_local") == 1) {
    read(inputtop, "output_file_is_local", input.output_file_is_local);
  }
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const Params::Param_t::Contract_t &input) {
  push(xml, path);

  write(xml, "num_vecs", input.num_vecs);
  write(xml, "t_sources", input.t_sources);
  write(xml, "decay_dir", input.decay_dir);
  write(xml, "Nt_forward", input.Nt_forward);
  write(xml, "Nt_backward", input.Nt_backward);
  write(xml, "mass_label", input.mass_label);
  write(xml, "max_rhs", input.max_rhs);
  write(xml, "phase", input.phase);
  write(xml, "phases", input.phases);
  write(xml, "use_superb_format", input.use_superb_format);
  write(xml, "output_file_is_local", input.output_file_is_local);
  write(xml, "use_device_for_contractions", input.use_device_for_contractions);

  pop(xml);
}

//! Propagator input
void read(XMLReader &xml, const std::string &path, Params::Param_t &input) {
  XMLReader inputtop(xml, path);
  int i = 4;

  read(inputtop, "Propagator", input.prop);
  read(inputtop, "Contractions", input.contract);
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const Params::Param_t &input) {
  push(xml, path);

  write(xml, "Propagator", input.prop);
  write(xml, "Contractions", input.contract);

  pop(xml);
}

//! Propagator input
void read(XMLReader &xml, const std::string &path, Params &input) {
  Params tmp(xml, path);
  input = tmp;
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path, const Params &input) {
  push(xml, path);

  write(xml, "Param", input.param);
  write(xml, "NamedObject", input.named_obj);

  pop(xml);
}

//----------------------------------------------------------------------------
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlineMeas(Params(xml_in, path));
}

//! Local registration flag
bool registered = false;

const std::string name = "PROP_AND_MATELEM_DISTILLATION_SUPERB";

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= WilsonTypeFermActsEnv::registerAll();
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);
    registered = true;
  }
  return success;
}

// Geometry helpers are shared with the child-local factorized measurements; see
// util/ferm/factorized_distillation.h
using FactorizedDistillation::FactorizedPropGeometry;
using FactorizedDistillation::writeFactorizedPropGeometry;

bool factorizedPropRequested(const Params::NamedObject_t &named_obj) {
  return !named_obj.frozen_boundary_sidecar_file.empty() ||
         named_obj.frozen_boundary_child_id >= 0 ||
         !named_obj.factorized_intermediate_file.empty() ||
         !named_obj.factorized_cross_file.empty();
}

FactorizedPropGeometry
loadFactorizedPropGeometry(const Params::NamedObject_t &named_obj,
                           int decay_dir) {
  if (!factorizedPropRequested(named_obj))
    return FactorizedPropGeometry();
  return FactorizedDistillation::loadFactorizedPropGeometry(
      named_obj.frozen_boundary_sidecar_file,
      named_obj.frozen_boundary_child_id, decay_dir);
}

//----------------------------------------------------------------------------
// Param stuff
Params::Params() { frequency = 0; }

Params::Params(XMLReader &xml_in, const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    // Parameters for source construction
    read(paramtop, "Param", param);

    // Read in the output propagator/source configuration info
    read(paramtop, "NamedObject", named_obj);

    // Possible alternate XML file pattern
    if (paramtop.count("xml_file") != 0) {
      read(paramtop, "xml_file", xml_file);
    }
  } catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Function call
void InlineMeas::operator()(unsigned long update_no, XMLWriter &xml_out) {
  // If xml file not empty, then use alternate
  if (params.xml_file != "") {
    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

    push(xml_out, "PropDistillation");
    write(xml_out, "update_no", update_no);
    write(xml_out, "xml_file", xml_file);
    pop(xml_out);

    XMLFileWriter xml(xml_file);
    func(update_no, xml);
  } else {
    func(update_no, xml_out);
  }
}

// Lattice-field helpers are shared with the child-local factorized measurements
using FactorizedDistillation::restrictToTimeslices;
using FactorizedDistillation::returnNLatticeFermions;
using FactorizedDistillation::toLatticeFermions;
using FactorizedDistillation::toSBTensor;
using FactorizedDistillation::zeroTemporalLinksOnSlice;

//! Per-spin distillation projector on the given time slices
/*!
 * out = sum_{t in slices} sum_m phi_m(t) phi_m(t)^dag psi, the identity in spin.
 * out is zero on every other slice. `colorvec_proj` has order cxyzXnt and holds
 * the projector vectors phi_m for all Lt slices.
 */
static SB::MultipleLatticeFermions
projectOnSlices(const SB::MultipleLatticeFermions &psi,
                const std::vector<int> &slices,
                SB::Tensor<Nd + 3, SB::Complex> colorvec_proj) {
  auto projected = returnNLatticeFermions(psi.size());
  for (int i = 0; i < psi.size(); i++) {
    for (int t : slices) {
      auto psi_sb = toSBTensor(
          SB::MultipleLatticeFermions(psi.begin() + i, psi.begin() + i + 1), t,
          1);
      auto proj = colorvec_proj.kvslice_from_size({{'t', t}}, {{'t', 1}});
      auto inner = SB::contract<3>(proj.conj().rename_dims({{'n', 'm'}}),
                                   psi_sb, "cxyzXt");
      auto outer =
          SB::contract<8>(proj.rename_dims({{'n', 'm'}}), inner, "m");
      outer.copyTo(SB::asTensorView(*projected[i])
                       .kvslice_from_size({{'t', t}}, {{'t', 1}}));
    }
  }
  return projected;
}

// Real work done here
void InlineMeas::func(unsigned long update_no, XMLWriter &xml_out) {
  START_CODE();

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  // Parse the phase
  std::vector<SB::Coor<3>> phasings; ///< list of all phasings
  std::map<int, std::set<int>>
      phasing_pairs; ///< phasing src -> list of phasing sink
  if (params.param.contract.phase.size() > 0) {
    phasings.push_back(SB::toCoor(params.param.contract.phase));
    phasing_pairs[0] = std::set<int>{0};
  } else if (params.param.contract.phases.size() > 0) {
    if (!params.param.contract.use_superb_format) {
      throw std::runtime_error(
          "only the new format has support for `phases'; please set "
          "'use_superb_format' to true");
    }
    for (const auto phasing_pair : params.param.contract.phases) {
      auto src = SB::toCoor(phasing_pair.src);
      auto snk_neg = SB::toCoor(phasing_pair.snk);
      auto snk = SB::Coor<3>{-snk_neg[0], -snk_neg[1], -snk_neg[2]};
      int src_idx =
          std::find(phasings.begin(), phasings.end(), src) - phasings.begin();
      if (src_idx == phasings.size())
        phasings.push_back(src);
      int snk_idx =
          std::find(phasings.begin(), phasings.end(), snk) - phasings.begin();
      if (snk_idx == phasings.size())
        phasings.push_back(snk);
      if (phasing_pairs.count(src_idx) == 0) {
        phasing_pairs[src_idx] = std::set<int>{snk_idx};
      } else {
        phasing_pairs[src_idx].insert(snk_idx);
      }
    }
  } else {
    phasings.push_back(SB::Coor<Nd - 1>{{}});
    phasing_pairs[0] = std::set<int>{0};
  }

  // Test and grab a reference to the gauge field
  multi1d<LatticeColorMatrix> u;
  XMLBufferWriter gauge_xml;
  try {
    u = TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(
        params.named_obj.gauge_id);
    TheNamedObjMap::Instance()
        .get(params.named_obj.gauge_id)
        .getRecordXML(gauge_xml);
  } catch (std::bad_cast) {
    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  } catch (const std::string &e) {
    QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
    QDP_abort(1);
  }

  push(xml_out, "PropDistillation");
  write(xml_out, "update_no", update_no);

  QDPIO::cout << name << ": propagator calculation" << std::endl;

  proginfo(xml_out); // Print out basic program info

  // Write out the input
  write(xml_out, "Input", params);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);

  // Will use TimeSliceSet-s a lot
  const int decay_dir = params.param.contract.decay_dir;
  const int Lt = Layout::lattSize()[decay_dir];

  // A sanity check
  if (decay_dir != Nd - 1) {
    QDPIO::cerr << name << ": TimeSliceIO only supports decay_dir= " << Nd - 1
                << "\n";
    QDP_abort(1);
  }

  const FactorizedPropGeometry factorized_geometry =
      loadFactorizedPropGeometry(params.named_obj, decay_dir);
  writeFactorizedPropGeometry(xml_out, factorized_geometry);

  if (params.named_obj.overlapping_factorization &&
      !factorized_geometry.enabled) {
    QDPIO::cerr << name
                << ": overlapping_factorization requires the factorized path "
                   "(frozen_boundary_sidecar_file and frozen_boundary_child_id)"
                << std::endl;
    QDP_abort(1);
  }

  //
  // Read in the source along with relevant information.
  //

  SB::ColorvecsStorage colorvecsSto =
      SB::openColorvecStorage(params.named_obj.colorvec_files);

  //
  // DB storage
  //
  std::vector<BinaryStoreDB<SerialDBKey<KeyPropElementalOperator_t>,
                            SerialDBData<ValPropElementalOperator_t>>>
      qdp_db{};
  SB::StorageTensor<8, SB::ComplexD> st;

  // Open the file, and write the meta-data and the binary for this operator
  if (!params.param.contract.use_superb_format) {
    qdp_db.resize(1);
    if (!qdp_db[0].fileExists(params.named_obj.prop_op_file)) {
      XMLBufferWriter file_xml;
      push(file_xml, "DBMetaData");
      write(file_xml, "id", std::string("propElemOp"));
      write(file_xml, "lattSize", QDP::Layout::lattSize());
      write(file_xml, "decay_dir", params.param.contract.decay_dir);
      proginfo(file_xml); // Print out basic program info
      write(file_xml, "Params", params.param);
      write(file_xml, "Config_info", gauge_xml);
      pop(file_xml);

      std::string file_str(file_xml.str());
      qdp_db[0].setMaxUserInfoLen(file_str.size());

      qdp_db[0].open(params.named_obj.prop_op_file, O_RDWR | O_CREAT, 0664);

      qdp_db[0].insertUserdata(file_str);
    } else {
      qdp_db[0].open(params.named_obj.prop_op_file, O_RDWR, 0664);
    }
  } else {
    // Read order; letter meaning:
    //   n/N: source/sink eigenvector index
    //   s/q: source/sink spin index
    //   p/P: time slice source/sink
    //   h/H: phasing source/sink
    const char *order = "sqnNpPhH";
    XMLBufferWriter metadata_xml;
    push(metadata_xml, "DBMetaData");
    write(metadata_xml, "id", std::string("propElemOp"));
    write(metadata_xml, "lattSize", QDP::Layout::lattSize());
    write(metadata_xml, "decay_dir", params.param.contract.decay_dir);
    proginfo(metadata_xml); // Print out basic program info
    write(metadata_xml, "Config_info", gauge_xml);
    write(metadata_xml, "Params", params.param);
    write(metadata_xml, "mass_label", params.param.contract.mass_label);
    write(metadata_xml, "tensorOrder", order);
    std::vector<multi1d<int>> phasings_xml;
    for (const auto &it : phasings)
      phasings_xml.push_back(SB::tomulti1d(it));
    write(metadata_xml, "phasings", phasings_xml);
    pop(metadata_xml);

    // NOTE: metadata_xml only has a valid value on Master node; so do a
    // broadcast
    std::string metadata = SB::broadcast(metadata_xml.str());

    st = SB::StorageTensor<8, SB::ComplexD>(
        params.named_obj.prop_op_file, metadata, order,
        SB::kvcoors<8>(order,
                       {
                           {'s', Ns},
                           {'q', Ns},
                           {'n', params.param.contract.num_vecs},
                           {'N', params.param.contract.num_vecs},
                           {'p', Lt},
                           {'P', Lt},
                           {'h', (int)phasings.size()},
                           {'H', (int)phasings.size()},
                       }),
        SB::Sparse, SB::checksum_type::BlockChecksum,
        params.param.contract.output_file_is_local ? SB::LocalFSFile
                                                   : SB::SharedFSFile);
  }

  QDPIO::cout << "Finished opening peram file" << std::endl;

  //
  // Try the factories
  //
  try {
    StopWatch swatch;
    swatch.reset();
    QDPIO::cout << "Try the various factories" << std::endl;

    //
    // Initialize fermion action and create the solver
    //
    SB::ChimeraSolver PP{params.param.prop.fermact, params.param.prop.invParam,
                         u};
    swatch.start();

    //
    // Loop over the source color and spin, creating the source
    // and calling the relevant propagator routines.
    //
    const int num_vecs = params.param.contract.num_vecs;
    const multi1d<int> &t_sources = params.param.contract.t_sources;
    const int max_rhs = params.param.contract.max_rhs;

    // Set place for doing the contractions
    SB::DeviceHost dev = params.param.contract.use_device_for_contractions
                             ? SB::OnDefaultDevice
                             : SB::OnHost;

    if (factorized_geometry.enabled) {

      FactorizedPropGeometry geometry =
          loadFactorizedPropGeometry(params.named_obj, decay_dir);

      // Overlapping (Giusti) factorization, see
      // specs/hier/factorized_prop_overlapping.md
      const bool overlapping = params.named_obj.overlapping_factorization;
      // Overlapping mode allows num_intermediate_projectors <= 0: no projector
      const bool use_projector =
          !overlapping || params.named_obj.num_intermediate_projectors > 0;

      // Overlapping geometry: projector slices w2 (source side of the hop) and
      // hop target slices w1 (sink side of the hop)
      std::vector<int> ov_proj_slices, ov_hop_slices;
      multi1d<int> ov_proj_slices_m(2), ov_hop_slices_m(2);

      // Non-overlapping: PP_sep on u_sep. Overlapping: PP_src, PP_snk.
      std::unique_ptr<SB::ChimeraSolver> PP_sep, PP_src, PP_snk;

      if (!overlapping) {
        // Enforce open boundary conditions through gauge link zeroing
        // Zero timeslices where the active zone reaches the boundary
        multi1d<LatticeColorMatrix> u_sep = u;

        // Assuming 2 regions:
        // Moving measurement D00 one slice in from frozen region
        // Eliminates clover term crossing into active regions
        // Note this assumes a boundary width of at least 3; enforced with error
        // handling
        for (int i = 0; i < geometry.plan.child0_frozen_local_intervals.size();
             i++) {
          // TODO:Handling for end of lattice below?
          u_sep = zeroTemporalLinksOnSlice(
              u_sep,
              SB::normalize_coor(
                  geometry.plan.child0_frozen_local_intervals[i].t_start, Lt));
          u_sep = zeroTemporalLinksOnSlice(
              u_sep, geometry.plan.child0_frozen_local_intervals[i].t_end - 1);
          // std::cout << "Freezing t links on slice: " <<
          // SB::normalize_coor(geometry.plan.child0_frozen_local_intervals[i].t_start
          // - 1, Lt)  << std::endl;
        }

        // SB::asTensorView(u_sep[3]).print("u_sep");
        // Printing without even/odd ordering
        SB::detail::toNaturalOrdering(SB::asTensorView(u_sep[3]).toComplex())
            .print("u_sep");
        // Compute prop based on open boundary conditions
        PP_sep.reset(new SB::ChimeraSolver(params.param.prop.fermact,
                                           params.param.prop.invParam, u_sep));
      } else {
        const multi1d<GaugeSubdomainSplitInterval> &fi =
            geometry.plan.child0_frozen_local_intervals;
        if (fi.size() != 2)
          throw std::runtime_error(
              "overlapping_factorization: expected 2 frozen intervals");
        const int fw = fi[0].t_end - fi[0].t_start + 1;
        const int c1 = fi[1].t_start;
        if (fi[1].t_end - fi[1].t_start + 1 != fw)
          throw std::runtime_error("overlapping_factorization: frozen "
                                   "intervals have different widths");
        if (fi[0].t_start != 0)
          throw std::runtime_error("overlapping_factorization: first frozen "
                                   "interval must start at t = 0 (cut0 = 0)");
        if (fw < 3)
          throw std::runtime_error(
              "overlapping_factorization: frozen width must be >= 3");
        if (c1 <= fw || c1 + fw > Lt - 1)
          throw std::runtime_error("overlapping_factorization: empty source "
                                   "or sink active region");
        for (int tt = 0; tt < t_sources.size(); ++tt) {
          if (t_sources[tt] < fw || t_sources[tt] > c1 - 1)
            throw std::runtime_error(
                "overlapping_factorization: t_source " +
                std::to_string(t_sources[tt]) +
                " is outside the source active region [" + std::to_string(fw) +
                ", " + std::to_string(c1 - 1) + "]");
        }

        // Projector on the first internal frozen slice from the source side,
        // hop one slice toward the sink
        ov_proj_slices = {fw - 2, c1 + 1};
        ov_hop_slices = {fw - 3, c1 + 2};
        for (int i = 0; i < 2; ++i) {
          ov_proj_slices_m[i] = ov_proj_slices[i];
          ov_hop_slices_m[i] = ov_hop_slices[i];
        }

        // Omega_src = [0, c1+fw-1]: cut F1 -> sink active and sink active -> F0
        multi1d<int> src_cut_links(2);
        src_cut_links[0] = c1 + fw - 1;
        src_cut_links[1] = Lt - 1;
        // Omega_snk = [c1+2, Lt-1] u [0, fw-3]: cut the hop links
        multi1d<int> snk_cut_links(2);
        snk_cut_links[0] = c1 + 1;
        snk_cut_links[1] = fw - 3;

        multi1d<LatticeColorMatrix> u_src = u;
        multi1d<LatticeColorMatrix> u_snk = u;
        for (int i = 0; i < 2; ++i) {
          u_src = zeroTemporalLinksOnSlice(u_src, src_cut_links[i], decay_dir);
          u_snk = zeroTemporalLinksOnSlice(u_snk, snk_cut_links[i], decay_dir);
        }

        PP_src.reset(new SB::ChimeraSolver(params.param.prop.fermact,
                                           params.param.prop.invParam, u_src));
        PP_snk.reset(new SB::ChimeraSolver(params.param.prop.fermact,
                                           params.param.prop.invParam, u_snk));

        QDPIO::cout << name << ": overlapping factorization, fw = " << fw
                    << ", projector slices " << ov_proj_slices[0] << " "
                    << ov_proj_slices[1] << ", hop target slices "
                    << ov_hop_slices[0] << " " << ov_hop_slices[1]
                    << ", projector "
                    << (use_projector ? "on" : "off (Pi = 1)") << std::endl;

        push(xml_out, "OverlappingFactorization");
        write(xml_out, "frozen_width", fw);
        write(xml_out, "cut1", c1);
        write(xml_out, "projector_slices", ov_proj_slices_m);
        write(xml_out, "hop_target_slices", ov_hop_slices_m);
        write(xml_out, "src_zeroed_link_slices", src_cut_links);
        write(xml_out, "snk_zeroed_link_slices", snk_cut_links);
        write(xml_out, "use_projector", use_projector);
        pop(xml_out);
      }

      // Build Dirac operator for boundary crossing
      std::istringstream xml_s(params.param.prop.fermact.xml);
      XMLReader fermacttop(xml_s);
      auto S_f(TheFermionActionFactory::Instance().createObject(
          params.param.prop.fermact.id, fermacttop,
          params.param.prop.fermact.path));
      auto state(S_f->createState(u));
      LinearOperator<LatticeFermion> *fLinOp = S_f->genLinOp(state);

      // Loop over time sources
      for (int tt = 0; tt < t_sources.size(); ++tt) {
        int t_source = t_sources[tt]; // This is the actual time-slice.
        QDPIO::cout << "t_source = " << t_source << std::endl;

        // ends of source subdomain
        // Measuring one timeslice in from frozen boundary
        // TODO: Needs case handling
        multi1d<int> inner_source_boundaries(4);
        inner_source_boundaries[0] = SB::normalize_coor(
            geometry.plan.child0_frozen_local_intervals[0].t_end - 1, Lt);
        inner_source_boundaries[1] = SB::normalize_coor(
            geometry.plan.child0_frozen_local_intervals[0].t_start + 1, Lt);
        inner_source_boundaries[2] = SB::normalize_coor(
            geometry.plan.child0_frozen_local_intervals[1].t_end - 1, Lt);
        inner_source_boundaries[3] = SB::normalize_coor(
            geometry.plan.child0_frozen_local_intervals[1].t_start + 1, Lt);

        // Simple case
        int num_tslices = Lt;

        // For now for simplicity
        int first_tslice = 0;

        // Create full fermionic vector for t_source

        // Get colorvecs
        SB::Tensor<Nd + 3, SB::Complex> colorvec =
            SB::getColorvecs<SB::Complex>(colorvecsSto, u, decay_dir, 0, Lt,
                                          num_vecs, "cxyzXnt", SB::Coor<3>{{}},
                                          dev);

        SB::Tensor<Nd + 3, SB::Complex> colorvec_proj =
            use_projector
                ? SB::getColorvecs<SB::Complex>(
                      colorvecsSto, u, decay_dir, 0, Lt,
                      params.named_obj.num_intermediate_projectors, "cxyzXnt",
                      SB::toCoor(params.named_obj.intermediate_projector_phase),
                      dev)
                : SB::Tensor<Nd + 3, SB::Complex>();

        // Loop over phasings/spin indexes
        for (const auto &it : phasing_pairs) {
          int phasing_src_idx = it.first;

          // Get all eigenvectors for `t_source`
          auto source_colorvec = SB::phaseColorvecs(
              colorvec.kvslice_from_size({{'t', t_source}}, {{'t', 1}}),
              t_source, phasings[phasing_src_idx]);
          for (int spin_source = 0; spin_source < Ns; ++spin_source) {
            SB::MultipleLatticeFermions sink_solns;
            if (overlapping) {
              // tau = - V^dag(t') D_snk^-1 D_{w1<-w2} Pi_{w2} D_src^-1 V(t)
              // (arXiv:1601.04587, Eq. 5.3)

              // psi_src = D_src^-1 V(t_source), all time slices
              auto psi_src_sb =
                  SB::doInversion(*PP_src, source_colorvec, t_source,
                                  first_tslice, Lt, {spin_source}, max_rhs,
                                  "cxyzXnSst")
                      .rename_dims({{'S', 's'}, {'s', 'S'}});
              auto psi_src = toLatticeFermions(psi_src_sb, 0);

              // Pi_{w2} psi_src; without projector only restrict to w2
              auto psi_w2 =
                  use_projector
                      ? projectOnSlices(psi_src, ov_proj_slices, colorvec_proj)
                      : restrictToTimeslices(psi_src, ov_proj_slices_m,
                                             decay_dir);
              psi_src.clear();

              // Hop D_{w1<-w2}: full D (uncut u), restricted to the w1 slices
              auto hop_full = returnNLatticeFermions(psi_w2.size());
              for (int i = 0; i < psi_w2.size(); i++) {
                (*fLinOp)(*hop_full[i], *psi_w2[i], PLUS);
              }
              psi_w2.clear();
              auto hop_w1 =
                  restrictToTimeslices(hop_full, ov_hop_slices_m, decay_dir);
              hop_full.clear();

              // Single block-inverse minus sign
              for (int i = 0; i < hop_w1.size(); i++) {
                *hop_w1[i] = -(*hop_w1[i]);
              }

              // psi_snk = D_snk^-1 hop_w1
              sink_solns = returnNLatticeFermions(hop_w1.size());
              SB::doInversion(*PP_snk, sink_solns,
                              Chroma::SB::ConstMultipleLatticeFermions(
                                  hop_w1.begin(), hop_w1.end()),
                              max_rhs);
              hop_w1.clear();
            } else {
              // Contract source timeslice with D_{11}^-1
              // return all timeslices for now
              auto contract1_sb =
                  SB::doInversion(*PP_sep, source_colorvec, t_source, first_tslice, Lt,
                                  {spin_source}, max_rhs, "cxyzXnSst")
                      .rename_dims({{'S', 's'}, {'s', 'S'}});

              auto contract1 = toLatticeFermions(contract1_sb, 0);
              // SB::detail::toNaturalOrdering(SB::asTensorView(*contract1[0]).toComplex()).print("contract1");

              // Multiply result with boundary Dirac matrix D_{01}
              /// y = Dslash * x (*fLinOp)(x, y, PLUS);
              auto y_boundary1 = returnNLatticeFermions(contract1.size());

              for (int i = 0; i < contract1.size(); i++) {
                (*fLinOp)(*y_boundary1[i], *contract1[i], PLUS);
              }

              contract1.clear();

              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*y_boundary1[0]).toComplex())
                  .print("y_boundary1");

              // restrict solution to frozen regions
              auto y_boundary1_rs = restrictToTimeslices(
                  y_boundary1, inner_source_boundaries, decay_dir);

              y_boundary1.clear();

              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*y_boundary1_rs[0]).toComplex())
                  .print("y_boundary1_rs");

              // contract2 = D_{00}^-1 y_boundary1
              auto contract2 = returnNLatticeFermions(y_boundary1_rs.size());

              SB::doInversion(*PP_sep, contract2,
                              Chroma::SB::ConstMultipleLatticeFermions(
                                  y_boundary1_rs.begin(), y_boundary1_rs.end()),
                              max_rhs);

              y_boundary1_rs.clear();

              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*contract2[0]).toComplex())
                  .print("contract2");

              // y_boundary2 = D_{10} contract2
              auto y_boundary2 = returnNLatticeFermions(contract2.size());

              for (int i = 0; i < contract2.size(); i++) {
                (*fLinOp)(*y_boundary2[i], *contract2[i], PLUS);
              }

              contract2.clear();

              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*y_boundary2[0]).toComplex())
                  .print("y_boundary2");

              // restrict solution again
              // Include frozen timeslices internal to sink region
              // Needs case handling
              multi1d<int> sink_active_region(
                  Lt + 1 - geometry.plan.child0_frozen_local_intervals[1].t_end);
              for (int i = 0; i < sink_active_region.size() - 1; i++) {
                sink_active_region[i] =
                    geometry.plan.child0_frozen_local_intervals[1].t_end + i;
              }
              sink_active_region[sink_active_region.size() - 1] = 0;
              auto y_boundary2_rs = restrictToTimeslices(
                  y_boundary2, sink_active_region, decay_dir);
              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*y_boundary2_rs[0]).toComplex())
                  .print("y_boundary2_rs");

              y_boundary2.clear();

              // Contract with set of projection vectors
              //  L(t) = y_boundary2_rs*phi_x,m
              // on the first and last slice of the sink block
              auto projected_fields = projectOnSlices(
                  y_boundary2_rs,
                  {sink_active_region[0],
                   sink_active_region[sink_active_region.size() - 1]},
                  colorvec_proj);

              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*projected_fields[0]).toComplex())
                  .print("projected_fields");

              // contract3 = D_{11}^-1 y_boundary2
              auto contract3 = returnNLatticeFermions(y_boundary2_rs.size());
              // 2nd copy of contraction for debugging between 'exact' and reduced
              // auto contract3_db = returnNLatticeFermions(y_boundary2.size());
              // Below is computed w/o intermediary projector
              // SB::doInversion(PP, contract3_db,
              // Chroma::SB::ConstMultipleLatticeFermions
              // (y_boundary2_rs.begin(),y_boundary2_rs.end()), max_rhs);

              SB::doInversion(
                  *PP_sep, contract3,
                  Chroma::SB::ConstMultipleLatticeFermions(
                      projected_fields.begin(), projected_fields.end()),
                  max_rhs);

              SB::detail::toNaturalOrdering(
                  SB::asTensorView(*contract3[0]).toComplex())
                  .print("contract3");

              projected_fields.clear();

              // SB::detail::toNaturalOrdering(SB::asTensorView(*contract3_db[0]).toComplex()).print("contract3_db");

              sink_solns = std::move(contract3);
            }

            auto quark_solns =
                toSBTensor(sink_solns, first_tslice, num_tslices);
            // t_sink * contract3

            // Rest of the code like usual
            StopWatch snarss1;
            snarss1.reset();
            snarss1.start();

            for (int phasing_snk_idx : it.second) {
              // Phase the sink
              auto colorvec_snk = SB::phaseColorvecs(colorvec, first_tslice,
                                                     phasings[phasing_snk_idx]);

              // Contract the distillation elements
              // NOTE: N: is colorvec in sink, and n is colorvec in source
              SB::Tensor<5, SB::Complex> elems(
                  "NnSst", {num_vecs, num_vecs, Ns, 1, num_tslices}, SB::OnHost,
                  !params.param.contract.use_superb_format ? SB::OnMaster
                                                           : SB::OnEveryone);
              elems.contract(colorvec_snk, {{'n', 'N'}}, SB::Conjugate,
                             quark_solns.append_dimension('S'),
                             {{'s', 'S'}, {'S', 's'}}, SB::NotConjugate);

              snarss1.stop();
              QDPIO::cout << "Time to contract for one spin source : "
                          << snarss1.getTimeInSeconds() << " secs" << std::endl;

              snarss1.reset();
              snarss1.start();

              if (!params.param.contract.use_superb_format) {
                ValPropElementalOperator_t val;
                val.mat.resize(num_vecs, num_vecs);
                val.mat = zero;
                auto local_elems = elems.getLocal();
                for (int i_tslice = 0; i_tslice < num_tslices; ++i_tslice) {
                  for (int spin_sink = 0; spin_sink < Ns; ++spin_sink) {
                    KeyPropElementalOperator_t key;
                    key.t_source = t_source;
                    key.t_slice =
                        SB::normalize_coor(i_tslice + first_tslice, Lt);
                    key.spin_src = spin_source;
                    key.spin_snk = spin_sink;
                    key.mass_label = params.param.contract.mass_label;
                    if (local_elems) {
                      for (int colorvec_sink = 0; colorvec_sink < num_vecs;
                           ++colorvec_sink) {
                        for (int colorvec_source = 0;
                             colorvec_source < num_vecs; ++colorvec_source) {
                          std::complex<REAL> e =
                              local_elems.get({colorvec_sink, colorvec_source,
                                               spin_sink, 0, i_tslice});
                          val.mat(colorvec_sink, colorvec_source)
                              .elem()
                              .elem()
                              .elem() = RComplex<REAL64>(e.real(), e.imag());
                        }
                      }
                    }
                    qdp_db[0].insert(key, val);
                  }
                }
              } else {
                st.kvslice_from_size({{'s', spin_source},
                                      {'p', t_source},
                                      {'P', first_tslice},
                                      {'h', phasing_src_idx},
                                      {'H', phasing_snk_idx}},
                                     {{'s', 1},
                                      {'p', 1},
                                      {'P', num_tslices},
                                      {'h', 1},
                                      {'H', 1}})
                    .copyFrom(elems.rename_dims({{'S', 'q'}, {'t', 'P'}}));
              }
            }

            snarss1.stop();
            QDPIO::cout << "Time to store the props : "
                        << snarss1.getTimeInSeconds() << " secs" << std::endl;
          } // for spin
        } // for phasing
      } // for timesource
    } // if factorized
    else {
      // Loop over each time source
      for (int tt = 0; tt < t_sources.size(); ++tt) {
        int t_source = t_sources[tt]; // This is the actual time-slice.
        QDPIO::cout << "t_source = " << t_source << std::endl;

        // Compute the first tslice and the number of tslices involved in the
        // contraction
        int first_tslice = SB::normalize_coor(
            t_source - params.param.contract.Nt_backward, Lt);
        int num_tslices =
            std::min(params.param.contract.Nt_backward +
                         std::max(1, params.param.contract.Nt_forward),
                     Lt);

        // Get `num_vecs` colorvecs, and `num_tslices` tslices starting from
        // time-slice `first_tslice`
        SB::Tensor<Nd + 3, SB::Complex> colorvec =
            SB::getColorvecs<SB::Complex>(colorvecsSto, u, decay_dir,
                                          first_tslice, num_tslices, num_vecs,
                                          "cxyzXnt", SB::Coor<3>{{}}, dev);

        for (const auto &it : phasing_pairs) {
          int phasing_src_idx = it.first;

          // Get all eigenvectors for `t_source`
          auto source_colorvec = SB::phaseColorvecs(
              colorvec.kvslice_from_size(
                  {{'t', SB::normalize_coor(t_source - first_tslice, Lt)}},
                  {{'t', 1}}),
              t_source, phasings[phasing_src_idx]);

          for (int spin_source = 0; spin_source < Ns; ++spin_source) {
            // Invert the source for `spin_source` spin and retrieve
            // `num_tslices` tslices starting from tslice `first_tslice` NOTE: s
            // is spin source, and S is spin sink
            SB::Tensor<Nd + 5, SB::Complex> quark_solns = SB::doInversion(
                PP, source_colorvec, t_source, first_tslice, num_tslices,
                {spin_source}, max_rhs, "cxyzXnSst");

            StopWatch snarss1;
            snarss1.reset();
            snarss1.start();

            for (int phasing_snk_idx : it.second) {
              // Phase the sink
              auto colorvec_snk = SB::phaseColorvecs(colorvec, first_tslice,
                                                     phasings[phasing_snk_idx]);

              // Contract the distillation elements
              // NOTE: N: is colorvec in sink, and n is colorvec in source
              SB::Tensor<5, SB::Complex> elems(
                  "NnSst", {num_vecs, num_vecs, Ns, 1, num_tslices}, SB::OnHost,
                  !params.param.contract.use_superb_format ? SB::OnMaster
                                                           : SB::OnEveryone);
              elems.contract(colorvec_snk, {{'n', 'N'}}, SB::Conjugate,
                             quark_solns, {}, SB::NotConjugate);

              snarss1.stop();
              QDPIO::cout << "Time to contract for one spin source : "
                          << snarss1.getTimeInSeconds() << " secs" << std::endl;

              snarss1.reset();
              snarss1.start();

              if (!params.param.contract.use_superb_format) {
                ValPropElementalOperator_t val;
                val.mat.resize(num_vecs, num_vecs);
                val.mat = zero;
                auto local_elems = elems.getLocal();
                for (int i_tslice = 0; i_tslice < num_tslices; ++i_tslice) {
                  for (int spin_sink = 0; spin_sink < Ns; ++spin_sink) {
                    KeyPropElementalOperator_t key;
                    key.t_source = t_source;
                    key.t_slice =
                        SB::normalize_coor(i_tslice + first_tslice, Lt);
                    key.spin_src = spin_source;
                    key.spin_snk = spin_sink;
                    key.mass_label = params.param.contract.mass_label;
                    if (local_elems) {
                      for (int colorvec_sink = 0; colorvec_sink < num_vecs;
                           ++colorvec_sink) {
                        for (int colorvec_source = 0;
                             colorvec_source < num_vecs; ++colorvec_source) {
                          std::complex<REAL> e =
                              local_elems.get({colorvec_sink, colorvec_source,
                                               spin_sink, 0, i_tslice});
                          val.mat(colorvec_sink, colorvec_source)
                              .elem()
                              .elem()
                              .elem() = RComplex<REAL64>(e.real(), e.imag());
                        }
                      }
                    }
                    qdp_db[0].insert(key, val);
                  }
                }
              } else {
                st.kvslice_from_size({{'s', spin_source},
                                      {'p', t_source},
                                      {'P', first_tslice},
                                      {'h', phasing_src_idx},
                                      {'H', phasing_snk_idx}},
                                     {{'s', 1},
                                      {'p', 1},
                                      {'P', num_tslices},
                                      {'h', 1},
                                      {'H', 1}})
                    .copyFrom(elems.rename_dims({{'S', 'q'}, {'t', 'P'}}));
              }
            }

            snarss1.stop();
            QDPIO::cout << "Time to store the props : "
                        << snarss1.getTimeInSeconds() << " secs" << std::endl;
          } // for spin_source
        } // for it
      } // for tt
    } // if not factorized

    swatch.stop();
    QDPIO::cout << "Propagators computed: time= " << swatch.getTimeInSeconds()
                << " secs" << std::endl;
  } catch (const std::exception &e) {
    QDP_error_exit("%s: caught exception: %s\n", name.c_str(), e.what());
  }

  // Close colorvecs storage
  SB::closeColorvecStorage(colorvecsSto);

  pop(xml_out); // prop_dist

  snoop.stop();
  QDPIO::cout << name << ": total time = " << snoop.getTimeInSeconds()
              << " secs" << std::endl;

  QDPIO::cout << name << ": ran successfully" << std::endl;

  END_CODE();
}

} // namespace InlinePropAndMatElemDistillationSuperbEnv

} // namespace Chroma

#endif // BUILD_SB
