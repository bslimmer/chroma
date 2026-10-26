/*! \file
 * \brief Edge-SVD boundary projector basis for the factorized perambulator
 *
 * See specs/hier/factorized_prop_svd_projectors.md.
 */

#include "meas/inline/hadron/inline_edge_svd_projectors_superb_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "fermact.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "util/ferm/factorized_distillation.h"
#include "util/ferm/mgproton.h"
#include "util/ferm/superb_contractions.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#include "actions/ferm/fermstates/fermstates.h"
#include "handle.h"
#include "linearop.h"

#include <algorithm>
#include <map>
#include <sstream>
#include <vector>

#ifdef BUILD_SB

namespace Chroma
{
  namespace InlineEdgeSVDProjectorsSuperbEnv
  {
    //! Parameter input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "rank", input.rank);

      input.oversampling = 16;
      if (inputtop.count("oversampling") == 1)
	read(inputtop, "oversampling", input.oversampling);

      input.power_iterations = 2;
      if (inputtop.count("power_iterations") == 1)
	read(inputtop, "power_iterations", input.power_iterations);

      input.seed = 11;
      if (inputtop.count("seed") == 1)
	read(inputtop, "seed", input.seed);

      input.overlapping_factorization = false;
      if (inputtop.count("overlapping_factorization") == 1)
	read(inputtop, "overlapping_factorization", input.overlapping_factorization);

      input.max_rhs = 8;
      if (inputtop.count("max_rhs") == 1)
	read(inputtop, "max_rhs", input.max_rhs);

      input.decay_dir = Nd - 1;
      if (inputtop.count("decay_dir") == 1)
	read(inputtop, "decay_dir", input.decay_dir);

      read(inputtop, "Propagator", input.prop);

      if (input.rank <= 0)
	throw std::runtime_error("EDGE_SVD_PROJECTORS_SUPERB: rank should be > 0");
      if (input.oversampling < 0)
	throw std::runtime_error("EDGE_SVD_PROJECTORS_SUPERB: oversampling should be >= 0");
      if (input.power_iterations < 0)
	throw std::runtime_error("EDGE_SVD_PROJECTORS_SUPERB: power_iterations should be >= 0");
    }

    //! Parameter output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);
      write(xml, "rank", input.rank);
      write(xml, "oversampling", input.oversampling);
      write(xml, "power_iterations", input.power_iterations);
      write(xml, "seed", input.seed);
      write(xml, "overlapping_factorization", input.overlapping_factorization);
      write(xml, "max_rhs", input.max_rhs);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "Propagator", input.prop);
      pop(xml);
    }

    //! Named object input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "frozen_boundary_sidecar_file", input.frozen_boundary_sidecar_file);
      read(inputtop, "frozen_boundary_child_id", input.frozen_boundary_child_id);
      read(inputtop, "svd_basis_file", input.svd_basis_file);
    }

    //! Named object output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);
      write(xml, "gauge_id", input.gauge_id);
      write(xml, "frozen_boundary_sidecar_file", input.frozen_boundary_sidecar_file);
      write(xml, "frozen_boundary_child_id", input.frozen_boundary_child_id);
      write(xml, "svd_basis_file", input.svd_basis_file);
      pop(xml);
    }

    //! Params output
    void write(XMLWriter& xml, const std::string& path, const Params& input)
    {
      push(xml, path);
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);
      pop(xml);
    }

    //----------------------------------------------------------------------------
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path)
    {
      return new InlineMeas(Params(xml_in, path));
    }

    //! Local registration flag
    bool registered = false;

    const std::string name = "EDGE_SVD_PROJECTORS_SUPERB";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }

    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path)
    {
      try
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	read(paramtop, "Param", param);
	read(paramtop, "NamedObject", named_obj);

	if (paramtop.count("xml_file") != 0)
	  read(paramtop, "xml_file", xml_file);
      } catch (const std::string& e)
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

    //----------------------------------------------------------------------------
    // Function call
    void InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "EdgeSVDProjectors");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }

    using FactorizedDistillation::computeEdgeSVDBasis;
    using FactorizedDistillation::EdgeSVDBasis;
    using FactorizedDistillation::EdgeSVDBasisInfo;
    using FactorizedDistillation::EdgeSVDFermionMap;
    using FactorizedDistillation::EdgeSVDSlices;
    using FactorizedDistillation::edgeSVDSlices;
    using FactorizedDistillation::FactorizedPropGeometry;
    using FactorizedDistillation::nonOverlapCutField;
    using FactorizedDistillation::restrictToTimeslices;
    using FactorizedDistillation::returnNLatticeFermions;
    using FactorizedDistillation::writeEdgeSVDBasis;
    using FactorizedDistillation::writeFactorizedPropGeometry;
    using FactorizedDistillation::zeroTemporalLinksOnSlice;

    // Real work done here
    void InlineMeas::func(unsigned long update_no, XMLWriter& xml_out)
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(
	  params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      } catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      } catch (const std::string& e)
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "EdgeSVDProjectors");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": edge-SVD projector basis" << std::endl;

      proginfo(xml_out); // Print out basic program info
      write(xml_out, "Input", params);
      write(xml_out, "Config_info", gauge_xml);
      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      const int decay_dir = params.param.decay_dir;
      if (decay_dir != Nd - 1)
      {
	QDPIO::cerr << name << ": only decay_dir = " << Nd - 1 << " is supported" << std::endl;
	QDP_abort(1);
      }
      const int Lt = Layout::lattSize()[decay_dir];
      const bool overlapping = params.param.overlapping_factorization;
      const int R = params.param.rank;
      const int V3 = Layout::vol() / Lt;
      if (R > 6 * V3)
	throw std::runtime_error(name + ": rank should be <= 6 V3 (the spin half that the hop "
				 "passes) = " +
				 std::to_string(6 * V3));

      // Geometry
      const FactorizedPropGeometry geometry = FactorizedDistillation::loadFactorizedPropGeometry(
	params.named_obj.frozen_boundary_sidecar_file, params.named_obj.frozen_boundary_child_id,
	decay_dir);
      writeFactorizedPropGeometry(xml_out, geometry);
      const EdgeSVDSlices sl = edgeSVDSlices(geometry, overlapping, Lt);

      // The overlapping layout requires an unsmeared Dirac operator (as the prop measurement)
      if (overlapping)
      {
	std::istringstream fs_xml(params.param.prop.fermact.xml);
	XMLReader fs_top(fs_xml);
	XMLReader fs_act(fs_top, params.param.prop.fermact.path);
	std::string ferm_state = "SIMPLE_FERM_STATE";
	if (fs_act.count("FermState") != 0)
	  read(fs_act, "FermState/Name", ferm_state);
	if (ferm_state != "SIMPLE_FERM_STATE")
	  throw std::runtime_error(name +
				   ": overlapping_factorization requires SIMPLE_FERM_STATE, got " +
				   ferm_state);
      }

      // Dirac operator on the uncut field, for the hops
      std::istringstream xml_s(params.param.prop.fermact.xml);
      XMLReader fermacttop(xml_s);
      auto S_f(TheFermionActionFactory::Instance().createObject(
	params.param.prop.fermact.id, fermacttop, params.param.prop.fermact.path));
      auto state(S_f->createState(u));
      Handle<LinearOperator<LatticeFermion>> fLinOp(S_f->genLinOp(state));

      // Slab solver: D00 (non-overlap, the prop measurement's u_sep) or D_F (overlap)
      multi1d<LatticeColorMatrix> u_slab;
      if (!overlapping)
      {
	u_slab = nonOverlapCutField(u, geometry, Lt);
      }
      else
      {
	u_slab = u;
	for (int t : sl.slab_cut_links)
	  u_slab = zeroTemporalLinksOnSlice(u_slab, t, decay_dir);
      }
      SB::ChimeraSolver PP{params.param.prop.fermact, params.param.prop.invParam, u_slab};

      const int max_rhs = params.param.max_rhs;
      auto to_m = [](const std::vector<int>& v) {
	multi1d<int> r(v.size());
	for (int i = 0; i < v.size(); ++i)
	  r[i] = v[i];
	return r;
      };
      auto restrict_to = [&](const SB::MultipleLatticeFermions& x,
			     const std::vector<int>& slices) {
	return restrictToTimeslices(x, to_m(slices), decay_dir);
      };
      // D (PLUS) or D^dag (MINUS) on the uncut field
      auto apply_D = [&](const SB::MultipleLatticeFermions& x, PlusMinus isign) {
	auto y = returnNLatticeFermions(x.size());
	for (int i = 0; i < x.size(); ++i)
	  (*fLinOp)(*y[i], *x[i], isign);
	return y;
      };
      auto apply_g5 = [&](const SB::MultipleLatticeFermions& x) {
	auto y = returnNLatticeFermions(x.size());
	for (int i = 0; i < x.size(); ++i)
	  *y[i] = Gamma(Ns * Ns - 1) * (*x[i]);
	return y;
      };
      auto solve_slab = [&](const SB::MultipleLatticeFermions& x) {
	auto y = returnNLatticeFermions(x.size());
	SB::doInversion(PP, y, Chroma::SB::ConstMultipleLatticeFermions(x.begin(), x.end()),
			max_rhs);
	return y;
      };

      // Basis = top right singular vectors on p of M : p -> q
      EdgeSVDFermionMap M, Mdag;
      if (!overlapping)
      {
	// T = R_p D R_{e_p} D00^-1 R_{e_s} D R_q, Pi acts on its output: M = T^dag
	// M y = R_q D^dag R_{e_s} g5 D00^-1 g5 R_{e_p} D^dag y
	M = [&](const SB::MultipleLatticeFermions& y) {
	  auto v = restrict_to(apply_D(y, MINUS), sl.a);
	  v = apply_g5(solve_slab(apply_g5(v)));
	  v = restrict_to(v, sl.b);
	  return restrict_to(apply_D(v, MINUS), sl.q);
	};
	// M^dag x = R_p D R_{e_p} D00^-1 R_{e_s} D x
	Mdag = [&](const SB::MultipleLatticeFermions& x) {
	  auto v = restrict_to(apply_D(x, PLUS), sl.b);
	  v = restrict_to(solve_slab(v), sl.a);
	  return restrict_to(apply_D(v, PLUS), sl.p);
	};
      }
      else
      {
	// G_F = R_q D_F^-1 R_{w1} D, Pi acts on its input (p = w2): M = G_F
	M = [&](const SB::MultipleLatticeFermions& x) {
	  auto v = restrict_to(apply_D(x, PLUS), sl.a);
	  return restrict_to(solve_slab(v), sl.q);
	};
	// M^dag y = R_p D^dag R_{w1} g5 D_F^-1 g5 y
	Mdag = [&](const SB::MultipleLatticeFermions& y) {
	  auto v = apply_g5(solve_slab(apply_g5(y)));
	  v = restrict_to(v, sl.a);
	  return restrict_to(apply_D(v, MINUS), sl.p);
	};
      }

      // Randomized SVD
      StopWatch sv_watch;
      sv_watch.reset();
      sv_watch.start();
      std::map<int, std::vector<double>> sv_values;
      std::map<int, double> sv_ortho_err;
      EdgeSVDBasis basis = computeEdgeSVDBasis(
	M, Mdag, sl.p, sl.q, R, params.param.oversampling, params.param.power_iterations,
	params.param.seed, std::max(1, max_rhs), sv_values, sv_ortho_err);
      sv_watch.stop();

      // Write the file
      EdgeSVDBasisInfo info;
      info.latt_size = Layout::lattSize();
      info.decay_dir = decay_dir;
      info.overlapping = overlapping;
      info.frozen_width = sl.frozen_width;
      info.cut1 = sl.cut1;
      info.projector_slices = sl.p;
      info.far_slices = sl.q;
      info.rank = R;
      info.oversampling = params.param.oversampling;
      info.power_iterations = params.param.power_iterations;
      info.seed = params.param.seed;
      for (int p : sl.p)
      {
	info.singular_values.push_back(sv_values.at(p));
	info.orthonormality_error.push_back(sv_ortho_err.at(p));
      }

      StopWatch io_watch;
      io_watch.reset();
      io_watch.start();
      writeEdgeSVDBasis(params.named_obj.svd_basis_file, basis, info, gauge_xml,
			params.param.prop);
      io_watch.stop();

      // Summary
      QDPIO::cout << name << ": rank " << R << ", SVD time " << sv_watch.getTimeInSeconds()
		  << " secs, write time " << io_watch.getTimeInSeconds() << " secs, file "
		  << params.named_obj.svd_basis_file << std::endl;
      push(xml_out, "Basis");
      write(xml_out, "svd_basis_file", params.named_obj.svd_basis_file);
      write(xml_out, "rank", R);
      write(xml_out, "svd_time_secs", sv_watch.getTimeInSeconds());
      write(xml_out, "write_time_secs", io_watch.getTimeInSeconds());
      for (int i = 0; i < sl.p.size(); ++i)
      {
	const std::vector<double>& sv = sv_values.at(sl.p[i]);
	multi1d<double> sv_m(sv.size());
	for (int j = 0; j < sv.size(); ++j)
	  sv_m[j] = sv[j];
	QDPIO::cout << name << ":   projector slice " << sl.p[i] << " (far slice " << sl.q[i]
		    << "): sigma_1 = " << (sv.size() > 0 ? sv[0] : 0.0) << ", sigma_R/sigma_1 = "
		    << (sv.size() >= R && sv[0] > 0 ? sv[R - 1] / sv[0] : 0.0)
		    << ", |phi^dag phi - 1|_max = " << sv_ortho_err.at(sl.p[i]) << std::endl;
	push(xml_out, "elem");
	write(xml_out, "projector_slice", sl.p[i]);
	write(xml_out, "far_slice", sl.q[i]);
	write(xml_out, "singular_values", sv_m);
	write(xml_out, "orthonormality_error", sv_ortho_err.at(sl.p[i]));
	pop(xml_out);
      }
      pop(xml_out);

      pop(xml_out); // EdgeSVDProjectors

      snoop.stop();
      QDPIO::cout << name << ": total time = " << snoop.getTimeInSeconds() << " secs"
		  << std::endl;
      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    }
  } // namespace InlineEdgeSVDProjectorsSuperbEnv
} // namespace Chroma

#endif // BUILD_SB
