/*! \file
 * \brief Shared helpers for factorized (domain-decomposed) distillation
 * measurements
 */

#include "util/ferm/factorized_distillation.h"

#ifdef BUILD_SB
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/key_timeslice_colorvec.h"
#endif

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Chroma {
namespace FactorizedDistillation {
void appendIntervalSlices(const multi1d<GaugeSubdomainSplitInterval> &intervals,
                          multi1d<int> &slices) {
  int total = 0;
  for (int i = 0; i < intervals.size(); ++i)
    total += intervals[i].t_end - intervals[i].t_start + 1;

  slices.resize(total);
  int p = 0;
  for (int i = 0; i < intervals.size(); ++i) {
    for (int t = intervals[i].t_start; t <= intervals[i].t_end; ++t)
      slices[p++] = t;
  }
}

multi1d<GaugeSubdomainSplitInterval> buildActiveIntervals(
    const multi1d<GaugeSubdomainSplitInterval> &frozen_intervals,
    int child_extent) {
  std::vector<GaugeSubdomainSplitInterval> active;
  int next_start = 0;

  for (int i = 0; i < frozen_intervals.size(); ++i) {
    if (next_start < frozen_intervals[i].t_start) {
      GaugeSubdomainSplitInterval interval;
      interval.t_start = next_start;
      interval.t_end = frozen_intervals[i].t_start - 1;
      active.push_back(interval);
    }

    next_start = frozen_intervals[i].t_end + 1;
  }

  if (next_start < child_extent) {
    GaugeSubdomainSplitInterval interval;
    interval.t_start = next_start;
    interval.t_end = child_extent - 1;
    active.push_back(interval);
  }

  multi1d<GaugeSubdomainSplitInterval> result(active.size());
  for (int i = 0; i < active.size(); ++i)
    result[i] = active[i];

  return result;
}

FactorizedPropGeometry
loadFactorizedPropGeometry(const std::string &sidecar_file, int child_id,
                           int decay_dir) {
  FactorizedPropGeometry geometry;

  if (sidecar_file.empty())
    throw std::runtime_error("factorized measurement requested without "
                             "frozen_boundary_sidecar_file");

  if (child_id < 0 || child_id > 1)
    throw std::runtime_error(
        "frozen_boundary_child_id must be 0 or 1 for a factorized "
        "measurement");

  GaugeSubdomainSplitPlan plan;
  XMLReader sidecar_in(sidecar_file);
  read(sidecar_in, "/GaugeSubdomainSplitInfo", plan);

  if (plan.param.t_dir != decay_dir) {
    std::ostringstream os;
    os << "factorized measurement decay_dir=" << decay_dir
       << " does not match sidecar t_dir=" << plan.param.t_dir;
    throw std::runtime_error(os.str());
  }

  // Check that the frozen width is at least 3 to allow space for measurement
  if (plan.param.frozen_width < 3) {
    std::ostringstream os;
    os << "frozen interval size " << plan.param.frozen_width
       << " insufficient; frozen intervals must be at least 3 timeslices";
    throw std::runtime_error(os.str());
  }

  geometry.enabled = true;
  geometry.sidecar_file = sidecar_file;
  geometry.child_id = child_id;
  geometry.t_dir = plan.param.t_dir;
  geometry.plan = plan;

  if (geometry.child_id == 0) {
    geometry.child_nrow = plan.child0_nrow;
    geometry.local_to_global_t = plan.child0_local_to_global_t;
    geometry.frozen_local_intervals = plan.child0_frozen_local_intervals;
  } else {
    geometry.child_nrow = plan.child1_nrow;
    geometry.local_to_global_t = plan.child1_local_to_global_t;
    geometry.frozen_local_intervals = plan.child1_frozen_local_intervals;
  }

  geometry.active_local_intervals = buildActiveIntervals(
      geometry.frozen_local_intervals, geometry.child_nrow[geometry.t_dir]);
  appendIntervalSlices(geometry.frozen_local_intervals,
                       geometry.gamma_star_local_slices);
  appendIntervalSlices(geometry.active_local_intervals,
                       geometry.gamma_local_slices);

  return geometry;
}

void writeFactorizedPropGeometry(XMLWriter &xml,
                                 const FactorizedPropGeometry &geometry) {
  push(xml, "FactorizedPropGeometry");
  write(xml, "enabled", geometry.enabled);
  if (geometry.enabled) {
    write(xml, "sidecar_file", geometry.sidecar_file);
    write(xml, "child_id", geometry.child_id);
    write(xml, "t_dir", geometry.t_dir);
    write(xml, "child_nrow", geometry.child_nrow);
    write(xml, "local_to_global_t", geometry.local_to_global_t);
    write(xml, "frozen_local_intervals", geometry.frozen_local_intervals);
    write(xml, "active_local_intervals", geometry.active_local_intervals);
    write(xml, "gamma_star_local_slices", geometry.gamma_star_local_slices);
    write(xml, "gamma_local_slices", geometry.gamma_local_slices);
  }
  pop(xml);
}

void getCutLeftRight(const GaugeSubdomainSplitPlan &plan, int &cut_left,
                     int &cut_right) {
  cut_left = std::min(plan.param.cut0, plan.param.cut1);
  cut_right = std::max(plan.param.cut0, plan.param.cut1);
}

multi1d<LatticeColorMatrix>
zeroTemporalLinksOnSlice(const multi1d<LatticeColorMatrix> &u, int t_slice,
                         int t_dir) {
  multi1d<LatticeColorMatrix> u_out(u.size());
  for (int mu = 0; mu < u.size(); ++mu)
    u_out[mu] = u[mu];

  // Boolean mask: true on sites whose coordinate in t_dir equals t_slice
  LatticeInteger t_coord = Layout::latticeCoordinate(t_dir);
  LatticeBoolean on_slice = (t_coord == t_slice);

  LatticeColorMatrix zero_mat = zero;

  u_out[t_dir] = where(on_slice, zero_mat, u[t_dir]);

  return u_out;
}

multi1d<LatticeColorMatrix>
buildCutLinkField(const multi1d<LatticeColorMatrix> &u,
                  const FactorizedPropGeometry &geometry) {
  const int t_dir = geometry.t_dir;
  const int Lt = Layout::lattSize()[t_dir];
  multi1d<LatticeColorMatrix> u_sep = u;
  for (int i = 0; i < geometry.frozen_local_intervals.size(); ++i) {
    int before = (geometry.frozen_local_intervals[i].t_start - 1 + Lt) % Lt;
    int last = geometry.frozen_local_intervals[i].t_end % Lt;
    u_sep = zeroTemporalLinksOnSlice(u_sep, before, t_dir);
    u_sep = zeroTemporalLinksOnSlice(u_sep, last, t_dir);
  }
  return u_sep;
}

int findParentTemporalCrossing(const FactorizedPropGeometry &geometry) {
  const int parent_Lt = geometry.plan.parent_nrow[geometry.t_dir];
  const multi1d<int> &l2g = geometry.local_to_global_t;
  for (int t = 0; t + 1 < l2g.size(); ++t)
    if (l2g[t] == parent_Lt - 1 && l2g[t + 1] == 0)
      return t;
  return -1;
}

multi1d<LatticeColorMatrix>
applyParentTemporalPhase(const multi1d<LatticeColorMatrix> &u,
                         const FactorizedPropGeometry &geometry,
                         const Real &phase) {
  multi1d<LatticeColorMatrix> u_out = u;
  const int t_star = findParentTemporalCrossing(geometry);
  if (t_star < 0)
    return u_out;

  LatticeInteger t_coord = Layout::latticeCoordinate(geometry.t_dir);
  LatticeBoolean on_slice = (t_coord == t_star);
  LatticeColorMatrix phased = phase * u[geometry.t_dir];
  u_out[geometry.t_dir] = where(on_slice, phased, u[geometry.t_dir]);
  return u_out;
}

#ifdef BUILD_SB
SB::MultipleLatticeFermions returnNLatticeFermions(int n) {
  SB::MultipleLatticeFermions r(n);
  for (int col = 0; col < n; col++) {
    r[col].reset(new LatticeFermion);
    *r[col] = zero;
  }
  return r;
}

SB::MultipleLatticeFermions
restrictToTimeslices(const SB::MultipleLatticeFermions &psi,
                     const multi1d<int> &keep_tslices, int t_dir) {
  LatticeInteger t_coord = Layout::latticeCoordinate(t_dir);

  LatticeBoolean keep = false; // starts false everywhere
  for (int i = 0; i < keep_tslices.size(); ++i)
    keep |= (t_coord == keep_tslices[i]);

  LatticeFermion zero_f = zero;

  SB::MultipleLatticeFermions r = returnNLatticeFermions(psi.size());
  for (int i = 0; i < psi.size(); i++)
    *r[i] = where(keep, *psi[i], zero_f);
  return r;
}

SB::MultipleLatticeFermions
toLatticeFermions(const SB::Tensor<Nd + 5, SB::Complex> &chi, int t_source) {
  const auto n = chi.kvdim().at('n');
  SB::MultipleLatticeFermions r(n);
  for (int col = 0; col < n; col++) {
    r[col].reset(new LatticeFermion);
    *r[col] = zero;
    chi.kvslice_from_size({{'n', col}}, {{'n', 1}})
        .copyTo(SB::asTensorView(*r[col]).kvslice_from_size({{'t', t_source}}));
  }
  return r;
}

SB::Tensor<Nd + 4, SB::Complex>
toSBTensor(const SB::MultipleLatticeFermions &chi, int t_source,
           int num_tslices) {
  SB::Tensor<Nd + 4, SB::Complex> r(
      "cxyztXns",
      SB::latticeSize<Nd + 4>("cxyztXns",
                              {{'t', num_tslices}, {'n', (int)chi.size()}}),
      SB::OnDefaultDevice, SB::OnEveryone);
  for (int col = 0; col < chi.size(); col++) {
    SB::asTensorView(*chi[col])
        .kvslice_from_size({{'t', t_source}}, {{'t', num_tslices}})
        .copyTo(r.kvslice_from_size({{'n', col}}, {{'n', 1}}));
  }
  return r;
}

void applyLinOp(const LinearOperator<LatticeFermion> &linop,
                SB::MultipleLatticeFermions &y,
                const SB::MultipleLatticeFermions &x) {
  if (y.size() != x.size())
    throw std::runtime_error("applyLinOp: input and output sizes differ");
  for (int i = 0; i < x.size(); i++)
    linop(*y[i], *x[i], PLUS);
}

namespace {
//! Deterministic, ordering-dependent checksum of a natural-ordering host tensor
//! on master
template <typename T>
double naturalChecksum(const SB::Tensor<Nd, T> &tnat, int v) {
  double r = 0;
  auto loc = tnat.getLocal();
  if (loc) {
    const T *p = loc.data();
    const std::size_t n = loc.volume();
    for (std::size_t i = 0; i < n; ++i)
      r += (double)((i % 97) + 1) * (double)(v + 1) *
           ((double)std::real(p[i]) + 2.0 * (double)std::imag(p[i]));
  }
  return r;
}
} // namespace

SB::Tensor<Nd + 3, SB::Complex>
readColorvecsOnSliceFromParentFile(const SB::ColorvecsStorage &sto,
                                   const multi1d<LatticeColorMatrix> &u,
                                   int global_t, int local_t, int n_colorvecs,
                                   int parent_Lt, double *checksum) {
  using namespace SB;

  double sum = 0;

  if (sto.is_dummy)
    throw std::runtime_error(
        "projector colorvecs cannot be read from a dummy colorvec file");
  if (global_t < 0 || global_t >= parent_Lt)
    throw std::runtime_error(
        "projector colorvec slice outside the parent time extent");

  const std::string order = "cxyzXnt";
  Tensor<Nd + 3, ComplexD> r(
      order, latticeSize<Nd + 3>(order, {{'t', 1}, {'n', n_colorvecs}}),
      OnHost);
  r.set_zero();

  // Single time-slice colorvec in red-black ordering for the child-local slice
  Tensor<Nd + 1, ComplexD> trb("cxyzX", latticeSize<Nd + 1>("cxyzX"), OnHost,
                               OnMaster);
  std::vector<Index> perm = ns_getColorvecs::getPermFromNatToRB(local_t);

  bool fingerprint = false;
  GroupXML_t link_smear;

  if (sto.s3t) {
    std::istringstream is(sto.s3t.metadata);
    XMLReader xml_buf(is);
    if (xml_buf.count("/MODMetaData/fingerprint") == 1)
      read(xml_buf, "/MODMetaData/fingerprint", fingerprint);
    if (fingerprint)
      link_smear = readXMLGroup(xml_buf, "/MODMetaData/LinkSmearing",
                                "LinkSmearingType");

    // Check that the file comes from a lattice with the same spatial extents
    // (only the stored corner for fingerprint files) and the parent time extent
    auto fdim = sto.s3t.kvdim();
    for (int d = 0; d < Nd - 1; ++d) {
      const char label = "xyz"[d];
      if ((!fingerprint && fdim.at(label) != Layout::lattSize()[d]) ||
          (fingerprint && fdim.at(label) > Layout::lattSize()[d]))
        throw std::runtime_error(
            "projector colorvec file spatial extents do not match the "
            "current lattice");
    }
    if (fdim.at('t') != parent_Lt)
      throw std::runtime_error(
          "projector colorvec file time extent does not match the parent "
          "lattice time extent of the split sidecar");
    if (fdim.at('n') < n_colorvecs)
      throw std::runtime_error(
          "projector colorvec file has fewer colorvecs than requested");

    Tensor<Nd, ComplexD> tnat(
        "cxyz", latticeSize<Nd>("cxyz", {{'x', Layout::lattSize()[0]}}), OnHost,
        OnMaster);
    for (int v = 0; v < n_colorvecs; ++v) {
      tnat.set_zero();
      sto.s3t
          .kvslice_from_size({{'t', global_t}, {'n', v}}, {{'t', 1}, {'n', 1}})
          .copyTo(tnat);
      sum += naturalChecksum(tnat, v);
      ns_getColorvecs::toRB(perm, tnat, trb);
      trb.copyTo(r.kvslice_from_size({{'n', v}}, {{'n', 1}}));
    }
  } else if (sto.mod) {
    Tensor<Nd, ComplexF> tnat(
        "cxyz", latticeSize<Nd>("cxyz", {{'x', Layout::lattSize()[0]}}), OnHost,
        OnMaster);
    Tensor<Nd + 1, ComplexF> trbf("cxyzX", latticeSize<Nd + 1>("cxyzX"), OnHost,
                                  OnMaster);
    for (int v = 0; v < n_colorvecs; ++v) {
      KeyTimeSliceColorVec_t key(global_t, v);
      if (!sto.mod->exist(key))
        throw std::runtime_error(
            "no projector colorvec exists with key t_slice= " +
            std::to_string(global_t) + " colorvec= " + std::to_string(v));
      sto.mod->get(key, tnat);
      sum += naturalChecksum(tnat, v);
      ns_getColorvecs::toRB(perm, tnat, trbf);
      trbf.copyTo(r.kvslice_from_size({{'n', v}}, {{'n', 1}}));
    }
  } else {
    throw std::runtime_error("projector colorvec storage is not open");
  }

  // Check that no vector is null
  auto norms = norm<1>(r, "n");
  for (int v = 0; v < n_colorvecs; ++v)
    if (norms.get({v}) == 0)
      throw std::runtime_error("null projector colorvec at parent t_slice= " +
                               std::to_string(global_t) +
                               " colorvec= " + std::to_string(v));

  if (fingerprint) {
    // Recompute the eigenvectors on the child-local slice from the (frozen)
    // spatial links and phase them so that they match the stored fingerprint,
    // as SB::getColorvecs does
    multi1d<LatticeColorMatrix> u_smr = u;
    {
      std::istringstream xml_l(link_smear.xml);
      XMLReader linktop(xml_l);
      Handle<LinkSmearing> linkSmearing(
          TheLinkSmearingFactory::Instance().createObject(
              link_smear.id, linktop, "/LinkSmearing"));
      (*linkSmearing)(u_smr);
    }
    auto computed = ns_getColorvecs::computeColorvecs(u_smr, local_t, 1,
                                                      n_colorvecs, order);
    Tensor<Nd + 3, ComplexD> colorvecs = computed.first;

    auto ip = contract<2>(r.conj(), colorvecs,
                          detail::remove_dimensions(colorvecs.order, "nt"),
                          OnHost, OnEveryoneReplicated)
                  .reorder("nt");
    auto phi = ip.like_this();
    for (int v = 0; v < n_colorvecs; ++v) {
      auto cv_norm = norms.get({v});
      auto phi_i = cv_norm * cv_norm / ip.get({v, 0});
      if (std::fabs(std::abs(phi_i) - 1) > 1e-4)
        throw std::runtime_error("the projector colorvec fingerprint does not "
                                 "correspond to the frozen links of "
                                 "this child at parent t_slice= " +
                                 std::to_string(global_t) + " vector " +
                                 std::to_string(v));
      phi.set({v, 0}, phi_i);
    }
    r.contract(colorvecs, {}, NotConjugate, phi, {}, NotConjugate);
  }

  // Only the master node read the data; share the checksum
  QDPInternal::broadcast(sum);
  if (checksum)
    *checksum = sum;

  return r.make_sure<SB::Complex>();
}
#endif

} // namespace FactorizedDistillation
} // namespace Chroma
