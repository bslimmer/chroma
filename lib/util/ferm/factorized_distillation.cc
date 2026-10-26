/*! \file
 * \brief Shared helpers for factorized (domain-decomposed) distillation
 * measurements
 */

#include "util/ferm/factorized_distillation.h"

#ifdef BUILD_SB
#include "io/xml_group_reader.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/key_timeslice_colorvec.h"
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
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

//----------------------------------------------------------------------------
// Edge-SVD boundary projectors, see specs/hier/factorized_prop_svd_projectors.md
// (moved from inline_prop_and_matelem_distillation_superb_w.cc, unchanged)
//----------------------------------------------------------------------------
namespace {
using EdgeSVDReal = SB::Complex::value_type;

//! New zeroed slice tensor with n columns
static EdgeSVDSliceTensor newEdgeSVDSliceTensor(int n) {
  EdgeSVDSliceTensor r(
      "cxyztXns",
      SB::latticeSize<Nd + 4>("cxyztXns", {{'t', 1}, {'n', n}}),
      SB::OnDefaultDevice, SB::OnEveryone);
  r.set_zero();
  return r;
}

//! Gram matrix (A^dag B)(m,n), contracting colour, space and spin
static SB::Tensor<2, SB::Complex> edgeSVDGram(const EdgeSVDSliceTensor &A,
                                               const EdgeSVDSliceTensor &B) {
  return SB::contract<2>(A.conj().rename_dims({{'n', 'm'}}), B, "cxyztXs");
}

//! Eigenpairs of a Hermitian PSD Gram matrix C (labels "mn") with SB::svd
/*!
 * Returns the eigenvectors on host and replicated (labels "m^") and the
 * eigenvalues in descending order. The factorization is done once by superbblas
 * and then replicated, so every rank uses the same vectors.
 */
static void edgeSVDGramEig(const SB::Tensor<2, SB::Complex> &C,
                           SB::Tensor<2, SB::Complex> &u_h,
                           std::vector<double> &lambda) {
  auto svd_res = SB::svd<2, 1, 2>(C, "m", "n", '^');
  u_h = std::get<0>(svd_res).make_sure(SB::none, SB::OnHost,
                                       SB::OnEveryoneReplicated);
  auto s_h = std::get<1>(svd_res).make_sure(SB::none, SB::OnHost,
                                            SB::OnEveryoneReplicated);
  const int k = s_h.kvdim().at('^');
  lambda.resize(k);
  for (int j = 0; j < k; ++j)
    lambda[j] = (double)s_h.get({j});
}

//! Y * X, with X a k x r matrix stored column-major; the result has r columns
/*!
 * `C_like` is a small tensor produced by SB::contract from Y; X is placed with
 * the same device and distribution before contracting with Y.
 */
static EdgeSVDSliceTensor
edgeSVDTimesSmall(const EdgeSVDSliceTensor &Y,
                  const SB::Tensor<2, SB::Complex> &C_like,
                  const std::vector<SB::Complex> &X, int k, int r) {
  if (Y.kvdim().at('n') != k)
    throw std::runtime_error("edgeSVDTimesSmall: inconsistent column count");
  SB::Tensor<2, SB::Complex> X_h("nm", SB::Coor<2>{k, r}, SB::OnHost,
                                 SB::OnEveryoneReplicated);
  for (int j = 0; j < r; ++j)
    for (int i = 0; i < k; ++i)
      X_h.set({i, j}, X[i + k * j]);
  auto Xd = C_like.like_this<2>("nm", {{'n', k}, {'m', r}});
  X_h.copyTo(Xd);
  return SB::contract<Nd + 4>(Y, Xd, "n").rename_dims({{'m', 'n'}});
}

//! Orthonormalize the columns of Y (Cholesky-QR through the Gram matrix, twice)
/*!
 * Directions with eigenvalue below rel_tol * largest are dropped, so the result
 * can have fewer columns than Y.
 */
static EdgeSVDSliceTensor edgeSVDOrthonormalize(const EdgeSVDSliceTensor &Y,
                                                double rel_tol) {
  EdgeSVDSliceTensor Q = Y;
  for (int pass = 0; pass < 2; ++pass) {
    auto C = edgeSVDGram(Q, Q);
    SB::Tensor<2, SB::Complex> u_h;
    std::vector<double> lambda;
    edgeSVDGramEig(C, u_h, lambda);
    const int k = lambda.size();
    if (k == 0 || !(lambda[0] > 0))
      throw std::runtime_error("edge SVD: the sketch is zero");
    int r = 0;
    while (r < k && lambda[r] > rel_tol * lambda[0])
      ++r;
    std::vector<SB::Complex> X(k * r);
    for (int j = 0; j < r; ++j) {
      const SB::Complex sc((EdgeSVDReal)(1.0 / std::sqrt(lambda[j])));
      for (int i = 0; i < k; ++i)
        X[i + k * j] = u_h.get({i, j}) * sc;
    }
    Q = edgeSVDTimesSmall(Q, C, X, k, r);
  }
  return Q;
}

//! out[i] = op(X[i] placed on in_slices[i]) read on out_slices[i]
/*!
 * The op is block-diagonal over the frozen intervals, so column j of every X[i]
 * goes into the same field and all intervals share the solves. When the column
 * counts differ, the shorter sets repeat earlier columns (outputs dropped), so
 * no source is ever zero. Fields are built `batch` columns at a time.
 */
static std::vector<EdgeSVDSliceTensor>
edgeSVDApply(const EdgeSVDFermionMap &op, const std::vector<EdgeSVDSliceTensor> &X,
             const std::vector<int> &in_slices,
             const std::vector<int> &out_slices, int batch) {
  const int ni = X.size();
  std::vector<int> ncols(ni);
  int nmax = 0;
  for (int i = 0; i < ni; ++i) {
    ncols[i] = X[i].kvdim().at('n');
    nmax = std::max(nmax, ncols[i]);
  }
  std::vector<EdgeSVDSliceTensor> out(ni);
  for (int i = 0; i < ni; ++i)
    out[i] = newEdgeSVDSliceTensor(ncols[i]);

  for (int b0 = 0; b0 < nmax; b0 += batch) {
    const int nb = std::min(batch, nmax - b0);
    auto x = returnNLatticeFermions(nb);
    for (int j = 0; j < nb; ++j)
      for (int i = 0; i < ni; ++i)
        X[i]
            .kvslice_from_size({{'n', (b0 + j) % ncols[i]}}, {{'n', 1}})
            .copyTo(SB::asTensorView(*x[j]).kvslice_from_size(
                {{'t', in_slices[i]}}, {{'t', 1}}));
    auto y = op(x);
    x.clear();
    for (int i = 0; i < ni; ++i) {
      if (b0 >= ncols[i])
        continue;
      const int nk = std::min(nb, ncols[i] - b0);
      toSBTensor(SB::MultipleLatticeFermions(y.begin(), y.begin() + nk),
                 out_slices[i], 1)
          .copyTo(out[i].kvslice_from_size({{'n', b0}}, {{'n', nk}}));
    }
  }
  return out;
}

} // namespace

//! Edge-SVD projector basis
/*!
 * Top-`rank` right singular vectors, on the slices p_slices[i], of an operator
 * M that maps p_slices[i] -> q_slices[i] block-diagonally over the frozen
 * intervals. Randomized SVD (Halko-Martinsson-Tropp):
 *   Y = M^dag Omega (Omega Gaussian on q, rank + oversampling columns),
 *   power_its times Y <- M^dag M orth(Y), Q = orth(Y),
 *   Z = M Q, Z^dag Z = W Sigma^2 W^dag, phi = Q W[:, :rank].
 * Returns phi per p slice; fills the Ritz singular values and the max deviation
 * of phi^dag phi from the identity.
 */
EdgeSVDBasis computeEdgeSVDBasis(
    const EdgeSVDFermionMap &M, const EdgeSVDFermionMap &Mdag,
    const std::vector<int> &p_slices, const std::vector<int> &q_slices,
    int rank, int oversampling, int power_its, int seed_value, int batch,
    std::map<int, std::vector<double>> &sing_vals,
    std::map<int, double> &ortho_err) {
  const int ni = p_slices.size();
  const int k = rank + oversampling;
  const double rel_tol = 100 * std::numeric_limits<EdgeSVDReal>::epsilon();

  // Gaussian probes on the q slices; the global RNG state is restored afterwards
  std::vector<EdgeSVDSliceTensor> Omega(ni);
  for (int i = 0; i < ni; ++i)
    Omega[i] = newEdgeSVDSliceTensor(k);
  {
    Seed saved_seed;
    QDP::RNG::savern(saved_seed);
    Seed svd_seed = seed_value;
    QDP::RNG::setrn(svd_seed);
    for (int j = 0; j < k; ++j) {
      auto eta = returnNLatticeFermions(1);
      gaussian(*eta[0]);
      for (int i = 0; i < ni; ++i)
        toSBTensor(eta, q_slices[i], 1)
            .copyTo(Omega[i].kvslice_from_size({{'n', j}}, {{'n', 1}}));
    }
    QDP::RNG::setrn(saved_seed);
  }

  // Range sketch on the p slices
  auto Y = edgeSVDApply(Mdag, Omega, q_slices, p_slices, batch);
  Omega.clear();
  for (int it = 0; it < power_its; ++it) {
    for (int i = 0; i < ni; ++i)
      Y[i] = edgeSVDOrthonormalize(Y[i], rel_tol);
    auto Z = edgeSVDApply(M, Y, p_slices, q_slices, batch);
    Y = edgeSVDApply(Mdag, Z, q_slices, p_slices, batch);
  }
  std::vector<EdgeSVDSliceTensor> Q(ni);
  for (int i = 0; i < ni; ++i) {
    Q[i] = edgeSVDOrthonormalize(Y[i], rel_tol);
    if (Q[i].kvdim().at('n') < rank)
      throw std::runtime_error(
          "edge SVD: the sketch on slice " + std::to_string(p_slices[i]) +
          " has numerical rank " + std::to_string(Q[i].kvdim().at('n')) +
          " < num_intermediate_projectors = " + std::to_string(rank));
  }
  Y.clear();

  // Rayleigh-Ritz
  auto Z = edgeSVDApply(M, Q, p_slices, q_slices, batch);
  EdgeSVDBasis basis;
  for (int i = 0; i < ni; ++i) {
    auto C = edgeSVDGram(Z[i], Z[i]);
    SB::Tensor<2, SB::Complex> u_h;
    std::vector<double> lambda;
    edgeSVDGramEig(C, u_h, lambda);
    const int kq = lambda.size();
    std::vector<SB::Complex> X(kq * rank);
    for (int j = 0; j < rank; ++j)
      for (int l = 0; l < kq; ++l)
        X[l + kq * j] = u_h.get({l, j});
    auto phi = edgeSVDTimesSmall(Q[i], C, X, kq, rank);

    std::vector<double> sv(kq);
    for (int j = 0; j < kq; ++j)
      sv[j] = std::sqrt(std::max(lambda[j], 0.0));
    sing_vals[p_slices[i]] = sv;

    auto G = edgeSVDGram(phi, phi).make_sure(SB::none, SB::OnHost,
                                             SB::OnEveryoneReplicated);
    double err = 0;
    for (int a = 0; a < rank; ++a)
      for (int b = 0; b < rank; ++b)
        err = std::max(err, (double)std::abs(G.get({a, b}) -
                                             SB::Complex(a == b ? 1 : 0)));
    ortho_err[p_slices[i]] = err;

    basis[p_slices[i]] = phi;
  }
  return basis;
}


EdgeSVDSlices edgeSVDSlices(const FactorizedPropGeometry &geometry,
                            bool overlapping, int Lt) {
  const multi1d<GaugeSubdomainSplitInterval> &fi =
      geometry.plan.child0_frozen_local_intervals;
  if (fi.size() != 2)
    throw std::runtime_error("edge SVD: expected 2 frozen intervals");
  if (fi[0].t_start != 0)
    throw std::runtime_error(
        "edge SVD: the first frozen interval must start at t = 0");
  const int fw0 = fi[0].t_end - fi[0].t_start + 1;
  const int fw1 = fi[1].t_end - fi[1].t_start + 1;
  if (fw0 < 3 || fw1 < 3)
    throw std::runtime_error("edge SVD: frozen width must be >= 3");

  EdgeSVDSlices s;
  s.frozen_width = fw0;
  s.cut1 = fi[1].t_start;
  if (!overlapping) {
    // M = T^dag; p = first sink slices of the non-overlapping chain
    s.p = {fi[1].t_end, fi[0].t_start};
    s.a = {fi[1].t_end - 1, fi[0].t_start + 1}; // e_p
    s.b = {fi[1].t_start + 1, fi[0].t_end - 1}; // e_s
    s.q = {fi[1].t_start, fi[0].t_end};
  } else {
    // M = G_F; same layout as the overlapping prop chain
    const int fw = fw0, c1 = fi[1].t_start;
    if (fw1 != fw)
      throw std::runtime_error(
          "edge SVD: frozen intervals have different widths");
    if (c1 <= fw || c1 + fw > Lt - 1)
      throw std::runtime_error(
          "edge SVD: empty source or sink active region");
    s.p = {fw - 1, c1};          // w2
    s.a = {fw - 2, c1 + 1};      // w1
    s.q = {0, c1 + fw - 1};      // sink-only buffer slices
    s.slab_cut_links = {c1, c1 + fw - 1, fw - 2, Lt - 1};
  }
  return s;
}

multi1d<LatticeColorMatrix>
nonOverlapCutField(const multi1d<LatticeColorMatrix> &u,
                   const FactorizedPropGeometry &geometry, int Lt) {
  multi1d<LatticeColorMatrix> u_sep = u;
  const multi1d<GaugeSubdomainSplitInterval> &fi =
      geometry.plan.child0_frozen_local_intervals;
  // D00 one slice in from each frozen boundary, so the clover term does not
  // cross into the active regions (needs a frozen width of at least 3)
  for (int i = 0; i < fi.size(); i++) {
    u_sep = zeroTemporalLinksOnSlice(u_sep,
                                     SB::normalize_coor(fi[i].t_start, Lt));
    u_sep = zeroTemporalLinksOnSlice(u_sep, fi[i].t_end - 1);
  }
  return u_sep;
}

namespace {
const std::string edge_svd_order = "cxyzXsnt";

std::string edgeSVDCanonicalXML(const std::string &s) {
  std::string r;
  for (char c : s)
    if (!std::isspace(static_cast<unsigned char>(c)))
      r.push_back(c);
  return r;
}

template <typename T> multi1d<T> toMulti1d(const std::vector<T> &v) {
  multi1d<T> r(v.size());
  for (int i = 0; i < v.size(); ++i)
    r[i] = v[i];
  return r;
}

template <typename T> std::vector<T> toStdVector(const multi1d<T> &v) {
  std::vector<T> r(v.size());
  for (int i = 0; i < v.size(); ++i)
    r[i] = v[i];
  return r;
}

SB::StorageTensor<Nd + 4, SB::ComplexD>
openEdgeSVDBasisStorage(const std::string &file) {
  return SB::StorageTensor<Nd + 4, SB::ComplexD>(
      file, true, std::string("/EdgeSVDBasis/order"));
}
} // namespace

void writeEdgeSVDBasis(const std::string &file, const EdgeSVDBasis &basis,
                       const EdgeSVDBasisInfo &info, XMLBufferWriter &config_xml,
                       const ChromaProp_t &prop) {
  const int ns = info.projector_slices.size();
  if (info.far_slices.size() != ns || info.singular_values.size() != ns ||
      info.orthonormality_error.size() != ns)
    throw std::runtime_error("writeEdgeSVDBasis: inconsistent metadata");
  const int nsv = ns > 0 ? info.singular_values[0].size() : 0;
  multi1d<double> sv_flat(ns * nsv);
  for (int i = 0; i < ns; ++i) {
    if (info.singular_values[i].size() != nsv)
      throw std::runtime_error(
          "writeEdgeSVDBasis: different number of singular values per slice");
    for (int j = 0; j < nsv; ++j)
      sv_flat[i * nsv + j] = info.singular_values[i][j];
  }

  XMLBufferWriter xml;
  push(xml, "EdgeSVDBasis");
  write(xml, "order", edge_svd_order);
  write(xml, "version", info.version);
  write(xml, "lattSize", info.latt_size);
  write(xml, "decay_dir", info.decay_dir);
  write(xml, "overlapping_factorization", info.overlapping);
  write(xml, "frozen_width", info.frozen_width);
  write(xml, "cut1", info.cut1);
  write(xml, "projector_slices", toMulti1d(info.projector_slices));
  write(xml, "far_slices", toMulti1d(info.far_slices));
  write(xml, "rank", info.rank);
  write(xml, "oversampling", info.oversampling);
  write(xml, "power_iterations", info.power_iterations);
  write(xml, "seed", info.seed);
  write(xml, "num_singular_values", nsv);
  write(xml, "singular_values", sv_flat); // slice-major, num_singular_values each
  write(xml, "orthonormality_error", toMulti1d(info.orthonormality_error));
  write(xml, "Propagator", prop);
  write(xml, "Config_info", config_xml);
  pop(xml);

  // NOTE: xml only has a valid value on the master node; so do a broadcast
  const std::string metadata = SB::broadcast(xml.str());

  SB::StorageTensor<Nd + 4, SB::ComplexD> st(
      file, metadata, edge_svd_order,
      SB::latticeSize<Nd + 4>(edge_svd_order, {{'n', info.rank}}), SB::Sparse,
      SB::checksum_type::BlockChecksum, SB::SharedFSFile);
  for (int p : info.projector_slices) {
    const EdgeSVDSliceTensor &phi = basis.at(p);
    if (phi.kvdim().at('n') != info.rank)
      throw std::runtime_error("writeEdgeSVDBasis: basis rank does not match");
    st.kvslice_from_size({{'t', p}}, {{'t', 1}}).copyFrom(phi);
  }
}

EdgeSVDBasisInfo readEdgeSVDBasisInfo(const std::string &file) {
  auto st = openEdgeSVDBasisStorage(file);
  std::istringstream is(st.metadata);
  XMLReader md(is);

  EdgeSVDBasisInfo info;
  read(md, "/EdgeSVDBasis/version", info.version);
  read(md, "/EdgeSVDBasis/lattSize", info.latt_size);
  read(md, "/EdgeSVDBasis/decay_dir", info.decay_dir);
  read(md, "/EdgeSVDBasis/overlapping_factorization", info.overlapping);
  read(md, "/EdgeSVDBasis/frozen_width", info.frozen_width);
  read(md, "/EdgeSVDBasis/cut1", info.cut1);
  multi1d<int> ps, fs;
  read(md, "/EdgeSVDBasis/projector_slices", ps);
  read(md, "/EdgeSVDBasis/far_slices", fs);
  info.projector_slices = toStdVector(ps);
  info.far_slices = toStdVector(fs);
  read(md, "/EdgeSVDBasis/rank", info.rank);
  read(md, "/EdgeSVDBasis/oversampling", info.oversampling);
  read(md, "/EdgeSVDBasis/power_iterations", info.power_iterations);
  read(md, "/EdgeSVDBasis/seed", info.seed);
  int nsv = 0;
  read(md, "/EdgeSVDBasis/num_singular_values", nsv);
  multi1d<double> sv_flat, oe;
  read(md, "/EdgeSVDBasis/singular_values", sv_flat);
  read(md, "/EdgeSVDBasis/orthonormality_error", oe);
  const int ns = info.projector_slices.size();
  if (sv_flat.size() != ns * nsv || oe.size() != ns)
    throw std::runtime_error("readEdgeSVDBasisInfo: inconsistent metadata in " +
                             file);
  info.singular_values.assign(ns, std::vector<double>(nsv));
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j < nsv; ++j)
      info.singular_values[i][j] = sv_flat[i * nsv + j];
  info.orthonormality_error = toStdVector(oe);
  info.fermact_xml =
      readXMLGroup(md, "/EdgeSVDBasis/Propagator/FermionAction", "FermAct")
          .xml;
  return info;
}

void checkEdgeSVDBasisInfo(const EdgeSVDBasisInfo &info,
                           const EdgeSVDSlices &slices, bool overlapping,
                           int decay_dir, int rank,
                           const std::string &fermact_xml) {
  auto fail = [](const std::string &what) {
    throw std::runtime_error("edge-SVD basis file does not match this run: " +
                             what);
  };
  const multi1d<int> &ls = Layout::lattSize();
  if (info.latt_size.size() != ls.size())
    fail("lattice size");
  for (int mu = 0; mu < ls.size(); ++mu)
    if (info.latt_size[mu] != ls[mu])
      fail("lattice size");
  if (info.decay_dir != decay_dir)
    fail("decay_dir");
  if (info.overlapping != overlapping)
    fail("overlapping_factorization (file " +
         std::string(info.overlapping ? "true" : "false") + ")");
  if (info.frozen_width != slices.frozen_width || info.cut1 != slices.cut1)
    fail("frozen geometry (file fw = " + std::to_string(info.frozen_width) +
         ", c1 = " + std::to_string(info.cut1) + "; run fw = " +
         std::to_string(slices.frozen_width) +
         ", c1 = " + std::to_string(slices.cut1) + ")");
  if (info.projector_slices != slices.p)
    fail("projector slices");
  if (info.rank < rank)
    fail("the file stores " + std::to_string(info.rank) +
         " vectors per slice, fewer than num_intermediate_projectors = " +
         std::to_string(rank));
  if (edgeSVDCanonicalXML(info.fermact_xml) != edgeSVDCanonicalXML(fermact_xml))
    fail("FermionAction differs");
}

EdgeSVDBasis readEdgeSVDBasis(const std::string &file,
                              const std::vector<int> &slices, int rank) {
  auto st = openEdgeSVDBasisStorage(file);

  // The file must use this build's layout (even-odd x/X split) and lattice
  const auto fdim = st.kvdim();
  for (char c : std::string("cxyzXst")) {
    const int expected = SB::latticeSize<1>(std::string(1, c))[0];
    if (fdim.at(c) != expected)
      throw std::runtime_error(
          "readEdgeSVDBasis: dimension `" + std::string(1, c) + "' of " + file +
          " is " + std::to_string(fdim.at(c)) + ", expected " +
          std::to_string(expected));
  }
  if (fdim.at('n') < rank)
    throw std::runtime_error("readEdgeSVDBasis: " + file + " stores only " +
                             std::to_string(fdim.at('n')) + " vectors");

  EdgeSVDBasis basis;
  for (int p : slices) {
    EdgeSVDSliceTensor phi = newEdgeSVDSliceTensor(rank);
    st.kvslice_from_size({{'t', p}, {'n', 0}}, {{'t', 1}, {'n', rank}})
        .copyTo(phi);

    // Sanity check: the vectors read must be orthonormal
    auto G = edgeSVDGram(phi, phi).make_sure(SB::none, SB::OnHost,
                                             SB::OnEveryoneReplicated);
    double err = 0;
    for (int a = 0; a < rank; ++a)
      for (int b = 0; b < rank; ++b)
        err = std::max(err, (double)std::abs(G.get({a, b}) -
                                             SB::Complex(a == b ? 1 : 0)));
    if (err > 1e-5)
      throw std::runtime_error("readEdgeSVDBasis: the vectors of slice " +
                               std::to_string(p) + " in " + file +
                               " are not orthonormal (max |phi^dag phi - 1| = " +
                               std::to_string(err) + ")");
    basis[p] = phi;
  }
  return basis;
}
#endif

} // namespace FactorizedDistillation
} // namespace Chroma
