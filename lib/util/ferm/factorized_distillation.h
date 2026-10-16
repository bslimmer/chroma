// -*- C++ -*-
/*! \file
 * \brief Shared helpers for factorized (domain-decomposed) distillation measurements
 *
 * These helpers were extracted from the parent-lattice prototype in
 * meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc so that the
 * prototype and the child-local factorized baryon tensor measurement
 * (meas/inline/hadron/inline_factorized_baryon_tensor_distillation_superb_w.cc)
 * share one implementation.
 *
 * See specs/hier/factorized_prop.md and specs/hier/factorized_baryon_tensors.md.
 */

#ifndef __factorized_distillation_h__
#define __factorized_distillation_h__

#include "chromabase.h"
#include "util/gauge/gauge_subdomain_split.h"

#ifdef BUILD_SB
#  include "util/ferm/mgproton.h"
#  include "util/ferm/superb_contractions.h"
#endif

#include <string>
#include <vector>

namespace Chroma
{
  namespace FactorizedDistillation
  {
    //! Child-local geometry of a split fixed-boundary lattice, read from the split sidecar
    struct FactorizedPropGeometry {
      bool enabled;
      std::string sidecar_file;
      int child_id;
      int t_dir;
      multi1d<int> child_nrow;
      multi1d<int> local_to_global_t;
      multi1d<GaugeSubdomainSplitInterval> frozen_local_intervals;
      multi1d<GaugeSubdomainSplitInterval> active_local_intervals;
      multi1d<int> gamma_star_local_slices; ///< frozen child-local slices
      multi1d<int> gamma_local_slices;      ///< active child-local slices
      GaugeSubdomainSplitPlan plan;

      FactorizedPropGeometry() : enabled(false), child_id(-1), t_dir(-1) {}
    };

    //! Expand a list of inclusive intervals into a list of slices
    void appendIntervalSlices(const multi1d<GaugeSubdomainSplitInterval>& intervals,
			      multi1d<int>& slices);

    //! Complement of the frozen intervals inside [0, child_extent)
    multi1d<GaugeSubdomainSplitInterval>
    buildActiveIntervals(const multi1d<GaugeSubdomainSplitInterval>& frozen_intervals,
			 int child_extent);

    //! Read the split sidecar and build the child-local geometry for `child_id`
    /*!
     * \param sidecar_file: split sidecar XML written by gauge_subdomain_split
     * \param child_id: 0 or 1
     * \param decay_dir: must match the sidecar t_dir
     */
    FactorizedPropGeometry loadFactorizedPropGeometry(const std::string& sidecar_file,
						      int child_id, int decay_dir);

    //! Write the geometry into the measurement XML output
    void writeFactorizedPropGeometry(XMLWriter& xml, const FactorizedPropGeometry& geometry);

    //! Parent-global time of the first slice of the two frozen blocks (cut_left, cut_right)
    void getCutLeftRight(const GaugeSubdomainSplitPlan& plan, int& cut_left, int& cut_right);

    //! Returns a copy of u with all t_dir links rooted on t_slice set to zero
    multi1d<LatticeColorMatrix> zeroTemporalLinksOnSlice(const multi1d<LatticeColorMatrix>& u,
							 int t_slice, int t_dir = Nd - 1);

    //! Cut-link field used for D_sep = diag(D00, D11, D00)
    /*!
     * Zero the t_dir links rooted on t_start - 1 and t_end (mod the lattice extent) of every
     * frozen local interval. On a child lattice this cuts the frozen|active interfaces and
     * the child's periodic wrap.
     */
    multi1d<LatticeColorMatrix> buildCutLinkField(const multi1d<LatticeColorMatrix>& u,
						  const FactorizedPropGeometry& geometry);

    //! Child-local slice t* whose forward t link is the parent link Lt-1 -> 0, or -1
    int findParentTemporalCrossing(const FactorizedPropGeometry& geometry);

    //! Multiply the t_dir links rooted on the parent-crossing slice by `phase`
    multi1d<LatticeColorMatrix> applyParentTemporalPhase(const multi1d<LatticeColorMatrix>& u,
							 const FactorizedPropGeometry& geometry,
							 const Real& phase);

#ifdef BUILD_SB
    //! Return n zero-initialized LatticeFermion's
    SB::MultipleLatticeFermions returnNLatticeFermions(int n);

    //! Zero every field except on the listed time-slices
    SB::MultipleLatticeFermions restrictToTimeslices(const SB::MultipleLatticeFermions& psi,
						     const multi1d<int>& keep_tslices, int t_dir);

    //! Return a vector of LatticeFermion from the given SB::Tensor
    /// \param chi: lattice color tensor, cxyztXn plus spin labels
    /// \param t_source: first time-slice in chi
    SB::MultipleLatticeFermions toLatticeFermions(const SB::Tensor<Nd + 5, SB::Complex>& chi,
						  int t_source);

    //! Return an SB::Tensor with order cxyztXns from a vector of LatticeFermion
    /// \param chi: vector of LatticeFermion
    /// \param t_source: first time-slice to return
    /// \param num_tslices: number of timeslices to return
    SB::Tensor<Nd + 4, SB::Complex> toSBTensor(const SB::MultipleLatticeFermions& chi,
					       int t_source, int num_tslices);

    //! y_i = D x_i with the full (unpreconditioned) operator returned by genLinOp
    void applyLinOp(const LinearOperator<LatticeFermion>& linop, SB::MultipleLatticeFermions& y,
		    const SB::MultipleLatticeFermions& x);

    //! Read `n_colorvecs` colorvecs of the parent-global slice `global_t` from a parent-lattice
    //! colorvec file and place them at child-local slice `local_t` of the current layout.
    /*!
     * The stock SB::getColorvecs indexes the file with the current layout's time; this helper
     * reads the parent-global slice instead. Spatial extents of the file must match the
     * current layout; the file time extent must match `parent_Lt` (s3t files only).
     *
     * Fingerprint s3t files (CREATE_COLORVECS_SUPERB with write_fingerprint = true) only store
     * a small corner of each vector. For those, the eigenvectors are recomputed from `u` on the
     * child-local slice (the frozen slices hold the parent's spatial links) and phased to match
     * the stored fingerprint, exactly as SB::getColorvecs does on the parent.
     *
     * \param u: gauge field of the current (child) lattice; used for fingerprint files only
     * \param checksum: (optional output) deterministic checksum of the stored data in natural
     *        ordering; it is the same on every child that reads the same parent slice
     * \return tensor with order cxyzXnt and a single time-slice
     */
    SB::Tensor<Nd + 3, SB::Complex> readColorvecsOnSliceFromParentFile(
      const SB::ColorvecsStorage& sto, const multi1d<LatticeColorMatrix>& u, int global_t,
      int local_t, int n_colorvecs, int parent_Lt, double* checksum = nullptr);
#endif

  } // namespace FactorizedDistillation
} // namespace Chroma

#endif
