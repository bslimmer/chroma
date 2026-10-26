// -*- C++ -*-
/*! \file
 * \brief Edge-SVD boundary projector basis for the factorized perambulator
 *
 * Computes the top right singular vectors of the edge-to-edge block that carries a
 * quark across the frozen region (randomized SVD) and writes them to a file that
 * PROP_AND_MATELEM_DISTILLATION_SUPERB reads with projector_type = singular_vectors.
 * See specs/hier/factorized_prop_svd_projectors.md.
 */

#ifndef __inline_edge_svd_projectors_superb_w_h__
#define __inline_edge_svd_projectors_superb_w_h__

#include "chromabase.h"
#include "io/qprop_io.h"
#include "meas/inline/abs_inline_measurement.h"

#ifdef BUILD_SB

namespace Chroma
{
  /*! \ingroup inlinehadron */
  namespace InlineEdgeSVDProjectorsSuperbEnv
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long frequency;

      struct Param_t {
	int rank;			/*!< Vectors stored per projector slice (R_max) */
	int oversampling;		/*!< Extra random probes beyond the rank (default 16) */
	int power_iterations;		/*!< Power iterations of the randomized SVD (default 2) */
	int seed;			/*!< QDP RNG seed for the random probes (default 11) */
	bool overlapping_factorization; /*!< Overlapping (true) or non-overlapping path */
	int max_rhs;			/*!< Columns per batch of solves (default 8) */
	int decay_dir;			/*!< Decay direction (must be Nd-1) */
	ChromaProp_t prop;		/*!< FermionAction + InvertParam */
      };

      struct NamedObject_t {
	std::string gauge_id;			  /*!< Gauge field */
	std::string frozen_boundary_sidecar_file; /*!< Split sidecar */
	int frozen_boundary_child_id;		  /*!< Child id selecting the geometry */
	std::string svd_basis_file;		  /*!< Output basis file */
      };

      Param_t param;
      NamedObject_t named_obj;
      std::string xml_file; /*!< Alternate XML file pattern */
    };

    //! Inline measurement
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const { return params.frequency; }

      //! Do the measurement
      void operator()(const unsigned long update_no, XMLWriter& xml_out);

    protected:
      //! Do the measurement
      void func(const unsigned long update_no, XMLWriter& xml_out);

    private:
      Params params;
    };

  } // namespace InlineEdgeSVDProjectorsSuperbEnv
}

#endif // BUILD_SB

#endif
