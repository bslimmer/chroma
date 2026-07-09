// -*- C++ -*-
/*! \file
 *  \brief Construct 0++, 2++ and 1+- glueball correlation functions from fuzzy links
 */

#ifndef __gluecor_h__
#define __gluecor_h__

#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Construct 0++, 2++ and 1+- glueball correlation functions from fuzzy links
  /*! 
   * \ingroup glue
   *
   * Construct 0++, 2++ and 1+- glueball correlation functions from
   * fuzzy links at blocking level bl_level and write them in XML
   * format unless simple_output is enabled.
   *
   * Warning: this works only for Nd = 4 !
   *
   * \param xml_out       xml file object ( Write )
   * \param xml_group     std::string used for writing xml data, or a CSV
   *                      basename when simple_output is true ( Read )
   * \param u             (blocked) gauge field ( Read )
   * \param bl_level      blocking level ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param simple_output if true, write op0 to CSV instead of XML ( Read )
   * \param simple_output_file optional CSV path override for simple_output
   *                           mode ( Read )
   */

  void gluecor(XMLWriter& xml_out, const std::string& xml_group,
	       const multi1d<LatticeColorMatrix>& u, 
	       const SftMom& phases,
	       int bl_level,
	       bool simple_output = false,
	       const std::string& simple_output_file = "");

}  // end namespace Chroma

#endif
