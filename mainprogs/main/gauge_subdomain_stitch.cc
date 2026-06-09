#include "chroma.h"
#include "util/gauge/gauge_subdomain_split.h"
#include "util/gauge/gauge_startup.h"
#include "io/cfgtype_io.h"
#include "io/gauge_io.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"

using namespace Chroma;

namespace
{
  struct GaugeSubdomainQIOOutput
  {
    GaugeSubdomainQIOOutput() : volfmt(QDPIO_SINGLEFILE) {}

    std::string cfg_file;
    QDP_volfmt_t volfmt;
  };

  struct StitchInput
  {
    std::string sidecar_file;
    Cfg_t child0_cfg;
    Cfg_t child1_cfg;
    GaugeSubdomainQIOOutput parent_out;
  };

  void readQIOOutput(XMLReader& xml, const std::string& path, GaugeSubdomainQIOOutput& p)
  {
    XMLReader out_xml(xml, path);
    read(out_xml, "cfg_file", p.cfg_file);
    if (out_xml.count("volfmt") != 0)
      read(out_xml, "volfmt", p.volfmt);
  }

  void readStitchInput(XMLReader& xml, const std::string& path, StitchInput& p)
  {
    XMLReader input_xml(xml, path);

    XMLReader sidecar_xml(input_xml, "Sidecar");
    read(sidecar_xml, "file", p.sidecar_file);

    read(input_xml, "Child0Cfg", p.child0_cfg);
    read(input_xml, "Child1Cfg", p.child1_cfg);
    readQIOOutput(input_xml, "ParentOut", p.parent_out);
  }

  void writeMinimalQIO(const std::string& cfg_file,
                       QDP_volfmt_t volfmt,
                       const multi1d<LatticeColorMatrix>& u)
  {
    XMLBufferWriter file_xml;
    XMLBufferWriter record_xml;

    push(file_xml, "gauge");
    write(file_xml, "id", int(0));
    pop(file_xml);

    push(record_xml, "gauge");
    pop(record_xml);

    writeGauge(file_xml, record_xml, u, cfg_file, volfmt, QDPIO_SERIAL);
  }
}

int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  int status = 0;

  try
  {
    XMLReader xml_in(Chroma::getXMLInputFileName());
    StitchInput input;
    readStitchInput(xml_in, "/gauge_subdomain_stitch", input);

    GaugeSubdomainSplitPlan plan;
    {
      XMLReader sidecar_in(input.sidecar_file);
      read(sidecar_in, "/GaugeSubdomainSplitInfo", plan);
    }

    multi1d<ColorMatrix> child0_links;
    {
      Layout::setLattSize(plan.child0_nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> child0_u;
      XMLReader gauge_file_xml;
      XMLReader gauge_xml;
      Cfg_t cfg = input.child0_cfg;
      gaugeStartup(gauge_file_xml, gauge_xml, child0_u, cfg);
      snapshotGaugeField(child0_u, child0_links);
    }

    multi1d<ColorMatrix> child1_links;
    {
      Layout::setLattSize(plan.child1_nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> child1_u;
      XMLReader gauge_file_xml;
      XMLReader gauge_xml;
      Cfg_t cfg = input.child1_cfg;
      gaugeStartup(gauge_file_xml, gauge_xml, child1_u, cfg);
      snapshotGaugeField(child1_u, child1_links);
    }

    multi1d<ColorMatrix> parent_links;
    stitchGaugeSubdomains(plan, child0_links, child1_links, parent_links);

    {
      Layout::setLattSize(plan.parent_nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> parent_u;
      materializeGaugeField(plan.parent_nrow, parent_links, parent_u);
      writeMinimalQIO(input.parent_out.cfg_file, input.parent_out.volfmt, parent_u);
    }

    QDPIO::cout << "gauge_subdomain_stitch: wrote reconstructed parent config "
                << input.parent_out.cfg_file << std::endl;
  }
  catch (const std::string& e)
  {
    QDPIO::cerr << e << std::endl;
    status = 1;
  }
  catch (std::exception& e)
  {
    QDPIO::cerr << "gauge_subdomain_stitch: standard exception: "
                << e.what() << std::endl;
    status = 1;
  }
  catch (...)
  {
    QDPIO::cerr << "gauge_subdomain_stitch: unknown exception" << std::endl;
    status = 1;
  }

  Chroma::finalize();
  return status;
}
