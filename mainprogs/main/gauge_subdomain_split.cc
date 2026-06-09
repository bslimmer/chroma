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

  struct SplitInput
  {
    multi1d<int> parent_nrow;
    Cfg_t parent_cfg;
    GaugeSubdomainSplitParams param;
    GaugeSubdomainQIOOutput child0_out;
    GaugeSubdomainQIOOutput child1_out;
    std::string sidecar_file;
  };

  void readQIOOutput(XMLReader& xml, const std::string& path, GaugeSubdomainQIOOutput& p)
  {
    XMLReader out_xml(xml, path);
    read(out_xml, "cfg_file", p.cfg_file);
    if (out_xml.count("volfmt") != 0)
      read(out_xml, "volfmt", p.volfmt);
  }

  void readSplitInput(XMLReader& xml, const std::string& path, SplitInput& p)
  {
    XMLReader input_xml(xml, path);

    XMLReader parent_xml(input_xml, "Parent");
    read(parent_xml, "nrow", p.parent_nrow);
    read(parent_xml, "Cfg", p.parent_cfg);

    read(input_xml, "Param", p.param);
    readQIOOutput(input_xml, "Child0", p.child0_out);
    readQIOOutput(input_xml, "Child1", p.child1_out);

    XMLReader sidecar_xml(input_xml, "Sidecar");
    read(sidecar_xml, "file", p.sidecar_file);
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
    SplitInput input;
    readSplitInput(xml_in, "/gauge_subdomain_split", input);

    QDPIO::cout << "gauge_subdomain_split: reading parent config "
                << input.parent_cfg.cfg_file << std::endl;

    GaugeSubdomainSplitPlan plan = makeGaugeSubdomainSplitPlan(input.parent_nrow, input.param);

    multi1d<ColorMatrix> parent_links;
    {
      Layout::setLattSize(input.parent_nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> parent_u;
      XMLReader gauge_file_xml;
      XMLReader gauge_xml;
      Cfg_t cfg = input.parent_cfg;
      gaugeStartup(gauge_file_xml, gauge_xml, parent_u, cfg);
      snapshotGaugeField(parent_u, parent_links);
    }

    multi1d<ColorMatrix> child0_links;
    multi1d<ColorMatrix> child1_links;
    extractGaugeSubdomain(parent_links, plan, 0, child0_links);
    extractGaugeSubdomain(parent_links, plan, 1, child1_links);

    {
      XMLFileWriter sidecar_out(input.sidecar_file);
      write(sidecar_out, "GaugeSubdomainSplitInfo", plan);
      sidecar_out.close();
    }

    {
      Layout::setLattSize(plan.child0_nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> child0_u;
      materializeGaugeField(plan.child0_nrow, child0_links, child0_u);
      writeMinimalQIO(input.child0_out.cfg_file, input.child0_out.volfmt, child0_u);
    }

    {
      Layout::setLattSize(plan.child1_nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> child1_u;
      materializeGaugeField(plan.child1_nrow, child1_links, child1_u);
      writeMinimalQIO(input.child1_out.cfg_file, input.child1_out.volfmt, child1_u);
    }

    QDPIO::cout << "gauge_subdomain_split: wrote child configs" << std::endl;
    QDPIO::cout << "  child0: " << input.child0_out.cfg_file << std::endl;
    QDPIO::cout << "  child1: " << input.child1_out.cfg_file << std::endl;
    QDPIO::cout << "  sidecar: " << input.sidecar_file << std::endl;
  }
  catch (const std::string& e)
  {
    QDPIO::cerr << e << std::endl;
    status = 1;
  }
  catch (std::exception& e)
  {
    QDPIO::cerr << "gauge_subdomain_split: standard exception: "
                << e.what() << std::endl;
    status = 1;
  }
  catch (...)
  {
    QDPIO::cerr << "gauge_subdomain_split: unknown exception" << std::endl;
    status = 1;
  }

  Chroma::finalize();
  return status;
}
