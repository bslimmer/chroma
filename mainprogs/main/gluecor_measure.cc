/*! \file
 *  \brief Measure blocked gluecor op0[t] from one saved gauge configuration
 */

#include "chroma.h"
#include "io/cfgtype_io.h"
#include "meas/glue/block.h"
#include "meas/glue/gluecor.h"
#include "util/ft/sftmom.h"
#include "util/gauge/gauge_startup.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace Chroma;

namespace
{
  void usage(const char* prog)
  {
    std::cerr << "Usage: " << prog
              << " --nrow n0 n1 n2 n3 [--cfg-type type]"
              << " [--decay-dir dir] [--bl-level n]"
              << " [--blk-accu eps] [--blk-max n]"
              << " <cfg_file> <csv_file>\n";
  }

  void fail(const std::string& message)
  {
    throw std::string("gluecor_measure: " + message);
  }

  bool endsWith(const std::string& value, const std::string& suffix)
  {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  CfgType parseCfgType(const std::string& cfg_type_name)
  {
    if (cfg_type_name == "MILC")
      return CFG_TYPE_MILC;
    if (cfg_type_name == "NERSC")
      return CFG_TYPE_NERSC;
    if (cfg_type_name == "SCIDAC")
      return CFG_TYPE_SCIDAC;
    if (cfg_type_name == "SZIN")
      return CFG_TYPE_SZIN;
    if (cfg_type_name == "SZINQIO")
      return CFG_TYPE_SZINQIO;
    if (cfg_type_name == "KYU")
      return CFG_TYPE_KYU;
    if (cfg_type_name == "DISORDERED")
      return CFG_TYPE_DISORDERED;
    if (cfg_type_name == "UNIT")
      return CFG_TYPE_UNIT;
    if (cfg_type_name == "CPPACS")
      return CFG_TYPE_CPPACS;
    if (cfg_type_name == "CERN")
      return CFG_TYPE_CERN;
    if (cfg_type_name == "WEAK_FIELD")
      return CFG_TYPE_WEAK_FIELD;
    if (cfg_type_name == "CLASSICAL_SF")
      return CFG_TYPE_CLASSICAL_SF;
    if (cfg_type_name == "WUP")
      return CFG_TYPE_WUPP;

    fail("unsupported cfg type: " + cfg_type_name);
    return CFG_TYPE_SCIDAC;
  }

  CfgType inferCfgType(const std::string& cfg_file)
  {
    if (endsWith(cfg_file, ".scidac"))
      return CFG_TYPE_SCIDAC;

    if (endsWith(cfg_file, ".lime"))
      return CFG_TYPE_SZINQIO;

    fail("could not infer cfg type from file extension for " + cfg_file +
         "; pass --cfg-type explicitly");
    return CFG_TYPE_SCIDAC;
  }

  int computeMaxBlockingLevel(const multi1d<int>& nrow, int decay_dir)
  {
    int bl_level_max = nrow[0];

    for (int mu = 0; mu < nrow.size(); ++mu)
    {
      if (mu == decay_dir)
        continue;

      int block_latt = nrow[mu];
      int bl_level = 0;
      while ((block_latt > 2) && ((block_latt & 1) == 0))
      {
        block_latt /= 2;
        bl_level += 1;
      }

      if (bl_level < bl_level_max)
        bl_level_max = bl_level;
    }

    return bl_level_max;
  }

  multi1d<LatticeColorMatrix> buildBlockedGaugeField(const multi1d<LatticeColorMatrix>& u,
                                                     int decay_dir,
                                                     int bl_level_selected,
                                                     double blk_accu,
                                                     int blk_max)
  {
    const int max_block_level =
      computeMaxBlockingLevel(Layout::lattSize(), decay_dir);

    if (bl_level_selected > max_block_level)
    {
      std::ostringstream os;
      os << "requested blocking level " << bl_level_selected
         << " exceeds geometry-supported maximum " << max_block_level;
      fail(os.str());
    }

    multi1d<LatticeColorMatrix> u_fuz = u;
    multi1d<LatticeColorMatrix> u_tmp(Nd);

    for (int bl_level = 0; bl_level < bl_level_selected; ++bl_level)
    {
      u_tmp = u_fuz;
      for (int mu = 0; mu < Nd; ++mu)
      {
        if (mu != decay_dir)
          block(u_tmp[mu], u_fuz, mu, bl_level, Real(blk_accu), blk_max, decay_dir);
      }
      u_fuz = u_tmp;
    }

    return u_fuz;
  }
}

int main(int argc, char* argv[])
{
  Chroma::initialize(&argc, &argv);

  {
    int exit_code = 0;

    try
    {
      multi1d<int> nrow(Nd);
      bool have_nrow = false;

      std::string cfg_file;
      std::string csv_file;
      std::string cfg_type_name;
      bool have_cfg_type = false;
      bool help_requested = false;

      int decay_dir = 3;
      int bl_level = 1;
      int blk_max = 50;
      double blk_accu = 1.0e-5;

      for (int i = 1; i < argc; ++i)
      {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h")
        {
          help_requested = true;
          break;
        }
        else if (arg == "--nrow")
        {
          if (i + Nd >= argc)
          {
            usage(argv[0]);
            fail("--nrow expects four lattice extents");
          }

          for (int mu = 0; mu < Nd; ++mu)
            nrow[mu] = std::atoi(argv[++i]);

          have_nrow = true;
        }
        else if (arg == "--cfg-type")
        {
          if (i + 1 >= argc)
          {
            usage(argv[0]);
            fail("--cfg-type expects a value");
          }

          cfg_type_name = argv[++i];
          have_cfg_type = true;
        }
        else if (arg == "--decay-dir")
        {
          if (i + 1 >= argc)
          {
            usage(argv[0]);
            fail("--decay-dir expects a value");
          }

          decay_dir = std::atoi(argv[++i]);
        }
        else if (arg == "--bl-level")
        {
          if (i + 1 >= argc)
          {
            usage(argv[0]);
            fail("--bl-level expects a value");
          }

          bl_level = std::atoi(argv[++i]);
        }
        else if (arg == "--blk-accu")
        {
          if (i + 1 >= argc)
          {
            usage(argv[0]);
            fail("--blk-accu expects a value");
          }

          blk_accu = std::atof(argv[++i]);
        }
        else if (arg == "--blk-max")
        {
          if (i + 1 >= argc)
          {
            usage(argv[0]);
            fail("--blk-max expects a value");
          }

          blk_max = std::atoi(argv[++i]);
        }
        else if (cfg_file.empty())
        {
          cfg_file = arg;
        }
        else if (csv_file.empty())
        {
          csv_file = arg;
        }
        else
        {
          usage(argv[0]);
          fail("unexpected extra argument: " + arg);
        }
      }

      if (help_requested)
      {
        usage(argv[0]);
        Chroma::finalize();
        return 0;
      }

      if (!have_nrow || cfg_file.empty() || csv_file.empty())
      {
        usage(argv[0]);
        fail("missing required arguments");
      }

      if (decay_dir < 0 || decay_dir >= Nd)
        fail("decay_dir must be in [0, Nd)");

      if (bl_level < 0)
        fail("bl-level must be >= 0");

      if (blk_max <= 0)
        fail("blk-max must be > 0");

      Cfg_t cfg;
      cfg.cfg_type = have_cfg_type ? parseCfgType(cfg_type_name) : inferCfgType(cfg_file);
      cfg.cfg_file = cfg_file;

      Layout::setLattSize(nrow);
      Layout::create();

      multi1d<LatticeColorMatrix> u(Nd);
      XMLReader gauge_file_xml, gauge_xml;
      gaugeStartup(gauge_file_xml, gauge_xml, u, cfg);

      const multi1d<LatticeColorMatrix> u_blocked =
        buildBlockedGaugeField(u, decay_dir, bl_level, blk_accu, blk_max);

      SftMom phases(0, true, decay_dir);
      XMLBufferWriter scratch_xml;
      gluecor(scratch_xml,
              "GlueCorr",
              u_blocked,
              phases,
              bl_level,
              true,
              csv_file);

      QDPIO::cout << "gluecor_measure: wrote " << csv_file << std::endl;
      QDPIO::cout << "  cfg_file: " << cfg_file << std::endl;
      QDPIO::cout << "  decay_dir: " << decay_dir << std::endl;
      QDPIO::cout << "  bl_level: " << bl_level << std::endl;
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << e << std::endl;
      exit_code = 1;
    }
    catch (std::exception& e)
    {
      QDPIO::cerr << "gluecor_measure: standard exception: "
                  << e.what() << std::endl;
      exit_code = 1;
    }
    catch (...)
    {
      QDPIO::cerr << "gluecor_measure: unknown exception" << std::endl;
      exit_code = 1;
    }

    Chroma::finalize();
    return exit_code;
  }
}
