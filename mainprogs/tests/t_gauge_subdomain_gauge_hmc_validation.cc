#include "chroma.h"
#include "util/gauge/gauge_subdomain_split.h"
#include "util/gauge/gauge_startup.h"
#include "io/cfgtype_io.h"

#include <cmath>
#include <sstream>

using namespace Chroma;

namespace
{
  struct ValidationInput
  {
    std::string sidecar_file;
    Cfg_t child0_initial;
    Cfg_t child0_final;
    Cfg_t child1_initial;
    Cfg_t child1_final;
    double nontrivial_epsilon;
    double consistency_epsilon;

    ValidationInput()
      : nontrivial_epsilon(1.0e-8),
        consistency_epsilon(1.0e-12)
    {
    }
  };

  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_gauge_subdomain_gauge_hmc_validation: "
                << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  int latticeVolume(const multi1d<int>& nrow)
  {
    int vol = 1;
    for (int mu = 0; mu < nrow.size(); ++mu)
      vol *= nrow[mu];
    return vol;
  }

  int linkIndex(int mu, int site, int vol)
  {
    return mu * vol + site;
  }

  bool isFrozenLocalTime(int local_t,
                         const multi1d<GaugeSubdomainSplitInterval>& intervals)
  {
    for (int i = 0; i < intervals.size(); ++i)
    {
      if (local_t >= intervals[i].t_start && local_t <= intervals[i].t_end)
        return true;
    }

    return false;
  }

  multi1d<int> shiftedCoord(const multi1d<int>& coord,
                            int dir,
                            int delta,
                            const multi1d<int>& nrow)
  {
    multi1d<int> shifted = coord;
    shifted[dir] = (shifted[dir] + delta) % nrow[dir];
    if (shifted[dir] < 0)
      shifted[dir] += nrow[dir];
    return shifted;
  }

  void readValidationInput(XMLReader& xml, const std::string& path, ValidationInput& input)
  {
    XMLReader input_xml(xml, path);

    XMLReader sidecar_xml(input_xml, "Sidecar");
    read(sidecar_xml, "file", input.sidecar_file);

    read(input_xml, "Child0Initial", input.child0_initial);
    read(input_xml, "Child0Final", input.child0_final);
    read(input_xml, "Child1Initial", input.child1_initial);
    read(input_xml, "Child1Final", input.child1_final);

    if (input_xml.count("Checks") != 0)
    {
      XMLReader checks_xml(input_xml, "Checks");
      if (checks_xml.count("nontrivial_epsilon") != 0)
        read(checks_xml, "nontrivial_epsilon", input.nontrivial_epsilon);
      if (checks_xml.count("consistency_epsilon") != 0)
        read(checks_xml, "consistency_epsilon", input.consistency_epsilon);
    }
  }

  void loadGaugeLinks(const multi1d<int>& nrow,
                      const Cfg_t& cfg,
                      multi1d<ColorMatrix>& links)
  {
    Layout::setLattSize(nrow);
    Layout::create();

    multi1d<LatticeColorMatrix> u;
    XMLReader gauge_file_xml;
    XMLReader gauge_xml;
    Cfg_t cfg_copy = cfg;
    gaugeStartup(gauge_file_xml, gauge_xml, u, cfg_copy);
    snapshotGaugeField(u, links);
  }

  void checkExpectedPlan(const GaugeSubdomainSplitPlan& plan)
  {
    multi1d<int> expected_parent(Nd);
    const int parent_arr[] = {4, 4, 4, 8};
    expected_parent = parent_arr;

    multi1d<int> expected_child(Nd);
    const int child_arr[] = {4, 4, 4, 5};
    expected_child = child_arr;

    const int expected_child0_map[] = {0, 1, 2, 3, 4};
    const int expected_child1_map[] = {4, 5, 6, 7, 0};

    check(plan.param.t_dir == 3, "expected t_dir = 3");
    check(plan.param.cut0 == 0, "expected cut0 = 0");
    check(plan.param.cut1 == 4, "expected cut1 = 4");
    check(plan.param.frozen_width == 1, "expected frozen_width = 1");
    check(plan.parent_nrow.size() == Nd, "parent_nrow size mismatch");
    check(plan.child0_nrow.size() == Nd, "child0_nrow size mismatch");
    check(plan.child1_nrow.size() == Nd, "child1_nrow size mismatch");

    for (int mu = 0; mu < Nd; ++mu)
    {
      check(plan.parent_nrow[mu] == expected_parent[mu],
            "unexpected parent lattice extent");
      check(plan.child0_nrow[mu] == expected_child[mu],
            "unexpected child0 lattice extent");
      check(plan.child1_nrow[mu] == expected_child[mu],
            "unexpected child1 lattice extent");
    }

    check(plan.child0_local_to_global_t.size() == 5,
          "unexpected child0 local_to_global_t size");
    check(plan.child1_local_to_global_t.size() == 5,
          "unexpected child1 local_to_global_t size");

    for (int i = 0; i < 5; ++i)
    {
      check(plan.child0_local_to_global_t[i] == expected_child0_map[i],
            "unexpected child0 local_to_global_t entry");
      check(plan.child1_local_to_global_t[i] == expected_child1_map[i],
            "unexpected child1 local_to_global_t entry");
    }

    check(plan.child0_frozen_local_intervals.size() == 2,
          "child0 should report two frozen local intervals");
    check(plan.child1_frozen_local_intervals.size() == 2,
          "child1 should report two frozen local intervals");

    for (int child = 0; child < 2; ++child)
    {
      const multi1d<GaugeSubdomainSplitInterval>& intervals =
        (child == 0) ? plan.child0_frozen_local_intervals
                     : plan.child1_frozen_local_intervals;

      check(intervals[0].t_start == 0 && intervals[0].t_end == 0,
            "first frozen interval should be [0,0]");
      check(intervals[1].t_start == 4 && intervals[1].t_end == 4,
            "second frozen interval should be [4,4]");
    }
  }

  int checkFrozenLinkInvariance(const std::string& label,
                                const multi1d<int>& nrow,
                                int t_dir,
                                const multi1d<GaugeSubdomainSplitInterval>& intervals,
                                const multi1d<ColorMatrix>& initial_links,
                                const multi1d<ColorMatrix>& final_links)
  {
    const int vol = latticeVolume(nrow);
    check(initial_links.size() == Nd * vol, label + ": initial link count mismatch");
    check(final_links.size() == Nd * vol, label + ": final link count mismatch");

    int checked_links = 0;
    for (int site = 0; site < vol; ++site)
    {
      multi1d<int> coord = crtesn(site, nrow);
      if (!isFrozenLocalTime(coord[t_dir], intervals))
        continue;

      for (int mu = 0; mu < Nd; ++mu)
      {
        const double diff =
          toDouble(norm2(final_links[linkIndex(mu, site, vol)] -
                         initial_links[linkIndex(mu, site, vol)]));
        if (diff != 0.0)
        {
          std::ostringstream os;
          os << label << ": frozen link changed at site " << site
             << " direction " << mu << " with norm2 diff " << diff;
          fail(os.str());
        }
        ++checked_links;
      }
    }

    return checked_links;
  }

  double computeActiveRootSpatialPlaquette(const std::string& label,
                                           const multi1d<int>& nrow,
                                           int t_dir,
                                           const multi1d<GaugeSubdomainSplitInterval>& intervals,
                                           const multi1d<ColorMatrix>& links)
  {
    const int vol = latticeVolume(nrow);
    check(links.size() == Nd * vol, label + ": link count mismatch");

    double sum_plaq = 0.0;
    int count = 0;

    for (int site = 0; site < vol; ++site)
    {
      multi1d<int> coord = crtesn(site, nrow);
      if (isFrozenLocalTime(coord[t_dir], intervals))
        continue;

      for (int mu = 0; mu < Nd; ++mu)
      {
        if (mu == t_dir)
          continue;

        for (int nu = 0; nu < mu; ++nu)
        {
          if (nu == t_dir)
            continue;

          const int site_plus_mu = local_site(shiftedCoord(coord, mu, 1, nrow), nrow);
          const int site_plus_nu = local_site(shiftedCoord(coord, nu, 1, nrow), nrow);

          const ColorMatrix plaq =
            links[linkIndex(mu, site, vol)] *
            links[linkIndex(nu, site_plus_mu, vol)] *
            adj(links[linkIndex(mu, site_plus_nu, vol)]) *
            adj(links[linkIndex(nu, site, vol)]);

          sum_plaq += toDouble(real(trace(plaq))) / double(Nc);
          ++count;
        }
      }
    }

    check(count > 0, label + ": no active rooted spatial plaquettes were counted");
    return sum_plaq / double(count);
  }
}

int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  int status = 0;

  try
  {
    XMLReader xml_in(Chroma::getXMLInputFileName());
    ValidationInput input;
    readValidationInput(xml_in, "/gauge_subdomain_gauge_hmc_validation", input);

    GaugeSubdomainSplitPlan plan;
    {
      XMLReader sidecar_in(input.sidecar_file);
      read(sidecar_in, "/GaugeSubdomainSplitInfo", plan);
    }

    checkExpectedPlan(plan);

    multi1d<ColorMatrix> child0_initial_links;
    multi1d<ColorMatrix> child0_final_links;
    multi1d<ColorMatrix> child1_initial_links;
    multi1d<ColorMatrix> child1_final_links;

    loadGaugeLinks(plan.child0_nrow, input.child0_initial, child0_initial_links);
    loadGaugeLinks(plan.child0_nrow, input.child0_final, child0_final_links);
    loadGaugeLinks(plan.child1_nrow, input.child1_initial, child1_initial_links);
    loadGaugeLinks(plan.child1_nrow, input.child1_final, child1_final_links);

    const int child0_frozen_checks =
      checkFrozenLinkInvariance("child0",
                                plan.child0_nrow,
                                plan.param.t_dir,
                                plan.child0_frozen_local_intervals,
                                child0_initial_links,
                                child0_final_links);
    const int child1_frozen_checks =
      checkFrozenLinkInvariance("child1",
                                plan.child1_nrow,
                                plan.param.t_dir,
                                plan.child1_frozen_local_intervals,
                                child1_initial_links,
                                child1_final_links);

    const double child0_active_plaq =
      computeActiveRootSpatialPlaquette("child0",
                                        plan.child0_nrow,
                                        plan.param.t_dir,
                                        plan.child0_frozen_local_intervals,
                                        child0_final_links);
    const double child1_active_plaq =
      computeActiveRootSpatialPlaquette("child1",
                                        plan.child1_nrow,
                                        plan.param.t_dir,
                                        plan.child1_frozen_local_intervals,
                                        child1_final_links);

    const double child0_shift = std::fabs(child0_active_plaq - 1.0);
    const double child1_shift = std::fabs(child1_active_plaq - 1.0);
    const double child_diff = std::fabs(child0_active_plaq - child1_active_plaq);

    check(child0_shift > input.nontrivial_epsilon,
          "child0 active-root spatial plaquette stayed too close to 1.0");
    check(child1_shift > input.nontrivial_epsilon,
          "child1 active-root spatial plaquette stayed too close to 1.0");
    check(child_diff <= input.consistency_epsilon,
          "child active-root spatial plaquettes disagree beyond tolerance");

    QDPIO::cout << "t_gauge_subdomain_gauge_hmc_validation: passed" << std::endl;
    QDPIO::cout << "  child0 frozen links checked: " << child0_frozen_checks << std::endl;
    QDPIO::cout << "  child1 frozen links checked: " << child1_frozen_checks << std::endl;
    QDPIO::cout << "  child0 active_root_spatial_plaquette: "
                << child0_active_plaq << std::endl;
    QDPIO::cout << "  child1 active_root_spatial_plaquette: "
                << child1_active_plaq << std::endl;
    QDPIO::cout << "  child plaquette difference: " << child_diff << std::endl;
  }
  catch (const std::string& e)
  {
    QDPIO::cerr << e << std::endl;
    status = 1;
  }
  catch (std::exception& e)
  {
    QDPIO::cerr << "t_gauge_subdomain_gauge_hmc_validation: standard exception: "
                << e.what() << std::endl;
    status = 1;
  }
  catch (...)
  {
    QDPIO::cerr << "t_gauge_subdomain_gauge_hmc_validation: unknown exception"
                << std::endl;
    status = 1;
  }

  Chroma::finalize();
  return status;
}
