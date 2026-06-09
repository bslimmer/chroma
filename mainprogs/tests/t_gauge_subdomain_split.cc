#include "chroma.h"
#include "util/gauge/gauge_subdomain_split.h"
#include "qdp_util.h"

#include <sstream>

using namespace Chroma;

namespace
{
  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_gauge_subdomain_split: " << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  int linkIndex(int mu, int site, int vol)
  {
    return mu * vol + site;
  }

  int latticeVolume(const multi1d<int>& nrow)
  {
    int vol = 1;
    for (int mu = 0; mu < nrow.size(); ++mu)
      vol *= nrow[mu];
    return vol;
  }

  void buildDeterministicGaugeField(const multi1d<int>& nrow,
                                    multi1d<LatticeColorMatrix>& u)
  {
    Layout::setLattSize(nrow);
    Layout::create();

    u.resize(Nd);
    for (int mu = 0; mu < Nd; ++mu)
      u[mu] = zero;

    for (int site = 0; site < Layout::vol(); ++site)
    {
      multi1d<int> coord = crtesn(site, nrow);

      for (int mu = 0; mu < Nd; ++mu)
      {
        ColorMatrix cm = zero;
        Real base = Real(site + 10 * mu + 1);
        for (int c = 0; c < Nc; ++c)
        {
          Complex value = cmplx(base + Real(c), base + Real(2 * c + 1));
          pokeColor(cm, value, c, c);
        }
        pokeSite(u[mu], cm, coord);
      }
    }
  }

  void checkPlanXMLRoundTrip(const GaugeSubdomainSplitPlan& plan)
  {
    XMLBufferWriter xml;
    write(xml, "GaugeSubdomainSplitInfo", plan);

    std::istringstream is(xml.str());
    XMLReader xml_in(is);

    GaugeSubdomainSplitPlan parsed;
    read(xml_in, "/GaugeSubdomainSplitInfo", parsed);

    check(parsed.param.t_dir == plan.param.t_dir, "plan XML round-trip lost t_dir");
    check(parsed.param.cut0 == plan.param.cut0, "plan XML round-trip lost cut0");
    check(parsed.param.cut1 == plan.param.cut1, "plan XML round-trip lost cut1");
    check(parsed.param.frozen_width == plan.param.frozen_width,
          "plan XML round-trip lost frozen_width");
    check(parsed.child0_local_to_global_t.size() == plan.child0_local_to_global_t.size(),
          "plan XML round-trip changed child0 map length");

    for (int i = 0; i < plan.child0_local_to_global_t.size(); ++i)
    {
      check(parsed.child0_local_to_global_t[i] == plan.child0_local_to_global_t[i],
            "plan XML round-trip changed child0 local_to_global_t");
    }
  }

  void verifyNonWrappingPlan()
  {
    multi1d<int> nrow(Nd);
    const int arr[] = {2, 2, 2, 8};
    nrow = arr;

    GaugeSubdomainSplitParams param;
    param.t_dir = Nd - 1;
    param.cut0 = 1;
    param.cut1 = 5;
    param.frozen_width = 1;

    GaugeSubdomainSplitPlan plan = makeGaugeSubdomainSplitPlan(nrow, param);

    check(plan.child0_local_to_global_t.size() == 5, "child0 size mismatch in non-wrapping plan");
    check(plan.child1_local_to_global_t.size() == 5, "child1 size mismatch in non-wrapping plan");

    const int expect_child0[] = {1, 2, 3, 4, 5};
    const int expect_child1[] = {5, 6, 7, 0, 1};

    for (int i = 0; i < 5; ++i)
    {
      check(plan.child0_local_to_global_t[i] == expect_child0[i],
            "unexpected child0 local_to_global_t in non-wrapping plan");
      check(plan.child1_local_to_global_t[i] == expect_child1[i],
            "unexpected child1 local_to_global_t in non-wrapping plan");
    }

    check(plan.child0_frozen_local_intervals[0].t_start == 0 &&
          plan.child0_frozen_local_intervals[0].t_end == 0,
          "unexpected child0 first frozen interval");
    check(plan.child0_frozen_local_intervals[1].t_start == 4 &&
          plan.child0_frozen_local_intervals[1].t_end == 4,
          "unexpected child0 second frozen interval");

    checkPlanXMLRoundTrip(plan);
  }

  void verifyWrappingPlan()
  {
    multi1d<int> nrow(Nd);
    const int arr[] = {2, 2, 2, 16};
    nrow = arr;

    GaugeSubdomainSplitParams param;
    param.t_dir = Nd - 1;
    param.cut0 = 4;
    param.cut1 = 14;
    param.frozen_width = 3;

    GaugeSubdomainSplitPlan plan = makeGaugeSubdomainSplitPlan(nrow, param);

    check(plan.child1_local_to_global_t[0] == 14, "wrapping plan child1 did not start at wrapped frozen slice");
    check(plan.child1_local_to_global_t[1] == 15, "wrapping plan child1 lost wrapped slice 15");
    check(plan.child1_local_to_global_t[2] == 0, "wrapping plan child1 lost wrapped slice 0");
  }

  void verifyRoundTripAndOwnedUpdates()
  {
    multi1d<int> nrow(Nd);
    const int arr[] = {2, 2, 2, 8};
    nrow = arr;

    GaugeSubdomainSplitParams param;
    param.t_dir = Nd - 1;
    param.cut0 = 1;
    param.cut1 = 5;
    param.frozen_width = 1;

    GaugeSubdomainSplitPlan plan = makeGaugeSubdomainSplitPlan(nrow, param);

    multi1d<LatticeColorMatrix> parent_u;
    buildDeterministicGaugeField(nrow, parent_u);

    multi1d<ColorMatrix> parent_links;
    snapshotGaugeField(parent_u, parent_links);

    multi1d<ColorMatrix> child0_links;
    multi1d<ColorMatrix> child1_links;
    extractGaugeSubdomain(parent_links, plan, 0, child0_links);
    extractGaugeSubdomain(parent_links, plan, 1, child1_links);

    multi1d<ColorMatrix> stitched_links;
    stitchGaugeSubdomains(plan, child0_links, child1_links, stitched_links);

    check(stitched_links.size() == parent_links.size(), "stitched gauge buffer size mismatch");
    for (int i = 0; i < parent_links.size(); ++i)
    {
      if (toDouble(norm2(stitched_links[i] - parent_links[i])) != 0.0)
        fail("exact round-trip changed the parent gauge field");
    }

    multi1d<ColorMatrix> child0_mod = child0_links;
    multi1d<ColorMatrix> child1_mod = child1_links;

    ColorMatrix child0_marker = 2.0;
    ColorMatrix child1_marker = 3.0;

    multi1d<int> child0_coord(Nd);
    child0_coord = 0;
    child0_coord[param.t_dir] = 1;
    const int child0_site = local_site(child0_coord, plan.child0_nrow);
    for (int mu = 0; mu < Nd; ++mu)
      child0_mod[linkIndex(mu, child0_site, latticeVolume(plan.child0_nrow))] = child0_marker;

    multi1d<int> child1_coord(Nd);
    child1_coord = 0;
    child1_coord[param.t_dir] = 1;
    const int child1_site = local_site(child1_coord, plan.child1_nrow);
    for (int mu = 0; mu < Nd; ++mu)
      child1_mod[linkIndex(mu, child1_site, latticeVolume(plan.child1_nrow))] = child1_marker;

    stitchGaugeSubdomains(plan, child0_mod, child1_mod, stitched_links);

    multi1d<int> parent_coord0(Nd);
    parent_coord0 = 0;
    parent_coord0[param.t_dir] = plan.child0_local_to_global_t[child0_coord[param.t_dir]];
    const int parent_site0 = local_site(parent_coord0, plan.parent_nrow);

    multi1d<int> parent_coord1(Nd);
    parent_coord1 = 0;
    parent_coord1[param.t_dir] = plan.child1_local_to_global_t[child1_coord[param.t_dir]];
    const int parent_site1 = local_site(parent_coord1, plan.parent_nrow);

    for (int mu = 0; mu < Nd; ++mu)
    {
      check(toDouble(norm2(stitched_links[linkIndex(mu, parent_site0, latticeVolume(plan.parent_nrow))] -
                            child0_marker)) == 0.0,
            "owned child0 update did not reach the stitched parent");
      check(toDouble(norm2(stitched_links[linkIndex(mu, parent_site1, latticeVolume(plan.parent_nrow))] -
                            child1_marker)) == 0.0,
            "owned child1 update did not reach the stitched parent");
    }
  }

  void verifyMismatchAndInvalidCases()
  {
    multi1d<int> nrow(Nd);
    const int arr[] = {2, 2, 2, 8};
    nrow = arr;

    GaugeSubdomainSplitParams param;
    param.t_dir = Nd - 1;
    param.cut0 = 1;
    param.cut1 = 5;
    param.frozen_width = 1;

    GaugeSubdomainSplitPlan plan = makeGaugeSubdomainSplitPlan(nrow, param);

    multi1d<LatticeColorMatrix> parent_u;
    buildDeterministicGaugeField(nrow, parent_u);

    multi1d<ColorMatrix> parent_links;
    snapshotGaugeField(parent_u, parent_links);

    multi1d<ColorMatrix> child0_links;
    multi1d<ColorMatrix> child1_links;
    extractGaugeSubdomain(parent_links, plan, 0, child0_links);
    extractGaugeSubdomain(parent_links, plan, 1, child1_links);

    ColorMatrix mismatch_marker = 4.0;
    for (int mu = 0; mu < Nd; ++mu)
      child1_links[linkIndex(mu, 0, latticeVolume(plan.child1_nrow))] = mismatch_marker;

    bool mismatch_caught = false;
    try
    {
      multi1d<ColorMatrix> parent_out;
      stitchGaugeSubdomains(plan, child0_links, child1_links, parent_out);
    }
    catch (const std::string&)
    {
      mismatch_caught = true;
    }

    check(mismatch_caught, "stitchGaugeSubdomains did not reject a frozen-boundary mismatch");

    GaugeSubdomainSplitParams invalid = param;
    invalid.cut0 = invalid.cut1;

    bool invalid_caught = false;
    try
    {
      (void)makeGaugeSubdomainSplitPlan(nrow, invalid);
    }
    catch (const std::string&)
    {
      invalid_caught = true;
    }

    check(invalid_caught, "makeGaugeSubdomainSplitPlan did not reject cut0 == cut1");
  }
}

int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_gauge_subdomain_split");

  verifyNonWrappingPlan();
  verifyWrappingPlan();
  verifyRoundTripAndOwnedUpdates();
  verifyMismatchAndInvalidCases();

  pop(xml_out);
  xml_out.close();

  Chroma::finalize();
  return 0;
}
