#include "chroma.h"
#include "actions/gauge/gaugebcs/temporal_zone_gaugebc.h"
#include "update/molecdyn/hmc/gauge_monomial_momentum_bc.h"
#include "update/molecdyn/hmc/lcm_hmc.h"
#include "update/molecdyn/hmc/const_lcm_hmc.h"

#include <sstream>
#include <vector>

using namespace Chroma;

namespace
{
  typedef multi1d<LatticeColorMatrix> LCM;

  const Double tolerance = Double(1.0e-12);

  struct MonomialSpec
  {
    std::string monomial_name;
    std::string monomial_id;
    bool masked;
    TemporalZoneGaugeBCParams params;
  };

  class RefreshHarness : public LatColMatHMCTrj
  {
  public:
    RefreshHarness(const Handle<GaugeMonomial>& mask_source)
      : LatColMatHMCTrj(dummy_hamiltonian, dummy_integrator, mask_source)
    {
    }

    void run(AbsFieldState<LCM, LCM>& state) const
    {
      refreshP(state);
    }

  private:
    static Handle< AbsHamiltonian<LCM, LCM> > dummy_hamiltonian;
    static Handle< AbsMDIntegrator<LCM, LCM> > dummy_integrator;
  };

  class ConstRefreshHarness : public ConstLatColMatHMCTrj
  {
  public:
    ConstRefreshHarness(const Handle<GaugeMonomial>& mask_source)
      : ConstLatColMatHMCTrj(dummy_hamiltonian, dummy_integrator, mask_source)
    {
    }

    void run(AbsFieldState<LCM, LCM>& state) const
    {
      refreshP(state);
    }

  private:
    static Handle< AbsHamiltonian<LCM, LCM> > dummy_hamiltonian;
    static Handle< AbsMDIntegrator<LCM, LCM> > dummy_integrator;
  };

  Handle< AbsHamiltonian<LCM, LCM> > RefreshHarness::dummy_hamiltonian;
  Handle< AbsMDIntegrator<LCM, LCM> > RefreshHarness::dummy_integrator;
  Handle< AbsHamiltonian<LCM, LCM> > ConstRefreshHarness::dummy_hamiltonian;
  Handle< AbsMDIntegrator<LCM, LCM> > ConstRefreshHarness::dummy_integrator;

  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_hmc_momentum_bc_autodiscovery: " << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  MonomialSpec makePeriodicSpec(const std::string& monomial_name,
                                const std::string& monomial_id)
  {
    MonomialSpec spec;
    spec.monomial_name = monomial_name;
    spec.monomial_id = monomial_id;
    spec.masked = false;
    return spec;
  }

  MonomialSpec makeMaskedSpec(const std::string& monomial_name,
                              const std::string& monomial_id,
                              int t_start,
                              int t_end)
  {
    MonomialSpec spec;
    spec.monomial_name = monomial_name;
    spec.monomial_id = monomial_id;
    spec.masked = true;
    spec.params.zero_intervals.resize(1);
    spec.params.zero_intervals[0].t_start = t_start;
    spec.params.zero_intervals[0].t_end = t_end;
    return spec;
  }

  void appendGaugeMonomial(XMLWriter& xml_out, const MonomialSpec& spec)
  {
    push(xml_out, "elem");
    write(xml_out, "Name", spec.monomial_name);

    push(xml_out, "GaugeAction");
    write(xml_out, "Name", std::string("WILSON_GAUGEACT"));
    write(xml_out, "beta", Real(5.7));
    write(xml_out, "gamma", Real(0));

    push(xml_out, "GaugeBC");
    if (spec.masked)
    {
      write(xml_out, "Name", std::string("TEMPORAL_ZONE_GAUGEBC"));
      write(xml_out, "t_dir", spec.params.t_dir);
      write(xml_out, "zero_intervals", spec.params.zero_intervals);
    }
    else
    {
      write(xml_out, "Name", std::string("PERIODIC_GAUGEBC"));
    }
    pop(xml_out);

    pop(xml_out);

    push(xml_out, "NamedObject");
    write(xml_out, "monomial_id", spec.monomial_id);
    pop(xml_out);
    pop(xml_out);
  }

  void loadNamedMonomials(const std::vector<MonomialSpec>& specs)
  {
    TheNamedObjMap::Instance().erase_all();

    XMLBufferWriter xml_buf;
    push(xml_buf, "Monomials");
    for (size_t i = 0; i < specs.size(); ++i)
      appendGaugeMonomial(xml_buf, specs[i]);
    pop(xml_buf);

    std::istringstream xml_stream(xml_buf.str());
    XMLReader xml_in(xml_stream);
    readNamedMonomialArray(xml_in, "/Monomials");
  }

  multi1d<std::string> buildIdList(const std::vector<MonomialSpec>& specs)
  {
    multi1d<std::string> ids(specs.size());
    for (int i = 0; i < ids.size(); ++i)
      ids[i] = specs[i].monomial_id;
    return ids;
  }

  LatticeBoolean buildSiteMask(const TemporalZoneGaugeBCParams& params)
  {
    LatticeInteger t = Layout::latticeCoordinate(params.t_dir);
    LatticeBoolean site_mask = false;

    for (int i = 0; i < params.zero_intervals.size(); ++i)
    {
      const TemporalZoneInterval& interval = params.zero_intervals[i];
      site_mask |= (t >= interval.t_start) && (t <= interval.t_end);
    }

    return site_mask;
  }

  void verifyMaskedMomentum(const std::string& label,
                            const multi1d<LatticeColorMatrix>& p,
                            const TemporalZoneGaugeBCParams& params)
  {
    LatticeBoolean site_mask = buildSiteMask(params);
    LatticeColorMatrix zero_matrix = zero;
    Double active_norm = zero;

    for (int mu = 0; mu < Nd; ++mu)
    {
      Double masked_norm = norm2(where(site_mask, p[mu], zero_matrix));
      if (toBool(masked_norm > tolerance))
        fail(label + ": masked momentum remained nonzero");

      active_norm += norm2(where(site_mask, zero_matrix, p[mu]));
    }

    check(toBool(active_norm > tolerance),
          label + ": expected nonzero momentum in the active region");
  }

  void testDiscovery()
  {
    {
      std::vector<MonomialSpec> specs;
      specs.push_back(makePeriodicSpec("GAUGE_MONOMIAL", "periodic"));
      loadNamedMonomials(specs);

      Handle<GaugeMonomial> discovered =
        discoverMomentumMaskingGaugeMonomial(buildIdList(specs));
      check(discovered.operator->() == 0,
            "periodic-only Hamiltonian should not discover a mask source");
    }

    {
      std::vector<MonomialSpec> specs;
      specs.push_back(makePeriodicSpec("GAUGE_MONOMIAL", "periodic"));
      specs.push_back(makeMaskedSpec("GAUGE_MONOMIAL", "masked", 0, 1));
      loadNamedMonomials(specs);

      Handle<GaugeMonomial> discovered =
        discoverMomentumMaskingGaugeMonomial(buildIdList(specs));
      check(discovered.operator->() != 0,
            "nontrivial gauge BC should be autodiscovered");
      check(discovered->hasNontrivialGaugeBC(),
            "discovered gauge monomial should report a nontrivial BC");
    }

    {
      std::vector<MonomialSpec> specs;
      specs.push_back(makeMaskedSpec("GAUGE_MONOMIAL", "masked_a", 0, 1));
      specs.push_back(makeMaskedSpec("CONST_GAUGE_MONOMIAL", "masked_b", 0, 1));
      loadNamedMonomials(specs);

      Handle<GaugeMonomial> discovered =
        discoverMomentumMaskingGaugeMonomial(buildIdList(specs));
      check(discovered.operator->() != 0,
            "compatible gauge and const-gauge monomials should share a mask source");
      check(discovered->hasNontrivialGaugeBC(),
            "compatible discovery should keep the nontrivial mask");
    }

    {
      std::vector<MonomialSpec> specs;
      specs.push_back(makeMaskedSpec("GAUGE_MONOMIAL", "masked_left", 0, 1));
      specs.push_back(makeMaskedSpec("GAUGE_MONOMIAL", "masked_right", 2, 3));
      loadNamedMonomials(specs);

      bool threw = false;
      try
      {
        discoverMomentumMaskingGaugeMonomial(buildIdList(specs));
      }
      catch (const std::string& e)
      {
        threw = true;
        check(e.find("masked_left") != std::string::npos,
              "incompatible-mask error should name the first monomial");
        check(e.find("masked_right") != std::string::npos,
              "incompatible-mask error should name the second monomial");
      }

      check(threw, "incompatible nontrivial masks should raise an error");
    }
  }

  void testMomentumMaskRefresh()
  {
    std::vector<MonomialSpec> specs;
    specs.push_back(makeMaskedSpec("GAUGE_MONOMIAL", "masked", 0, 1));
    loadNamedMonomials(specs);

    Handle<GaugeMonomial> discovered =
      discoverMomentumMaskingGaugeMonomial(buildIdList(specs));
    check(discovered.operator->() != 0,
          "refresh test requires a discovered momentum mask source");

    multi1d<LatticeColorMatrix> p(Nd);
    multi1d<LatticeColorMatrix> q(Nd);
    for (int mu = 0; mu < Nd; ++mu)
    {
      p[mu] = zero;
      q[mu] = Real(1);
    }

    GaugeFieldState hmc_state(p, q);
    RefreshHarness hmc_refresh(discovered);
    hmc_refresh.run(hmc_state);
    verifyMaskedMomentum("LatColMatHMCTrj", hmc_state.getP(), specs[0].params);

    GaugeFieldState const_hmc_state(p, q);
    ConstRefreshHarness const_hmc_refresh(discovered);
    const_hmc_refresh.run(const_hmc_state);
    verifyMaskedMomentum("ConstLatColMatHMCTrj",
                         const_hmc_state.getP(),
                         specs[0].params);
  }
}

int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  const int nrow_arr[] = {2, 2, 2, 4};
  multi1d<int> nrow(Nd);
  nrow = nrow_arr;
  Layout::setLattSize(nrow);
  Layout::create();

  bool ok = true;
  ok &= GaugeMonomialEnv::registerAll();
  ok &= ConstGaugeMonomialEnv::registerAll();
  QDPIO::cout << "Linkage = " << ok << std::endl;

  testDiscovery();
  testMomentumMaskRefresh();

  TheNamedObjMap::Instance().erase_all();

  Chroma::finalize();
  return 0;
}
