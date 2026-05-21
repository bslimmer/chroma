#include "chroma.h"
#include "actions/gauge/gaugebcs/temporal_zone_gaugebc.h"

#include <sstream>

using namespace Chroma;

namespace
{
  const Double tolerance = Double(1.0e-12);

  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_temporal_zone_gaugebc: " << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  std::string dirLabel(int mu)
  {
    std::ostringstream os;
    os << mu;
    return os.str();
  }

  int entrySlice(const TemporalZoneInterval& interval, int t_extent)
  {
    return (interval.t_start + t_extent - 1) % t_extent;
  }

  LatticeBoolean buildDirMask(const TemporalZoneGaugeBCParams& params, int mu)
  {
    LatticeInteger t = Layout::latticeCoordinate(params.t_dir);
    LatticeBoolean mask = false;
    const int t_extent = Layout::lattSize()[params.t_dir];

    for (int i = 0; i < params.zero_intervals.size(); ++i)
    {
      const TemporalZoneInterval& interval = params.zero_intervals[i];
      mask |= (t >= interval.t_start) && (t <= interval.t_end);

      if (mu == params.t_dir)
        mask |= (t == entrySlice(interval, t_extent));
    }

    return mask;
  }

  void verifyZeroing(const TemporalZoneGaugeBCParams& params,
                     const std::string& label,
                     XMLWriter& xml_out)
  {
    TemporalZoneGaugeBC bc(params);
    check(bc.nontrivialP(), label + ": expected nontrivialP() to be true");

    multi1d<LatticeColorMatrix> u(Nd);
    multi1d<LatticeColorMatrix> u_before(Nd);
    multi1d<LatticeColorMatrix> ds_u(Nd);
    multi1d<LatticeColorMatrix> ds_u_before(Nd);

    for (int mu = 0; mu < Nd; ++mu)
    {
      gaussian(u[mu]);
      gaussian(ds_u[mu]);
      u_before[mu] = u[mu];
      ds_u_before[mu] = ds_u[mu];
    }

    bc.modify(u);
    bc.zero(ds_u);

    LatticeColorMatrix z = zero;

    push(xml_out, label);
    write(xml_out, "t_dir", params.t_dir);
    write(xml_out, "zero_intervals", params.zero_intervals);

    for (int mu = 0; mu < Nd; ++mu)
    {
      LatticeBoolean dir_mask = buildDirMask(params, mu);
      Double modify_diff = norm2(u[mu] - u_before[mu]);
      Double masked_norm = norm2(where(dir_mask, ds_u[mu], z));
      Double unmasked_diff = norm2(where(dir_mask, z, ds_u[mu] - ds_u_before[mu]));

      push(xml_out, "Direction");
      write(xml_out, "mu", mu);
      write(xml_out, "modify_diff", modify_diff);
      write(xml_out, "masked_norm", masked_norm);
      write(xml_out, "unmasked_diff", unmasked_diff);
      pop(xml_out);

      if (toBool(modify_diff > tolerance))
        fail(label + ": modify() changed gauge links in direction " + dirLabel(mu));

      if (toBool(masked_norm > tolerance))
        fail(label + ": zero() left nonzero values in masked region for direction " + dirLabel(mu));

      if (toBool(unmasked_diff > tolerance))
        fail(label + ": zero() changed unmasked values in direction " + dirLabel(mu));
    }

    pop(xml_out);
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

  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_temporal_zone_gaugebc");

  TemporalZoneGaugeBCParams params;
  params.zero_intervals.resize(1);
  params.zero_intervals[0].t_start = 0;
  params.zero_intervals[0].t_end = 1;

  XMLBufferWriter xml_buf;
  push(xml_buf, "TemporalZoneGaugeBCTest");
  write(xml_buf, "Params", params);
  pop(xml_buf);

  std::istringstream xml_stream(xml_buf.str());
  XMLReader xml_in(xml_stream);

  TemporalZoneGaugeBCParams parsed_params;
  read(xml_in, "/TemporalZoneGaugeBCTest/Params", parsed_params);

  check(parsed_params.t_dir == params.t_dir, "serialized t_dir did not round-trip");
  check(parsed_params.zero_intervals.size() == params.zero_intervals.size(),
        "serialized zero_intervals size did not round-trip");
  check(parsed_params.zero_intervals[0].t_start == params.zero_intervals[0].t_start,
        "serialized t_start did not round-trip");
  check(parsed_params.zero_intervals[0].t_end == params.zero_intervals[0].t_end,
        "serialized t_end did not round-trip");

  XMLBufferWriter default_buf;
  push(default_buf, "TemporalZoneGaugeBCTest");
  push(default_buf, "Params");
  write(default_buf, "zero_intervals", params.zero_intervals);
  pop(default_buf);
  pop(default_buf);

  std::istringstream default_stream(default_buf.str());
  XMLReader default_in(default_stream);
  TemporalZoneGaugeBCParams default_params(default_in, "/TemporalZoneGaugeBCTest/Params");
  check(default_params.t_dir == Nd - 1, "default t_dir did not fall back to Nd - 1");

  TemporalZoneGaugeBCParams non_wrap_params;
  non_wrap_params.zero_intervals.resize(1);
  non_wrap_params.zero_intervals[0].t_start = 1;
  non_wrap_params.zero_intervals[0].t_end = 2;
  verifyZeroing(non_wrap_params, "SingleIntervalWithEntryLink", xml_out);

  TemporalZoneGaugeBCParams overlap_params;
  overlap_params.zero_intervals.resize(2);
  overlap_params.zero_intervals[0].t_start = 1;
  overlap_params.zero_intervals[0].t_end = 2;
  overlap_params.zero_intervals[1].t_start = 2;
  overlap_params.zero_intervals[1].t_end = 3;
  verifyZeroing(overlap_params, "OverlappingIntervals", xml_out);

  verifyZeroing(parsed_params, "WraparoundEntryLink", xml_out);

  pop(xml_out);
  xml_out.close();

  Chroma::finalize();
  return 0;
}
