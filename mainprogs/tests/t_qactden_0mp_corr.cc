#include "chroma.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace Chroma;

namespace
{
  struct ExpectedSourceCount
  {
    int dt;
    int value;

    ExpectedSourceCount() : dt(-1), value(-1) {}
  };

  void read(XMLReader& xml, const std::string& path, ExpectedSourceCount& expected)
  {
    XMLReader top(xml, path);
    read(top, "dt", expected.dt);
    read(top, "value", expected.value);
  }

  void write(XMLWriter& xml, const std::string& path, const ExpectedSourceCount& expected)
  {
    push(xml, path);
    write(xml, "dt", expected.dt);
    write(xml, "value", expected.value);
    pop(xml);
  }

  struct CheckerInput
  {
    std::string hmc_file;
    std::string summary_file;
    std::string csv_file;
    multi1d<int> nrow;
    int decay_dir;
    std::string mode;
    int support_guard;
    int max_dt;
    int expected_measurements;
    int expected_timeslices;
    multi1d<ExpectedSourceCount> expected_num_sources;

    CheckerInput()
      : decay_dir(-1),
        support_guard(0),
        max_dt(-1),
        expected_measurements(-1),
        expected_timeslices(-1)
    {
    }
  };

  struct MeasurementResult
  {
    unsigned long update_no;
    multi1d<Double> timeslice_operator;
    multi1d<Double> correlator;
    multi1d<int> num_sources;

    MeasurementResult() : update_no(0) {}
  };

  void fail(const std::string& message)
  {
    QDPIO::cerr << "t_qactden_0mp_corr: " << message << std::endl;
    QDP_abort(1);
  }

  void check(bool condition, const std::string& message)
  {
    if (!condition)
      fail(message);
  }

  void readCheckerInput(XMLReader& xml, const std::string& path, CheckerInput& input)
  {
    XMLReader input_xml(xml, path);

    {
      XMLReader file_xml(input_xml, "Input");
      read(file_xml, "file", input.hmc_file);
    }

    if (input_xml.count("Output") != 0)
    {
      XMLReader output_xml(input_xml, "Output");
      if (output_xml.count("summary_file") != 0)
        read(output_xml, "summary_file", input.summary_file);
      if (output_xml.count("csv_file") != 0)
        read(output_xml, "csv_file", input.csv_file);
    }

    {
      XMLReader geometry_xml(input_xml, "Geometry");
      read(geometry_xml, "nrow", input.nrow);
      read(geometry_xml, "decay_dir", input.decay_dir);
    }

    {
      XMLReader reducer_xml(input_xml, "Reducer");
      read(reducer_xml, "mode", input.mode);
      read(reducer_xml, "support_guard", input.support_guard);
      read(reducer_xml, "max_dt", input.max_dt);
    }

    if (input_xml.count("Checks") != 0)
    {
      XMLReader checks_xml(input_xml, "Checks");
      if (checks_xml.count("expected_measurements") != 0)
        read(checks_xml, "expected_measurements", input.expected_measurements);
      if (checks_xml.count("expected_timeslices") != 0)
        read(checks_xml, "expected_timeslices", input.expected_timeslices);
      if (checks_xml.count("expected_num_sources") != 0)
      {
        XMLReader sources_xml(checks_xml, "expected_num_sources");
        input.expected_num_sources.resize(sources_xml.count("elem"));
        for (int i = 0; i < input.expected_num_sources.size(); ++i)
        {
          std::ostringstream xpath;
          xpath << "elem[" << (i + 1) << "]";
          read(sources_xml, xpath.str(), input.expected_num_sources[i]);
        }
      }
    }
  }

  void checkGeometry(const CheckerInput& input)
  {
    check(input.mode == "PERIODIC_ALL_T",
          "only PERIODIC_ALL_T mode is implemented in the first checker");
    check(input.nrow.size() == Nd, "nrow size must match Nd");
    check(input.decay_dir >= 0 && input.decay_dir < Nd,
          "decay_dir must be in [0, Nd)");
    check(input.max_dt >= 0, "max_dt must be nonnegative");
    check(input.max_dt <= input.nrow[input.decay_dir] / 2,
          "max_dt must be <= floor(T/2)");
  }

  void checkHmcGeometry(XMLReader& hmc_xml, const CheckerInput& input)
  {
    multi1d<int> hmc_nrow;
    read(hmc_xml, "/hmc/Input/Params/HMCTrj/nrow", hmc_nrow);
    check(hmc_nrow.size() == input.nrow.size(), "HMC nrow size mismatch");

    for (int mu = 0; mu < input.nrow.size(); ++mu)
    {
      if (hmc_nrow[mu] != input.nrow[mu])
      {
        std::ostringstream os;
        os << "HMC nrow mismatch at direction " << mu
           << ": expected " << input.nrow[mu]
           << ", got " << hmc_nrow[mu];
        fail(os.str());
      }
    }
  }

  int latticeVolume(const multi1d<int>& nrow)
  {
    int volume = 1;
    for (int mu = 0; mu < nrow.size(); ++mu)
      volume *= nrow[mu];
    return volume;
  }

  MeasurementResult reducePeriodicQactden(unsigned long update_no,
                                          XMLReader& qactden_xml,
                                          const CheckerInput& input)
  {
    const int length = input.nrow[input.decay_dir];
    const int volume = latticeVolume(input.nrow);

    MeasurementResult result;
    result.update_no = update_no;
    result.timeslice_operator.resize(length);
    result.timeslice_operator = zero;
    result.correlator.resize(input.max_dt + 1);
    result.num_sources.resize(input.max_dt + 1);

    XMLReader olattice_xml(qactden_xml, "naiveTopCharge/OLattice");
    const int n_sites = olattice_xml.count("elem");
    check(n_sites == volume, "unexpected OLattice site count in naiveTopCharge");

    for (int site = 0; site < n_sites; ++site)
    {
      std::ostringstream xpath;
      xpath << "elem[" << (site + 1) << "]";

      double q = 0.0;
      read(olattice_xml, xpath.str(), q);

      const multi1d<int> coord = crtesn(site, input.nrow);
      result.timeslice_operator[coord[input.decay_dir]] += Double(q);
    }

    for (int dt = 0; dt <= input.max_dt; ++dt)
    {
      double accum = 0.0;
      for (int t0 = 0; t0 < length; ++t0)
      {
        const int t1 = (t0 + dt) % length;
        accum += toDouble(result.timeslice_operator[t0]) *
                 toDouble(result.timeslice_operator[t1]);
      }

      result.correlator[dt] = Double(accum / double(length));
      result.num_sources[dt] = length;
    }

    return result;
  }

  std::vector<MeasurementResult> collectMeasurements(XMLReader& hmc_xml,
                                                     const CheckerInput& input)
  {
    std::vector<MeasurementResult> results;

    XMLReader updates_xml(hmc_xml, "/hmc/doHMC/MCUpdates");
    const int n_updates = updates_xml.count("./elem");

    for (int i = 1; i <= n_updates; ++i)
    {
      std::ostringstream update_xpath;
      update_xpath << "./elem[" << i << "]/Update";
      XMLReader update_xml(updates_xml, update_xpath.str());

      if (update_xml.count("InlineObservables") == 0)
        continue;

      unsigned long update_no = 0;
      read(update_xml, "update_no", update_no);

      XMLReader observables_xml(update_xml, "InlineObservables");
      const int n_observables = observables_xml.count("./elem");

      for (int obs = 1; obs <= n_observables; ++obs)
      {
        std::ostringstream observable_xpath;
        observable_xpath << "./elem[" << obs << "]";
        XMLReader observable_xml(observables_xml, observable_xpath.str());

        if (observable_xml.count("QActDen") == 0)
          continue;

        XMLReader qactden_xml(observable_xml, "QActDen");

        unsigned long measurement_update_no = 0;
        read(qactden_xml, "update_no", measurement_update_no);
        check(measurement_update_no == update_no,
              "QActDen update number disagrees with enclosing HMC update");

        results.push_back(reducePeriodicQactden(measurement_update_no,
                                                qactden_xml,
                                                input));
      }
    }

    return results;
  }

  multi1d<Double> computeMeanCorrelator(const std::vector<MeasurementResult>& results,
                                        int max_dt)
  {
    multi1d<Double> mean(max_dt + 1);
    mean = zero;

    if (results.empty())
      return mean;

    for (std::size_t i = 0; i < results.size(); ++i)
    {
      for (int dt = 0; dt <= max_dt; ++dt)
        mean[dt] += results[i].correlator[dt];
    }

    for (int dt = 0; dt <= max_dt; ++dt)
      mean[dt] /= Double(results.size());

    return mean;
  }

  void runChecks(const CheckerInput& input, const std::vector<MeasurementResult>& results)
  {
    if (input.expected_measurements >= 0)
    {
      check(int(results.size()) == input.expected_measurements,
            "unexpected number of QActDen measurements");
    }

    check(!results.empty(), "no QActDen measurements found in the HMC XML");

    for (std::size_t i = 0; i < results.size(); ++i)
    {
      const MeasurementResult& result = results[i];

      if (input.expected_timeslices >= 0)
      {
        check(result.timeslice_operator.size() == input.expected_timeslices,
              "unexpected timeslice operator length");
      }

      check(result.correlator.size() == input.max_dt + 1,
            "unexpected correlator length");
      check(result.num_sources.size() == input.max_dt + 1,
            "unexpected num_sources length");

      for (int j = 0; j < input.expected_num_sources.size(); ++j)
      {
        const int dt = input.expected_num_sources[j].dt;
        check(dt >= 0 && dt < result.num_sources.size(),
              "expected_num_sources requested an out-of-range dt");
        check(result.num_sources[dt] == input.expected_num_sources[j].value,
              "num_sources(dt) check failed");
      }
    }
  }

  void writeSummary(const CheckerInput& input,
                    const std::vector<MeasurementResult>& results)
  {
    if (input.summary_file.empty())
      return;

    XMLFileWriter xml_out(input.summary_file);
    push(xml_out, "qactden_0mp_corr_summary");

    push(xml_out, "Input");
    write(xml_out, "hmc_file", input.hmc_file);
    write(xml_out, "nrow", input.nrow);
    write(xml_out, "decay_dir", input.decay_dir);
    write(xml_out, "mode", input.mode);
    write(xml_out, "support_guard", input.support_guard);
    write(xml_out, "max_dt", input.max_dt);
    pop(xml_out);

    write(xml_out, "measurement_count", int(results.size()));
    write(xml_out, "mean_correlator", computeMeanCorrelator(results, input.max_dt));

    push(xml_out, "Measurements");
    for (std::size_t i = 0; i < results.size(); ++i)
    {
      push(xml_out, "elem");
      write(xml_out, "update_no", results[i].update_no);
      write(xml_out, "timeslice_operator", results[i].timeslice_operator);
      write(xml_out, "correlator", results[i].correlator);
      write(xml_out, "num_sources", results[i].num_sources);
      pop(xml_out);
    }
    pop(xml_out);

    pop(xml_out);
    xml_out.close();
  }

  void writeCsv(const CheckerInput& input,
                const std::vector<MeasurementResult>& results)
  {
    if (input.csv_file.empty())
      return;

    std::ofstream out(input.csv_file.c_str());
    check(out.good(), "failed to open CSV output file");

    out << "update_no,series,index,value,num_sources\n";
    out << std::setprecision(17);

    for (std::size_t i = 0; i < results.size(); ++i)
    {
      for (int t = 0; t < results[i].timeslice_operator.size(); ++t)
      {
        out << results[i].update_no
            << ",O_q,"
            << t
            << ","
            << toDouble(results[i].timeslice_operator[t])
            << ",\n";
      }

      for (int dt = 0; dt < results[i].correlator.size(); ++dt)
      {
        out << results[i].update_no
            << ",C_cfg,"
            << dt
            << ","
            << toDouble(results[i].correlator[dt])
            << ","
            << results[i].num_sources[dt]
            << "\n";
      }
    }
  }
}

int main(int argc, char *argv[])
{
  Chroma::initialize(&argc, &argv);

  int status = 0;

  try
  {
    XMLReader xml_in(Chroma::getXMLInputFileName());
    CheckerInput input;
    readCheckerInput(xml_in, "/qactden_0mp_corr_check", input);
    checkGeometry(input);

    XMLReader hmc_xml(input.hmc_file);
    checkHmcGeometry(hmc_xml, input);

    const std::vector<MeasurementResult> results =
      collectMeasurements(hmc_xml, input);

    runChecks(input, results);
    writeSummary(input, results);
    writeCsv(input, results);

    QDPIO::cout << "t_qactden_0mp_corr: passed" << std::endl;
    QDPIO::cout << "  measurements: " << results.size() << std::endl;
    QDPIO::cout << "  hmc_xml: " << input.hmc_file << std::endl;
    if (!input.summary_file.empty())
      QDPIO::cout << "  summary_xml: " << input.summary_file << std::endl;
    if (!input.csv_file.empty())
      QDPIO::cout << "  csv: " << input.csv_file << std::endl;
  }
  catch (const std::string& e)
  {
    QDPIO::cerr << e << std::endl;
    status = 1;
  }
  catch (std::exception& e)
  {
    QDPIO::cerr << "t_qactden_0mp_corr: standard exception: "
                << e.what() << std::endl;
    status = 1;
  }
  catch (...)
  {
    QDPIO::cerr << "t_qactden_0mp_corr: unknown exception" << std::endl;
    status = 1;
  }

  Chroma::finalize();
  return status;
}
