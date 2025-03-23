/*
 * cell_tests.hpp
 *
 * Example cell test functions;
 *
 * Created on: 04 Apr 2022
 * Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../src/slide.hpp"

#include <string>
#include <memory>
#include <fstream>
#include <span>

namespace slide::examples {

inline auto GITT_test()
{
  // This function implements a Galvanostatic Intermittent Titration Technique (GITT) test
  // GITT is used to determine diffusion coefficients in battery electrodes
  // It applies a series of current pulses followed by rest periods to observe voltage relaxation
  
  std::string ID = "temp";  // Identifier for the cell
  Clock clk;  // Create a clock object for timing

  // Temperature check - this test requires T_MODEL=0
  // The GITT data was collected at approximately 24.2°C (23.5-25.9°C range)
  if (settings::T_MODEL != 0) {
    std::cerr << "GITT_test example works with T_MODEL=0 but it is not!\n";
    throw 1234;
  }
  slide::DEG_ID deg{};  // Initialize degradation ID object with default values

  // Create a single particle model (SPM) cell with default parameters
  auto c = Cell_SPM("cell_ancillary", deg, 1, 1, 1, 1);
  c.setBlockDegAndTherm(true);  // Disable degradation and thermal calculations
  c.setT(21.0_degC);  // Set cell temperature to 21°C

  // Maximum concentrations for positive and negative electrodes (mol/m³)
  double Cmaxpos{ 51385 };
  double Cmaxneg{ 30555 };
  double cps{}, cns{};  // Variables to track electrode concentrations (unused)

  auto &st = c.getStateObj();  // Get reference to cell state object
  auto cyc = Cycler(&c, "charge");  // Create cycler for charging operations

  auto d = c;  // Create a copy of the cell for discharge tests
  auto dcyc = Cycler(&d, "discharge");  // Create cycler for discharging operations

  // Initialize cell to a known state by fully discharging and then charging
  ThroughputData th{};  // Structure to track cycling throughput data
  cyc.CCCV(1, 2.7, 0.0001, 1, 0, th);  // Discharge to 2.7V at 1C with 0.0001C cutoff
  cyc.rest(100, 1, 0, th);  // Rest for 100 seconds

  dcyc.CCCV(1, 4.2, 0.0001, 1, 0, th);  // Charge to 4.2V at 1C with 0.0001C cutoff
  dcyc.rest(100, 1, 0, th);  // Rest for 100 seconds

  // Define GITT test parameters
  const auto N_repeat{ 20 };     // Number of pulse-rest cycles
  const auto t_pulse = 1 * 3600; // Duration of each current pulse (1 hour in seconds)
  const auto t_rest = 2 * 3600;  // Duration of each rest period (2 hours in seconds)
  const auto dt = 1;             // Time step for simulation (1 second)
  const auto Crate = 0.05;       // C-rate for pulses (0.05C)
  auto current = Crate * c.Cap(); // Calculate current based on C-rate and cell capacity

  // Create output file for recording GITT test results
  std::ofstream out_GITT{ PathVar::results / "GITT_20x0.05C_1h_rest_2h.csv" };
  out_GITT << "Time [s],"
           << "Current [A],"
           << "Terminal voltage [V],"
           << "Current [A],"
           << "Terminal voltage [V]\n";  // Write CSV header

  // Initialize time counter and record initial state
  double t_all{};
  out_GITT << t_all << ',' << c.I() << ',' << c.V() << ','
           << d.I() << ',' << d.V() << '\n';

  // Execute GITT test cycles
  for (int i{}; i < N_repeat; i++) {
    // Apply current pulse (charging for cell c, discharging for cell d)
    c.setCurrent(-current);  // Negative current for charging
    d.setCurrent(current);   // Positive current for discharging

    // Run pulse for 1 hour (3600 seconds)
    for (int j{}; j < 3600; j++) {
      // Record data at each time step
      out_GITT << t_all << ',' << c.I() << ',' << c.V() << ','
               << d.I() << ',' << d.V() << '\n';

      // Advance simulation by 1 second under constant current
      c.timeStep_CC(1, 1);
      d.timeStep_CC(1, 1);

      t_all += 1;  // Update total elapsed time
    }

    // Set cells to rest (zero current)
    c.setCurrent(0);
    d.setCurrent(0);

    // Run rest period for 2 hours (7200 seconds)
    for (int j{}; j < 7200; j++) {
      // Record data at each time step
      out_GITT << t_all << ',' << c.I() << ',' << c.V() << ','
               << d.I() << ',' << d.V() << '\n';

      // Advance simulation by 1 second under no current
      c.timeStep_CC(1, 1);
      d.timeStep_CC(1, 1);

      t_all += 1;  // Update total elapsed time
    }
  }
  out_GITT.close();  // Close output file
  std::cout << "GITT_test Output folder: results" << '\n';  // Print location of results
}

} // namespace slide::examples