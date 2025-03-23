/*
 * drive_cycles.hpp
 *
 *  Example drive cycle applications for battery simulation
 *  This file contains utilities for testing batteries under different
 *  driving profiles (like the Artemis cycle)
 *
 *  Created on: 02 Nov 2022
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../src/slide.hpp"

#include <string>
#include <memory>
#include <fstream>
#include <span>

namespace slide::examples {
/**
 * @brief Creates and initializes a standard example cell with default parameters
 * 
 * This function creates a Single Particle Model (SPM) cell with predefined
 * degradation settings, including SEI (Solid Electrolyte Interface) growth, 
 * cracking and LAM (Loss of Active Material).
 * 
 * @return A unique pointer to the initialized Cell_SPM object
 */
inline auto get_exampleCell()
{
  // Make one cell with standard parameters// Degradation settings
  slide::DEG_ID deg;
  // Kinetic SEI + porosity + Dai/Laresgoiti LAM

  deg.SEI_id.add_model(4); // add model 4 - Kinetic SEI growth model
  deg.SEI_porosity = 1;    // Enable porosity effects for SEI

  deg.CS_id.add_model(0);  // add cracking/stress model 0
  deg.CS_diffusion = 0;    // Disable diffusion effects for cracking/stress

  deg.LAM_id.add_model(1); // add Loss of Active Material model 1 (Dai/Laresgoiti)
  deg.pl_id = 0;           // No plating model

  auto example_cell = make<Cell_SPM>("cell_ancillary", deg, 1, 1, 1, 1);

  // Specify the OCV (Open Circuit Voltage) parameters 
  // (calculated previously by determineOCV::estimateOCVparameters)
  OCVparam ocvfit;
  ocvfit.elec_surf = (3.4 / 2.65) * 0.0982; // electrode surface Cap_ratio (effective surface area)
  ocvfit.ep = 0.5;                          // volume fraction of active material in the cathode
  ocvfit.en = 0.5;                          // volume fraction of active material in the anode
  ocvfit.thickp = 70e-6;                    // thickness of the cathode [m]
  ocvfit.thickn = 73.5e-6;                  // thickness of the anode [m]
  ocvfit.lifracpini = 0.6862;               // lithium fraction in the cathode at 50% SOC
  ocvfit.lifracnini = 0.4843;               // lithium fraction in the anode at 50% SOC
  ocvfit.cmaxp = 51385;                     // maximum lithium concentration in the cathode [mol m-3]
  ocvfit.cmaxn = 30555;                     // maximum lithium concentration in the anode [mol m-3]
  ocvfit.cap = 3.4;                         // the capacity of the cell [Ah]
  ocvfit.Vmax = 4.2;                        // maximum voltage of the cell [V]
  ocvfit.Vmin = 2.6;                        // minimum voltage of the cell [V]

  // Initialize the cell with the concentration parameters
  example_cell->setInitialConcentration(ocvfit.cmaxp, ocvfit.cmaxn, ocvfit.lifracpini, ocvfit.lifracnini);
  example_cell->setGeometricParameters(ocvfit.cap, ocvfit.elec_surf, ocvfit.ep, ocvfit.en, ocvfit.thickp, ocvfit.thickn);

  example_cell->setT(22.0_degC);    // set the temperature of the cell to 22°C
  example_cell->setTenv(22.0_degC); // set the environmental temperature to 22°C

  // Print the initial state of the cell
  std::cout << "Voltage: " << example_cell->V() << " V.\n";
  std::cout << "Current: " << example_cell->I() << " A.\n";
  std::cout << "Capacity: " << example_cell->Cap() << " Ah.\n";

  return example_cell;
}

/**
 * @brief Initializes a cell to a specific voltage and current state
 * 
 * This function brings the cell to a target voltage and current by performing 
 * charge/discharge operations, while tracking the accumulated Ah throughput and energy.
 * During initialization, degradation and thermal effects are blocked to ensure
 * consistent starting conditions.
 *
 * @param su Reference to the cell to initialize
 * @param I_0 Target current [A]
 * @param V_0 Target voltage [V]
 * @param Ah Reference to variable that will track Ah throughput
 * @param ttot Reference to variable that will track total time [s]
 * @param dWh Reference to variable that will track energy throughput [Wh]
 */
auto inline init_cell_manual(auto &su, double I_0, double V_0, double &Ah, double &ttot, double &dWh)
{
  // Initialize counters
  Ah = 0;
  ttot = 0;
  dWh = 0;

  // Disable degradation and thermal models during initialization
  su->setBlockDegAndTherm(true);

  // Case 1: Very low/no current (voltage control mode)
  if (std::abs(I_0) < 0.01) {
    const double V_lim = V_0;
    const auto testV = su->V();
    
    // Need to discharge if target voltage is below current voltage
    if (V_lim > su->V()) {
      su->setCurrent(-0.02 * su->Cap(), true); // Small negative current (charging)
      while (V_lim > su->V()) {
        su->timeStep_CC(0.1, 1);          // Take small time steps (0.1s)
        Ah += 0.5 * 0.1;                  // Track charge throughput
        ttot += 0.1;                      // Track time elapsed
        dWh += 0.5 * 0.1 * su->V();       // Track energy throughput
      }
    } 
    // Need to charge if target voltage is above current voltage
    else {
      su->setCurrent(0.02 * su->Cap(), true); // Small positive current (discharging)
      while (V_lim < su->V()) {
        su->timeStep_CC(0.1, 1);          // Take small time steps (0.1s)
        Ah += 0.5 * 0.1;                  // Track charge throughput
        ttot += 0.1;                      // Track time elapsed
        dWh += 0.5 * 0.1 * su->V();       // Track energy throughput
      }
    }

    // Set final target current
    su->setCurrent(I_0, true);
  } 
  // Case 2: Significant current (current control mode)
  else {
    // Charging operation (negative current in SLIDE convention)
    if (I_0 < 0) 
    {
      // First discharge cell to make room for charging
      const double V_lim = std::max(su->Vmin(), (V_0 - 0.2)); // Lower voltage limit with safety margin
      
      // Only discharge if not already at low enough voltage
      if (V_lim < su->V()) 
        su->setCurrent(5, true); // High discharge current to empty cell faster

      while (V_lim < su->V()) {
        su->timeStep_CC(2, 1);            // Take larger time steps (2s)
        Ah += 5 * 2;                      // Track charge throughput
        ttot += 2;                        // Track time elapsed
        dWh += 5 * 2 * su->V();           // Track energy throughput
      }

      // Apply target charging current
      su->setCurrent(I_0, true);

      // Charge until target voltage is reached
      while (V_0 > su->V()) {
        su->timeStep_CC(1, 1);
        Ah += std::abs(I_0);
        ttot += 1;
        dWh += std::abs(I_0) * su->V();
      }
    } 
    // Discharging operation (positive current in SLIDE convention)
    else {
      // First charge cell to prepare for discharge
      const double V_lim = std::min(su->Vmax(), (V_0 + 0.2)); // Upper voltage limit with safety margin
      
      // Only charge if not already at high enough voltage
      if (V_lim > su->V()) 
        su->setCurrent(-5, true); // High charge current to fill cell faster

      while (V_lim > su->V()) {
        su->timeStep_CC(2, 1);
        Ah += 5 * 2;
        ttot += 2;
        dWh += 5 * 2 * su->V();
      }

      // Apply target discharge current
      su->setCurrent(I_0, true);

      // Discharge until target voltage is reached
      while (V_0 < su->V()) {
        su->timeStep_CC(1, 1);
        Ah += std::abs(I_0);
        ttot += 1;
        dWh += std::abs(I_0) * su->V();
      }
    }
  }

  // Convert tracking variables to hour-based units
  Ah /= 3600.0;     // Convert Ampere-seconds to Ampere-hours
  dWh /= 3600.0;    // Convert Watt-seconds to Watt-hours
  
  // Re-enable degradation and thermal models for actual simulation
  su->setBlockDegAndTherm(false);
  
  // Reset counters for the main simulation
  Ah = 0;
  ttot = 0;
  dWh = 0;
}

/**
 * @brief Runs a simulation of the Artemis drive cycle
 * 
 * This function:
 * 1. Generates an OCV curve for reference
 * 2. Tests the cell under the Artemis drive cycle profile
 * 3. Tracks voltage, SOC and other parameters during cycling
 * 4. Performs experiments starting from different SOC levels
 * 
 * The Artemis profile is a standardized drive cycle used to evaluate 
 * automotive/battery performance under realistic driving conditions.
 */
inline void drive_cycle_artemis()
{
  // ID for file naming and tracking
  std::string ID = "temp";
  Clock clk;      // Clock for timing the simulation

  // Check that temperature model is disabled as required for this simulation
  // Our data is between 23.5 and 25.9 with mean 24.2 C temperature. 26.44 for test data.
  if (settings::T_MODEL != 0) {
    std::cerr << "drive_cycle_artemis works with T_MODEL=0 but it is not!\n";
    throw 1234;
  }

  // Path to the scaled Artemis drive cycle profile
  std::string profile_path{ "profiles/drive_cycles/ArtemisM_scaled.csv" };

  // Create and configure a cell with standard degradation settings
  slide::DEG_ID deg;
  
  // Configure degradation models:
  deg.SEI_id.add_model(4); // Kinetic SEI growth model
  deg.SEI_porosity = 1;    // Enable porosity effects for SEI

  deg.CS_id.add_model(0);  // Basic cracking/stress model
  deg.CS_diffusion = 0;    // Disable diffusion effects for cracking/stress

  deg.LAM_id.add_model(1); // Dai/Laresgoiti Loss of Active Material model
  deg.pl_id = 0;           // No lithium plating model

  // Create the cell and initialize parameters
  auto c = Cell_SPM("cell_ancillary", deg, 1, 1, 1, 1);
  c.setBlockDegAndTherm(true);  // Block degradation during setup
  c.setT(21.0_degC);            // Set cell temperature to 21°C

  // Maximum concentration values for positive and negative electrodes
  double Cmaxpos{ 51385 };  // Maximum Li concentration in positive electrode [mol/m³]
  double Cmaxneg{ 30555 };  // Maximum Li concentration in negative electrode [mol/m³]
  double cps, cns;          // Variables to store current surface concentrations

  // Get reference to the cell's state object for detailed monitoring
  auto &st = c.getStateObj();
  
  // Create a cycler to handle standardized cycling procedures
  auto cyc = Cycler(&c, "discharge");
  ThroughputData th{};      // To track charge throughput

  // First fully discharge the cell and let it rest to reach equilibrium
  cyc.CCCV(1, 2.7, 0.0001, 1, 0, th);  // Discharge to 2.7V with C/1 rate
  cyc.rest(100, 1, 0, th);             // Rest for 100 seconds with 1s timesteps

  // Create a CSV file to store the lithium concentration fractions during the OCV curve generation
  // This file will track how lithium concentrations change as the cell is charged
  std::ofstream out_fraction{ PathVar::results / "Kokam_NMC_16Ah_OCV.csv" };

  // Write the CSV header with column descriptions:
  // - Current applied to the cell (in Amperes)
  // - Cell voltage (in Volts)
  // - Positive electrode (cathode) lithium fraction (dimensionless) - calculated as surface concentration divided by max concentration
  // - Negative electrode (anode) lithium fraction (dimensionless) - calculated as surface concentration divided by max concentration
  // - Normalized lithium concentration at 3rd node of positive electrode discretization (dimensionless) 
  // - Normalized lithium concentration at 3rd node of negative electrode discretization (dimensionless)
  out_fraction << "Current [A],"
               << "Voltage [V],"
               << "Li+ fraction [-],"
               << "Li- fraction [-],"
               << "zp(3),"
               << "zn(3)\n";

  // Get the initial surface concentrations in both electrodes
  c.getCSurf(cps, cns, false);
  
  // Store initial values to track changes and determine when to write data
  double V_old{ c.V() }, fp_old{ cps / Cmaxpos }, fn_old{ cns / Cmaxneg };
  
  // Write the first data point to the file (initial state)
  out_fraction << c.I() << ',' << c.V() << ',' << fp_old << ',' << fn_old << ',' << st.zp(3) << ',' << st.zn(3) << '\n';

  // Charge the cell slowly (0.02C rate) while recording OCV curve data
  c.setCurrent(-0.02);
  double threshold{ 1e-4 };  // Threshold for recording data points (when changes exceed this value)
  
  // Continue charging until upper voltage limit (4.2V) is reached
  while (c.V() < 4.2) {
    c.getCSurf(cps, cns, false);
    auto fp{ cps / Cmaxpos }, fn{ cns / Cmaxneg };

    // Record data when voltage or concentration changes significantly
    if (std::abs(V_old - c.V()) > threshold || std::abs(fp_old - fp) > threshold || std::abs(fn_old - fn) > threshold) {
      V_old = c.V();
      fp_old = fp;
      fn_old = fn;
      out_fraction << c.I() << ',' << c.V() << ',' << fp << ',' << fn << ',' << st.zp(3) << ',' << st.zn(3) << '\n';
    }
    c.timeStep_CC(1, 1);  // 1-second time step with constant current
  }

  out_fraction.close();  // Close the OCV data file

  // Print final state information
  c.getCSurf(cps, cns, false);
  std::cout << "V: " << c.V() << " cps, cns : " << cps / Cmaxpos << ", " << cns / Cmaxneg << ',' << st.zp(3) << ',' << st.zn(3) << "\n";
  
  // Print all concentration values across discretized electrode nodes
  for (auto z_i : st.z())
    std::cout << z_i << ' ';
  std::cout << '\n';

  // Lookup tables for SOC-voltage-concentration relationships
  // These arrays provide reference points for the relationship between
  // voltage, SOC, and normalized concentrations at different electrodes
  double SOC_vs_Volt[] = { 2.7, 3.4872, 3.5606, 3.6206, 3.6492, 3.6845, 3.7635, 3.8388, 3.9289, 4.0439, 4.2 };
  double SOC_vs_zn3[] = { 0.017456, 0.072543, 0.127650, 0.182757, 0.237865, 0.292971, 0.348078, 0.403186, 0.458292, 0.513400, 0.568507 };
  double SOC_vs_zp3[] = { 0.670246, 0.630914, 0.591568, 0.552221, 0.512874, 0.473528, 0.434182, 0.394835, 0.355489, 0.316142, 0.276795 };

  double SOC_vs_fp[] = { 0.970530, 0.913563, 0.856588, 0.799614, 0.742639, 0.685665, 0.628691, 0.571716, 0.514742, 0.457767, 0.400792 };
  double SOC_vs_fn[] = { 0.028906, 0.120170, 0.211423, 0.302676, 0.393928, 0.485181, 0.576434, 0.667687, 0.758939, 0.850192, 0.941444 };

  // Verify lookup tables by setting concentrations and checking resulting voltages
  for (int i = 0; i < 11; i++) {
    c.setC(SOC_vs_fp[i], SOC_vs_fn[i]);
    std::cout << "i = " << i << " V = " << c.V() << '\n';
  }

  // Load the Artemis drive cycle profile from CSV file
  DynamicMatrix<double> profile;
  loadCSV_Ncol(PathVar::data / profile_path, profile);

  // Define the experiment function that will run drive cycle from different SOC levels
  auto experiment = [&](std::string name) {
    std::vector<double> voltage;             // Store voltage measurements
    std::vector<State_SPM> states;           // Store complete cell states
    std::vector<double> estimated_SOC;       // Store estimated SOC values
    double cap = c.Cap();                    // Cell capacity in Ah
    
    // Loop through each current value in the profile
    for (auto c_rate : profile.data) {
      c.setCurrent(c_rate * cap, true);      // Set current as C-rate * capacity
      voltage.push_back(c.V());              // Record voltage
      states.push_back(c.getStateObj());     // Record full cell state

      // Estimate SOC based on anode concentration, normalized to 0-100%
      estimated_SOC.push_back(100.0 * ((st.zn(3) - SOC_vs_zn3[0]) / (SOC_vs_zn3[10] - SOC_vs_zn3[0])));

      c.timeStep_CC(1, 1);                   // Step simulation forward by 1 second

      // Stop if SOC drops below 10%
      if (st.zn(3) < SOC_vs_zn3[1])          // SOC_vs_zn3[1] corresponds to 10% SOC
        break;
    }

    // Save voltage and estimated SOC data to CSV
    std::string Vname = std::string("drive_voltage_from_") + name + "_percent.csv";
    std::ofstream out(PathVar::results / Vname, std::ios_base::out);

    out << "Current [A],Voltage[V],Estimated SOC[%]\n";

    for (size_t i = 0; i < voltage.size(); i++)
      out << profile.data[i] * cap << ',' << voltage[i] << ',' << estimated_SOC[i] << '\n';

    out.close();

    // Save detailed state data to CSV for deeper analysis
    // e.g.: drive_states_from_20_percent.csv
    std::string Sname = std::string("drive_states_from_") + name + "_percent.csv";
    std::ofstream Sout(PathVar::results / Sname, std::ios_base::out);

    // Write CSV header with state variable descriptions:
    // I[A]      - Current applied to the cell [Amperes]
    // V[V]      - Cell voltage [Volts]
    // T[T]      - Cell temperature [Kelvin or Celsius]
    // delta     - SEI layer thickness [m]
    // LLI       - Loss of Lithium Inventory [%]
    // thickp    - Positive electrode thickness [m]
    // thickn    - Negative electrode thickness [m]
    // ep        - Positive electrode porosity [-]
    // en        - Negative electrode porosity [-]
    // ap        - Positive electrode specific surface area [m²/m³]
    // an        - Negative electrode specific surface area [m²/m³]
    // CS        - Crack surface area [m²]
    // Dp        - Positive electrode diffusion coefficient [m²/s]
    // Dn        - Negative electrode diffusion coefficient [m²/s]
    // delta_pl  - Plated lithium layer thickness [m]
    // zp_0..4   - Normalized Li concentration at nodes 0-4 in positive electrode [-]
    // zn_0..4   - Normalized Li concentration at nodes 0-4 in negative electrode [-]
    // rDCp      - Positive electrode DC resistance [Ohm]
    // rDCn      - Negative electrode DC resistance [Ohm]
    // rDCcc     - Current collector DC resistance [Ohm]
    Sout << "I[A],V[V],T[T],delta,LLI,thickp,thickn,ep,en,ap,an,CS,Dp,Dn,delta_pl,"
         << "zp_0,zp_1,zp_2,zp_3,zp_4,zn_0,zn_1,zn_2,zn_3,zn_4,rDCp,rDCn,rDCcc\n";

    // Write all states to file
    for (auto st_i : states) {
      for (auto s : std::span<double>(st_i.begin(), st_i.begin() + 28))
        Sout << s << ',';

      Sout << '\n';
    }

    Sout.close();
  };

  // Run experiments starting from different SOC levels (90% down to 20%)
  // The following block is commented out but kept for reference
  // c.setCurrent(0.01);
  // while (c.V() > SOC_vs_Volt[9]) // 9 -> 90%
  //   c.timeStep_CC(1, 1);

  // Run experiments from 90% down to 20% SOC in steps of 10%
  for (int i = 9; i > 1; i--) {
    c.setC(SOC_vs_fp[i], SOC_vs_fn[i]);  // Set concentration based on lookup table
    experiment(std::to_string(i) + "0");  // Run experiment and name it by SOC percentage
  }

  // Print final results
  std::cout << "V: " << c.V() << '\n';
  std::cout << "SOC: " << 100 * c.SOC() << '\n';
  std::cout << "Finished drive_cycle example in " << clk << ".\n";
  // Print the output folder location
  std::cout << "drive_cycle_artemis Output folder: results" << '\n';
}

} // namespace slide::examples