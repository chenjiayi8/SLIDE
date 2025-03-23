import os
import socket
import datetime
import pybamm
import numpy as np
import concurrent.futures
import pandas as pd
import matplotlib.pyplot as plt

def run_simulation(config):
    """Run battery simulation with given configuration"""
    # Initialize model based on config
    model_class = pybamm.lithium_ion.DFN if config["model"] == "DFN" else pybamm.lithium_ion.SPM
    model = model_class({
        "SEI": config["sei_type"],
        "lithium plating": config["plating"],
        "loss of active material": config["lam"]
    })

    # Load drive cycle from CSV skipping header
    drive_cycle = pd.read_csv("/app/drive_cycle.csv", header=None, skiprows=1).values
    
    # Create current interpolant - fix the dimension mismatch
    time_points = drive_cycle[:, 0]  # Assuming first column is time
    current_values = drive_cycle[:, 1]  # Assuming second column is current

    current_interpolant = pybamm.Interpolant(
        time_points,
        current_values,
        pybamm.t
    )

    params = pybamm.ParameterValues("Chen2020")
    params.update({
        "SEI kinetic rate constant [m.s-1]": config["sei_k"],
        "Lithium plating kinetic rate constant [m.s-1]": config["plating_k"],
        "Positive electrode LAM constant [s-1]": config["lam_p"],
        "Negative electrode LAM constant [s-1]": config["lam_n"],
        "Positive electrode LAM constant chemical term [s-1]": 2e-20,
        "Positive electrode OCP entropic change [V.K-1]": config["dUdT"],
        "SEI porosity cut-off": config["sei_porosity"],
        "Lithium plating porosity cut-off": config["plating_porosity"],
        "LAM porosity cut-off": config["lam_porosity"],
        "Use OCV aging": config["OCV_aging"],
        "Current function [A]": current_interpolant,
        "Exchange-current density for stripping [A.m-2]": config["stripping_j0"]
    }, check_already_exists=False)

    sim = pybamm.Simulation(model, parameter_values=params)
    solution = sim.solve([0, 3600])
    
    # Create output directory matching C++ structure
    output_dir = f"/app/python/output_{config['id']:04d}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save results to CSV files instead of HDF5
    # Save main simulation data
    time_data = solution["Time [s]"].entries
    voltage_data = solution["Voltage [V]"].entries
    temperature_data = solution["Temperature [K]"].entries
    
    # Core data CSV
    main_df = pd.DataFrame({
        'time': time_data,
        'voltage': voltage_data,
        'current': drive_cycle[:, 0],
        'temperature': temperature_data
    })
    main_df.to_csv(f"{output_dir}/simulation_results.csv", index=False)
    
    # Degradation data CSV
    degradation_data = {
        'time': time_data,
        'sei_thickness': solution["SEI thickness [m]"].entries
    }
    
    if "Lithium plating thickness [m]" in solution:
        degradation_data['li_plating'] = solution["Lithium plating thickness [m]"].entries
    
    if "Positive electrode LAM [m]" in solution:
        degradation_data['lam_positive'] = solution["Positive electrode LAM [m]"].entries
    
    if "Negative electrode LAM [m]" in solution:
        degradation_data['lam_negative'] = solution["Negative electrode LAM [m]"].entries
    
    degradation_df = pd.DataFrame(degradation_data)
    degradation_df.to_csv(f"{output_dir}/degradation_results.csv", index=False)
    
    # Save metadata to a separate CSV
    metadata = {
        'parameter': ['thread_count', 'hostname', 'timestamp'] + list(config.keys()),
        'value': [
            config.get("threads", os.cpu_count()),
            socket.gethostname(),
            datetime.datetime.now().isoformat()
        ] + [str(config[key]) for key in config.keys()]
    }
    
    metadata_df = pd.DataFrame(metadata)
    metadata_df.to_csv(f"{output_dir}/metadata.csv", index=False)
    
    return solution

if __name__ == "__main__":
    # Configuration parameters matching C++ benchmark cases
    configurations = [
        {
            "id": 1,
            "dUdT": 0.0001,
            "model": "DFN",
            "sei_type": "reaction limited",
            "plating": "reversible",
            "lam": "none",
            "sei_k": 1e-12,
            "plating_k": 1e-7,
            "lam_p": 0,
            "lam_n": 0,
            "sei_porosity": 0.05,
            "plating_porosity": 0.1,
            "lam_porosity": 0.2,
            "OCV_aging": False,
            "stripping_j0": 1.0
        },
        {
            "id": 2,
            "dUdT": 0.0002,
            "model": "DFN",
            "sei_type": "solvent-diffusion limited",
            "plating": "irreversible",
            "lam": "cathode",
            "sei_k": 5e-13,
            "plating_k": 5e-8,
            "lam_p": 1e-6,
            "lam_n": 0,
            "sei_porosity": 0.1,
            "plating_porosity": 0.15,
            "lam_porosity": 0.25,
            "OCV_aging": True,
            "stripping_j0": 1.0
        },
        {
            "id": 3,
            "dUdT": 0.0003,
            "model": "SPM",
            "sei_type": "electron-migration limited",
            "plating": "reversible",
            "lam": "anode",
            "sei_k": 1e-13,
            "plating_k": 1e-8,
            "lam_p": 0,
            "lam_n": 2e-6,
            "sei_porosity": 0.15,
            "plating_porosity": 0.2,
            "lam_porosity": 0.3,
            "OCV_aging": False,
            "stripping_j0": 1.0
        }
    ]

    # Run simulations in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(run_simulation, configurations))

    # Generate consolidated plot
    plt.figure(figsize=(12, 6))
    for i, solution in enumerate(results):
        time = solution["Time [h]"].entries
        voltage = solution["Voltage [V]"].entries
        plt.plot(time, voltage, label=f'Config {i+1}')
    
    plt.xlabel("Time [h]")
    plt.ylabel("Voltage [V]")
    plt.title("Comparative Voltage Profiles")
    plt.legend()
    plt.grid(True)
    plt.savefig("/app/python/simulation_plot.png")
    
    # Calculate and print degradation summary
    print("\nDegradation Summary:")
    for i, solution in enumerate(results):
        config = configurations[i]
        sei_avg = np.mean(solution["SEI thickness [m]"].entries)
        plating_avg = np.mean(solution["Lithium plating thickness [m]"].entries) if "Lithium plating thickness [m]" in solution else 0.0
        lam_p_avg = np.mean(solution["Positive electrode LAM [m]"].entries) if "Positive electrode LAM [m]" in solution else 0.0
        lam_n_avg = np.mean(solution["Negative electrode LAM [m]"].entries) if "Negative electrode LAM [m]" in solution else 0.0
        
        print(f"Config {config['id']}:")
        print(f"  SEI Thickness Avg: {sei_avg:.2e} m")
        print(f"  Li Plating Avg: {plating_avg:.2e} m")
        print(f"  Pos LAM Avg: {lam_p_avg:.2e} m")
        print(f"  Neg LAM Avg: {lam_n_avg:.2e} m")

    print("All simulations completed. Results saved to CSV files and plot generated.")