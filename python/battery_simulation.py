import pybamm
import pandas as pd
import matplotlib.pyplot as plt

# Load drive cycle data
drive_cycle_path = "drive_cycle.csv"
data = pd.read_csv(drive_cycle_path)
time = data['time'].values
current = data['current'].values

# Define the interpolant for current
current_interpolant = pybamm.Interpolant(time, current, pybamm.t)

# Load the default parameter values
parameter_values = pybamm.ParameterValues("Ai2020")

# Create model
model = pybamm.lithium_ion.DFN()

# Replace the input current with the interpolant
parameter_values.update({"Current function [A]": current_interpolant})

# Set up simulation
sim = pybamm.Simulation(model, parameter_values=parameter_values)

# Run simulation with t_eval set to the time points in the drive cycle
t_eval = time
sim.solve(t_eval=t_eval)

# Save results to a file
sim.save("simulation_results.h5")

# Extract variables for plotting
time = sim.solution.t
voltage = sim.solution["Terminal voltage [V]"].data

# Plot using Matplotlib and save the plot to a file
plt.plot(time, voltage)
plt.xlabel("Time [s]")
plt.ylabel("Voltage [V]")
plt.title("Battery Voltage over Time")
plt.savefig("simulation_plot.png")