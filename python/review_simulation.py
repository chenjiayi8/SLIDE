import pybamm
import matplotlib.pyplot as plt

# Load the simulation results
solution = pybamm.load("simulation_results.h5")

# Extract variables for analysis
time = solution.variables["Time [h]"]
voltage = solution.variables["Terminal voltage [V]"]

# Print basic statistics
print(f"Minimum voltage: {min(voltage)} V")
print(f"Maximum voltage: {max(voltage)} V")

# Plot the voltage over time
plt.plot(time, voltage)
plt.xlabel("Time [h]")
plt.ylabel("Voltage [V]")
plt.title("Battery Voltage over Time")
plt.savefig("voltage_over_time.png")