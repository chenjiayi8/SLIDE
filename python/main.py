import pybamm

# Define the battery model
model = pybamm.lithium_ion.DFN()

# Set up the simulation
parameter_values = pybamm.ParameterValues("Chen2020")
sim = pybamm.Simulation(model, parameter_values=parameter_values)

# Run the simulation
sim.solve([0, 3600])  # Simulate for 1 hour

# Extract variables from the solution
solution = sim.solution
time = solution["Time [h]"].entries
voltage = solution["Voltage [V]"].entries

# Create a plot using Matplotlib and save it to a file
import matplotlib.pyplot as plt

plt.figure()
plt.plot(time, voltage)
plt.xlabel("Time [h]")
plt.ylabel("Voltage [V]")
plt.title("Battery Voltage over Time")
plt.savefig("simulation_plot.png")

print("Simulation completed and plot saved.")