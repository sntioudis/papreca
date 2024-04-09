import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Define the function
def func(ra, rd, t):
    return ra / (ra + rd) * (1 - np.exp( -(ra + rd) * t) )

# Define a function to save x and y values to a text file
def save_to_file(t_values, y_values, filename):
    with open(filename, 'w') as file:
        file.write("t\ty\n")
        for t, y in zip(t_values, y_values):
            file.write(f"{t}\t{y}\n")

# Check if command-line arguments are provided
if len(sys.argv) != 4:
    print("Usage: python3 compare_plots.py ra rd filename")
    sys.exit(1)

# Get ra, rd, and filename from command-line arguments
try:
    ra = float(sys.argv[1])
    rd = float(sys.argv[2])
    filename = sys.argv[3]
except ValueError:
    print("Error: ra and rd must be numeric values.")
    sys.exit(1)

# Create a results directory if it does not exist
results_dir = 'results'
os.makedirs(results_dir, exist_ok=True)

# Define parameters for the function
t_values = np.linspace(0, 4, 1000)  # Generating 100 points between 0 and 5 for t

# Calculate function values for the given ra and rd
y_values = func(ra, rd, t_values)

# Save x and y values to a text file for the given ra and rd
output_txtfilename = f"{results_dir}/output_analytic_ra{ra}_rd{rd}.txt"
save_to_file(t_values, y_values , output_txtfilename )

# Read data from the specified file, skipping the first 5 lines
try:
    data = np.loadtxt(filename, skiprows=5)
    t_values_file = data[:, 0]
    y_values_file = data[:, 1]
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    sys.exit(1)
except ValueError:
    print(f"Error: Invalid data format in file '{filename}'.")
    sys.exit(1)

# Plot both functions
plt.plot(t_values, y_values, label='Analytic Solution')
plt.plot(t_values_file, y_values_file, label='KMC (PAPRECA) solution')
plt.title('Adsorption/Desorption Model')
plt.xlabel('Time (sec)')
plt.ylabel('Surface Coverage (-)')
plt.legend()
plt.grid(True)

# Save plot as a JPEG image
output_jpgfilename = f"{results_dir}/comparison_plot_ra{ra}_rd{rd}.jpg"
plt.savefig(output_jpgfilename , dpi=1200)
