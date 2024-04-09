import os
import re
import sys
import matplotlib.pyplot as plt

# Function to read and plot data from log files
def plot_log_files(folder_path, file_numbers):
    # Get list of log files with names in the format "distribution XXX.log"
    log_files = sorted([file for file in os.listdir(folder_path) if re.match(r'^distribution \d+\.log$', file)])

    for number in file_numbers:
        filename = f"distribution {number}.log"
        if filename in log_files:
            file_path = os.path.join(folder_path, filename)
            x_vals, y_vals = [], []

            # Read data from log file
            with open(file_path, 'r') as file:
                lines = file.readlines()[6:]  # Skip first 6 lines
                for line in lines:
                    data = line.split()
                    x_vals.append(float(data[0]))
                    y_vals.append(float(data[1]))

            # Plot data
            plt.plot(x_vals, y_vals, label=filename)
        else:
            print(f"File 'distribution {number}.log' not found. Skipping...")

    # Add labels and legend
    plt.xlabel('Height (LAMMPS units) of x-y slice')
    plt.ylabel('#Atoms')
    plt.title('Brownian Motion Model')
    plt.legend()

    # Save plot as JPEG with DPI=1200
    plt.savefig('distributions.jpg', dpi=1200)

    # Show plot
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) >= 3:
        folder_path = sys.argv[1]
        file_numbers = [int(number) for number in sys.argv[2:]]
    else:
        folder_path = input("Enter the path to the results folder: ")
        file_numbers = [int(number) for number in input("Enter the file numbers (space-separated): ").split()]
        
    plot_log_files(folder_path, file_numbers)
