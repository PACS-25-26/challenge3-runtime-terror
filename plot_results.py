import pandas as pd
import matplotlib.pyplot as plt
import os

# 1. Load data from CSV
csv_path = 'test/data/scalability_results.csv'
print(f"Reading data from {csv_path}...")
df = pd.read_csv(csv_path)

df.columns = df.columns.str.strip()

# 2. Create folder for plots if it doesn't exist
os.makedirs('plots', exist_ok=True)

plt.figure(figsize=(8, 6))

# Draw the serial time line (we take the mean for each grid)
df_serial = df.groupby('n')['serial_time_s'].mean().reset_index()
plt.plot(df_serial['n'], df_serial['serial_time_s'], marker='x', linestyle='--', color='black', label='Serial')

# Draw lines for each MPI rank configuration
for ranks, group in df.groupby('mpi_ranks'):
    plt.plot(group['n'], group['parallel_time_s'], marker='o', label=f'Parallel ({ranks} MPI Ranks)')

plt.title('Execution Times (Serial vs Parallel Hybrid)')
plt.xlabel('Grid Size (n)')
plt.ylabel('Time (seconds)')
plt.legend()
plt.grid(True)
plt.savefig('plots/1_execution_times.png')
plt.close()

plt.figure(figsize=(8, 6))

for n, group in df.groupby('n'):
    plt.plot(group['mpi_ranks'], group['speedup'], marker='s', label=f'Grid {n}x{n}')

plt.title('Scalability: Speedup vs MPI Processes')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Speedup')
plt.xticks([1, 2, 4]) # Force x-axis to show only 1, 2, 4
plt.legend()
plt.grid(True)
plt.savefig('plots/2_speedup.png')
plt.close()

plt.figure(figsize=(8, 6))

df_error = df[df['mpi_ranks'] == 1].sort_values('n')

plt.plot(df_error['n'], df_error['l2_error'], marker='^', color='red')
plt.yscale('log')
plt.xscale('log')
plt.title('Analytical L2 Error')
plt.xlabel('Grid Size (n)')
plt.ylabel('L2 Error (Log Scale)')
plt.grid(True, which="both", ls="--")
plt.savefig('plots/3_convergence_error.png')
plt.close()

print("Done! The 3 plots have been generated in the 'plots/' folder.")