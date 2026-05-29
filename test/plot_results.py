import pandas as pd
import matplotlib.pyplot as plt
import os

# 1. Load data from CSV
csv_path = 'test/data/scalability_results.csv'
plot_dir = 'test/output/plots'

print(f"Reading data from {csv_path}...")
df = pd.read_csv(csv_path)

df.columns = df.columns.str.strip()

# 2. Create folder for plots if it doesn't exist
os.makedirs(plot_dir, exist_ok=True)

grids = sorted(df['n'].unique())
ranks = sorted(df['mpi_ranks'].unique())
colors = ['blue', 'orange', 'green', 'red', 'purple']
grid_colors = {g : colors[i % len(colors)] for i, g in enumerate(grids)}
rank_colors = {r : colors[i % len(colors)] for i, r in enumerate(ranks)}


# EXECUTION TIMES
plt.figure(figsize=(10, 6))
# Draw the serial time line (we take the mean for each grid)
df_serial = df.groupby('n')['serial_time_s'].mean().reset_index()
plt.plot(df_serial['n'], df_serial['serial_time_s'], marker='x', linestyle='--', color='black', label='Serial')

# Draw lines for each MPI rank configuration
for rank, group in df.groupby('mpi_ranks'):
    c = rank_colors[rank]
    plt.plot(group['n'], group['parallel_time_s'], marker='o', linestyle='-', color=c, label=f'Jacobi ({rank} Ranks)')
    plt.plot(group['n'], group['schwarz_time_s'], marker='s', linestyle='--', color=c, label=f'Schwarz ({rank} Ranks)')

plt.title('Execution Times (Serial vs Jacobi vs Schwarz)')
plt.xlabel('Grid Size (n)')
plt.ylabel('Time (seconds)')
plt.xscale('log', base=2)
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.savefig(f'{plot_dir}/1_execution_times.png')
plt.close()

# SPEEDUP vs MPI PROCESSES
plt.figure(figsize=(10, 6))

for n, group in df.groupby('n'):
    c = grid_colors[n]
    plt.plot(group['mpi_ranks'], group['speedup'], marker='o', linestyle='-', color=c, label=f'Jacobi Grid {n}x{n}')
    plt.plot(group['mpi_ranks'], group['schwarz_speedup'], marker='s', linestyle='--', color=c, label=f'Schwarz Grid {n}x{n}')

plt.title('Scalability: Speedup vs MPI Processes')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Speedup')
plt.xticks(ranks) # Force x-axis to show only effective ranks (1, 2, 4)
plt.legend()
plt.grid(True)
plt.savefig(f'{plot_dir}/2_speedup.png')
plt.close()

# L2 ERROR vs GRID SIZE
plt.figure(figsize=(10, 6))

df_error = df[df['mpi_ranks'] == 1].sort_values('n')

plt.plot(df_error['n'], df_error['l2_error'], marker='o', linestyle='-', color='blue', label='Jacobi L2 Error')
plt.plot(df_error['n'], df_error['schwarz_error'], marker='x', linestyle='--', color='red', markersize=10, label='Schwarz L2 Error')
plt.yscale('log')
plt.xscale('log', base=2)
plt.title('Analytical L2 Error')
plt.xlabel('Grid Size (n)')
plt.ylabel('L2 Error (Log Scale)')
plt.grid(True, which="both", ls="--")
plt.savefig(f'{plot_dir}/3_convergence_error.png')
plt.close()

print("Done! The 3 plots have been generated in the 'test/plots/' folder.")