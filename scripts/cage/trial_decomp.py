import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

# Configuration
num_runs = 10
# mu_values = [0.0, 0.001, 0.00316227766, 0.01, 0.0316227766, 0.1, 0.31622776601, 1]
mu_values = sys.argv[1]
base_output_dir = "../runs/cage"
command_base = "./bin/iter_decomp"
max_parallel = 4  # Number of parallel threads per mu

# Fixed command arguments (excluding --mu and --output)
base_args = [
    "--bumpy",
    "--nv", "8",
    "--np", "16",
    "--points", "100000",
    "--dt", "0.001",
    "--dphi", "0.01",
    "--steps", "6",
    "--angles", "1000",
]

def run_command(index, mu, output_dir):
    output_file = os.path.join(output_dir, f"{index:06d}.h5")
    cmd = [command_base, "--output", output_file] + base_args + ["--mu", str(mu)]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd)

# Iterate over mu values sequentially
for mu in mu_values:
    output_dir = os.path.join(base_output_dir, f"mu_{mu}")
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n=== Running mu = {mu} ===")

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = [
            executor.submit(run_command, i + 1, mu, output_dir)
            for i in range(num_runs)
        ]
        for future in as_completed(futures):
            future.result()  # Propagate exceptions

print("All simulations completed.")
