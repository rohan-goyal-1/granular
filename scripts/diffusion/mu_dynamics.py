import subprocess
import os
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

num_runs = 1
base_output_dir = "runs/dynamics"
command_base = "./bin/dynamics"
max_parallel_runs = 4

base_args = [
    "--nv", "1",
    "--np", "16",
    "--dphi", "0.001",
    "--steps", "3",
    "--dt", "0.001",
    "--dt_min", "1.0",
    "--dt_max", "1000000.0",
]

def run_simulation(run_index, mu_value, output_dir):
    """
    Constructs and executes a single simulation command.

    Args:
        run_index (int): The index of the run (e.g., 1, 2, 3...).
        mu_value (float): The friction coefficient.
        output_dir (str): The directory to save the output file in.
    """
    # Format the output filename with leading zeros (e.g., 000001.h5).
    output_file = os.path.join(output_dir, f"{run_index:06d}.h5")

    # Combine all parts of the command.
    command = [
        command_base,
        "--output", output_file,
        "--mu", str(mu_value)
    ] + base_args

    print(f"Starting Run {run_index} (mu={mu_value}): {' '.join(command)}")

    subprocess.run(command)

    print(f"Finished Run {run_index} (mu={mu_value})")
    return f"Run {run_index} completed successfully."

def main():
    """
    Main function to parse arguments and manage simulation runs.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Run multiple dynamics simulations for a given friction value (mu).")
    parser.add_argument(
        "--mu",
        type=float,
        required=True,
        help="The friction coefficient (mu) for the simulation."
    )
    args = parser.parse_args()
    mu_value = args.mu

    # --- Directory Setup ---
    output_dir = os.path.join(base_output_dir, f"mu_{mu_value}")
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n=== Running {num_runs} simulations for mu = {mu_value} ===")
    print(f"Output will be saved in: {output_dir}")
    print(f"Running up to {max_parallel_runs} simulations in parallel.")

    # --- Execute Runs in Parallel ---
    with ThreadPoolExecutor(max_workers=max_parallel_runs) as executor:
        # Submit all simulation tasks to the thread pool.
        future_to_run = {
            executor.submit(run_simulation, i + 1, mu_value, output_dir): i + 1
            for i in range(num_runs)
        }

        # Process results as they complete.
        for future in as_completed(future_to_run):
            run_num = future_to_run[future]
            try:
                result = future.result()
            except subprocess.CalledProcessError as e:
                print(f"Run {run_num} failed with exit code {e.returncode}")
                print(f"  --> Stderr: {e.stderr.strip()}")
                print(f"  --> Stdout: {e.stdout.strip()}")
            except Exception as exc:
                print(f"Run {run_num} generated an exception: {exc}")

    print("\nAll simulations completed.")

if __name__ == "__main__":
    main()
