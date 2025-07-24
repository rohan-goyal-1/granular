import subprocess
import threading

# Configuration
start_index = 1
num_runs = 1
output_dir = "runs/cage2/bumpy2"
command_base = "./bin/iter_decomp"
output_template = f"{output_dir}/{{:06d}}.h5"

# Fixed command arguments
args = [
    "--bumpy",
    "--nv", "8",
    "--np", "20",
    "--points", "100000",
    "--dt", "0.001",
    "--dphi", "0.01",
    "--steps", "6",
    "--angles", "1000"
]

def run_command(index):
    output_file = output_template.format(index)
    cmd = [command_base, "--output", output_file] + args
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd)

# Create and start threads
threads = []
for i in range(start_index, start_index + num_runs):
    t = threading.Thread(target=run_command, args=(i,))
    t.start()
    threads.append(t)

# Wait for all threads to finish
for t in threads:
    t.join()

print("All subprocesses finished.")

