import sys
import subprocess

if sys.argv[1] == "rbfe":
    subprocess.run([
        "sbatch",
        "./scripts/run_simulations.sh",
        sys.argv[2],
    ])