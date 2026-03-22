import subprocess
import sys
from pathlib import Path


def submit(config_path):
    config = Path(config_path)

    if not config.exists():
        print(f"Config not found: {config}")
        sys.exit(1)

    slurm_script = Path("scripts/rbfe/run.slurm")

    cmd = ["sbatch", str(slurm_script), str(config)]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("Submission failed:")
        print(result.stderr)
        sys.exit(1)

    print(result.stdout)
