#!/usr/bin/env bash
#SBATCH --job-name=check-orb-env
#SBATCH --partition=volta-gpu
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=check_orb_env.%j.out
#SBATCH --error=check_orb_env.%j.err
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1

set -euo pipefail

module load anaconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate orb

python - <<'PY'
import sys
try:
    import orb_models
    import nvalchemiops
    import warp
except Exception as exc:
    print(f"ERROR: import failed: {exc!r}")
    raise SystemExit(1)
print("OK: imports succeeded")
print("Python:", sys.executable)
print("orb_models:", getattr(orb_models, "__version__", "unknown"))
print("nvalchemiops:", getattr(nvalchemiops, "__version__", "unknown"))
print("warp:", getattr(warp, "__version__", "unknown"))
PY
