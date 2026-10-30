from pathlib import Path

# repository root, two levels above this file (src/data_analysis/); works in
# a git clone and in the Zenodo archive alike
REPO_ROOT = str(Path(__file__).resolve().parents[2])
