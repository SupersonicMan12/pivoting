# Plotting scaffold generated with AI assistance.
# Numerical results were computed locally from gaussian_growth.csv.

## Credits: ChatGPT-5.6 Sol Medium Thinking
## Prompt: Given the format of the gaussian_growth.csv
## visualize the average growth and compare the
## theoretical n^{2/3} for partial and n^{1/2} for complete
 
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "gaussian_growth.csv"
FIGURE_DIR = ROOT.parent / "figures"

FIGURE_DIR.mkdir(exist_ok=True)

EXPECTED_COUNTS = {
    16: 2000,
    32: 2000,
    64: 1500,
    128: 1000,
    256: 500,
    512: 200,
    1024: 40,
}

data = pd.read_csv(DATA_PATH)

assert len(data) == sum(EXPECTED_COUNTS.values())
assert data.groupby("n").size().to_dict() == EXPECTED_COUNTS
assert data["matrix_seed"].is_unique
assert data["partial_growth"].notna().all()
assert data["complete_growth"].notna().all()
assert (data["partial_growth"] >= 1).all()
assert (data["complete_growth"] >= 1).all()

medians = data.groupby("n")[
    ["partial_growth", "complete_growth"]
].median()

n = medians.index.to_numpy(dtype=float)
partial = medians["partial_growth"].to_numpy()
complete = medians["complete_growth"].to_numpy()


def reference_power(values, exponent):
    log_coefficient = np.mean(
        np.log(values) - exponent * np.log(n)
    )
    coefficient = np.exp(log_coefficient)
    return coefficient * n**exponent


partial_reference = reference_power(partial, 2 / 3)
complete_reference = reference_power(complete, 1 / 2)

fig, ax = plt.subplots(figsize=(7, 5))

ax.loglog(
    n,
    partial,
    marker="o",
    label="Partial pivoting median",
)

ax.loglog(
    n,
    complete,
    marker="o",
    label="Complete pivoting median",
)

ax.loglog(
    n,
    partial_reference,
    linestyle="--",
    label=r"Reference slope $n^{2/3}$",
)

ax.loglog(
    n,
    complete_reference,
    linestyle="--",
    label=r"Reference slope $n^{1/2}$",
)

ax.set_xlabel("Matrix dimension $n$")
ax.set_ylabel("Median growth factor")
ax.set_title("Median Gaussian-elimination growth")
ax.grid(True, which="both", linestyle=":", alpha=0.5)
ax.legend()

fig.tight_layout()
fig.savefig(FIGURE_DIR / "median_growth.png", dpi=200)
fig.savefig(FIGURE_DIR / "median_growth.pdf")
plt.close(fig)

print(medians)
print("Saved median-growth figures")