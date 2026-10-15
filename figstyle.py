"""Shared figure style; importing this module applies it."""
from pathlib import Path

import matplotlib.pyplot as plt

FIGURES = Path("figures")

plt.rcParams.update({
    # figures are drawn at the report's text width (6.5 in), so these sizes are the printed sizes
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5, "lines.linewidth": 1.2,
    "lines.markersize": 4, "axes.linewidth": 0.6,
    "legend.frameon": False, "figure.constrained_layout.use": True, "savefig.dpi": 200,
})


def save(fig, name):
    """PNG for the web page, PDF (vector) for the report."""
    (FIGURES / "pdf").mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{name}.png")
    fig.savefig(FIGURES / "pdf" / f"{name}.pdf")
    plt.close(fig)
