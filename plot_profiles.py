"""Plot phonon temperature, spinon (magnetic) temperature and the comparison
T_uncoupled - T_coupled, as in the thesis presentation. Run the solver first.

Usage: python plot_profiles.py [title]
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

ph = np.loadtxt("phtemp.txt")   # x, T_phonon, T_uncoupled - T_phonon, t
mg = np.loadtxt("magtemp.txt")  # x, T_magnetic, t
times = np.unique(ph[:, 3])
shown = [times[np.argmin(abs(times - t))] for t in np.arange(0.1, 1.01, 0.1)]

fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
cmap = plt.cm.viridis
for i, t in enumerate(shown):
    c = cmap(i / max(len(shown) - 1, 1))
    p = ph[np.isclose(ph[:, 3], t)]
    m = mg[np.isclose(mg[:, 2], t)]
    ax[0].plot(p[:, 0], p[:, 1], color=c, label=f"t={t:.2f}")
    ax[1].plot(m[:, 0], m[:, 1], color=c)
    ax[2].plot(p[:, 0], p[:, 2], color=c)
ax[0].set_title("Phonon temperature")
ax[1].set_title("Spinon temperature")
ax[2].set_title("Comparison: T_uncoupled - T_coupled")
ax[2].axhline(0, color="grey", lw=0.8)
for a in ax:
    a.set_xlabel("x")
    a.grid(alpha=0.3)
ax[0].legend(fontsize=7)
if len(sys.argv) > 1:
    fig.suptitle(sys.argv[1])
fig.tight_layout()
fig.savefig("profiles.png", dpi=130)

e = np.loadtxt("energy.txt")  # t, phonon energy, magnetic energy, total
print(f"Total energy drift: {100 * (e[-1, 3] - e[0, 3]) / e[0, 3]:.2f}%")
