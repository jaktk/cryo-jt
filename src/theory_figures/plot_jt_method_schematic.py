import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')

BLUE = '#4A90E2'
RED = '#C8102E'
K_CURV = 2.2


def T_iso(p):
    """Schematic concave isenthalp T(p) in arbitrary units."""
    return 0.08 + 0.80 * (1 - np.exp(-K_CURV * p)) / (1 - np.exp(-K_CURV))


def dT_iso(p):
    return 0.80 * K_CURV * np.exp(-K_CURV * p) / (1 - np.exp(-K_CURV))


def _axes(ax):
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.0)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.minorticks_off()
    ax.plot(1.05, 0, '>k', ms=4, clip_on=False, transform=ax.get_yaxis_transform())
    ax.plot(0, 1.0, '^k', ms=4, clip_on=False, transform=ax.get_xaxis_transform())
    ax.text(1.06, -0.02, r'$p$', ha='left', va='top', transform=ax.transData)
    ax.text(-0.02, 1.0, r'$T$', ha='right', va='center')
    p = np.linspace(0, 1.0, 200)
    ax.plot(p, T_iso(p), color=BLUE, lw=1.4)
    ax.text(0.97, T_iso(0.97) + 0.03, r'$h = \mathrm{const}$', color=BLUE,
            ha='right', va='bottom', rotation=8)


def _point(ax, p, label, color='k', lw=0.8):
    T = T_iso(p)
    ax.plot([p, p], [0, T], ':', color=color, lw=lw)
    ax.plot([0, p], [T, T], ':', color=color, lw=lw)
    ax.text(p, -0.03, label[0], color=color, ha='center', va='top')
    ax.text(-0.02, T, label[1], color=color, ha='right', va='center')
    ax.plot(p, T, 'o', color=color, mec='k', mew=0.5, ms=4, zorder=5)


def _angle(ax, x0, y0, slope, r, color, label):
    deg = np.degrees(np.arctan(slope / ax_aspect(ax)))
    ax.add_patch(Arc((x0, y0), 2 * r, 2 * r * ax_aspect(ax), theta1=0,
                     theta2=deg, color=color, lw=1.0))
    mid = np.radians(deg / 2)
    f = 0.6 if r > 0.15 else 1.6  # label inside large arcs, outside small ones
    ax.text(x0 + f * r * np.cos(mid), y0 + f * r * ax_aspect(ax) * np.sin(mid),
            label, color=color, ha='center', va='center')


def ax_aspect(ax):
    """Data-to-display ratio so arcs look circular on the page."""
    w, h = ax.get_window_extent().width, ax.get_window_extent().height
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    return (h / (y1 - y0)) ** -1 * (w / (x1 - x0))


def _tangent(ax, p_r, p_lo, p_hi):
    T_r, s = T_iso(p_r), dT_iso(p_r)
    p = np.array([p_lo, p_hi])
    ax.plot(p, T_r + s * (p - p_r), color=RED, lw=1.1)
    return T_r + s * (p_lo - p_r), s


def plot_small():
    fig, ax = plt.subplots(figsize=(3.0, 2.45))
    fig.canvas.draw()
    _axes(ax)
    p2, p1 = 0.27, 0.70
    secant = (T_iso(p1) - T_iso(p2)) / (p1 - p2)
    # reference point where the tangent is parallel to the secant
    p_r = -np.log(secant * (1 - np.exp(-K_CURV)) / (0.80 * K_CURV)) / K_CURV
    _point(ax, p2, (r'$p_2$', r'$T_2$'))
    _point(ax, p1, (r'$p_1$', r'$T_1$'))
    _point(ax, p_r, (r'$p_{\mathrm{r}}$', r'$T_{\mathrm{r}}$'), color=RED)
    ax.plot([p2, p1], [T_iso(p2), T_iso(p1)], 'k-', lw=1.1)
    x0 = p2 - 0.13
    y0, s = _tangent(ax, p_r, x0, p1 - 0.04)
    ax.plot([x0, x0 + 0.24], [y0, y0], color=RED, lw=1.1)
    _angle(ax, x0, y0, s, 0.09, RED, r'$\alpha$')
    ax.plot([p2, p2 + 0.30], [T_iso(p2)] * 2, 'k-', lw=0.6)
    _angle(ax, p2, T_iso(p2), secant, 0.24, 'k', r'$\beta$')
    fig.savefig(os.path.join(OUT_DIR, 'JT_theory_small_dp.pdf'))


def plot_large():
    fig, ax = plt.subplots(figsize=(3.0, 2.45))
    fig.canvas.draw()
    _axes(ax)
    for p, lab in [(0.22, '4'), (0.38, '3'), (0.62, '2'), (0.86, '1')]:
        _point(ax, p, (rf'$p_{lab}$', rf'$T_{lab}$'))
    p_r = 0.45
    _point(ax, p_r, (r'$p_{\mathrm{r}}$', r'$T_{\mathrm{r}}$'), color=RED)
    x0 = 0.26
    y0, s = _tangent(ax, p_r, x0, 0.66)
    ax.plot([x0, x0 + 0.26], [y0, y0], color=RED, lw=1.1)
    _angle(ax, x0, y0, s, 0.09, RED, r'$\alpha$')
    fig.savefig(os.path.join(OUT_DIR, 'JT_theory_large_dp.pdf'))


def main():
    plt.style.use(STYLE)
    os.makedirs(OUT_DIR, exist_ok=True)
    plot_small()
    plot_large()


if __name__ == '__main__':
    main()
