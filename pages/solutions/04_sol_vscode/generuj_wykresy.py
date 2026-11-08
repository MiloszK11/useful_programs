"""Generate the figures used by raport.md.

Run from any directory with: python generuj_wykresy.py
Requires NumPy and Matplotlib.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parent / "wykresy"
OUT.mkdir(parents=True, exist_ok=True)

G = 9.81
V0 = 36.0
ALPHA_DEG = 70.0
MASS = 0.145
RHO = 1.225
CD = 0.47
RADIUS = 0.0366
AREA = np.pi * RADIUS**2
DRAG_FACTOR = RHO * CD * AREA / (2 * MASS)

BLUE = "#2463A6"
ORANGE = "#D47722"
GREEN = "#168578"
PURPLE = "#8455A4"
INK = "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def bez_oporu(czas, alpha_rad):
    """Analytic coordinates and velocities for a projectile without drag."""
    vx0 = V0 * np.cos(alpha_rad)
    vy0 = V0 * np.sin(alpha_rad)
    x = vx0 * czas
    y = vy0 * czas - 0.5 * G * czas**2
    vx = np.full_like(czas, vx0)
    vy = vy0 - G * czas
    return x, y, vx, vy


def pochodna(stan):
    """State derivative for quadratic air resistance."""
    x, y, vx, vy = stan
    speed = np.hypot(vx, vy)
    ax = -DRAG_FACTOR * speed * vx
    ay = -G - DRAG_FACTOR * speed * vy
    return np.array([vx, vy, ax, ay])


def krok_rk4(stan, dt):
    """One fourth-order Runge-Kutta step."""
    k1 = pochodna(stan)
    k2 = pochodna(stan + 0.5 * dt * k1)
    k3 = pochodna(stan + 0.5 * dt * k2)
    k4 = pochodna(stan + dt * k3)
    return stan + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def z_oporem(alpha_deg, dt=0.002, t_max=30.0):
    """Integrate until the projectile returns to the launch height."""
    alpha = np.deg2rad(alpha_deg)
    stan = np.array([
        0.0,
        0.0,
        V0 * np.cos(alpha),
        V0 * np.sin(alpha),
    ])
    czasy = [0.0]
    stany = [stan.copy()]
    czas = 0.0

    while czas < t_max:
        nastepny = krok_rk4(stan, dt)
        if czas > 0 and nastepny[1] <= 0:
            udzial = stan[1] / (stan[1] - nastepny[1])
            czas += udzial * dt
            stan = stan + udzial * (nastepny - stan)
            stan[1] = 0.0
            czasy.append(czas)
            stany.append(stan.copy())
            break
        czas += dt
        stan = nastepny
        czasy.append(czas)
        stany.append(stan.copy())
    else:
        raise RuntimeError("The projectile did not return to launch height.")

    return np.asarray(czasy), np.asarray(stany)


def zapisz(fig, nazwa):
    fig.savefig(OUT / f"{nazwa}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    alpha = np.deg2rad(ALPHA_DEG)
    vx0 = V0 * np.cos(alpha)
    vy0 = V0 * np.sin(alpha)
    czas_lotu_idealny = 2 * vy0 / G
    zasieg_idealny = vx0 * czas_lotu_idealny
    wysokosc_idealna = vy0**2 / (2 * G)

    czas_idealny = np.linspace(0, czas_lotu_idealny, 800)
    x0, y0, vx0_t, vy0_t = bez_oporu(czas_idealny, alpha)
    czas_opor, stan_opor = z_oporem(ALPHA_DEG)
    x_opor, y_opor, vx_opor, vy_opor = stan_opor.T
    predkosc_opor = np.hypot(vx_opor, vy_opor)

    print(f"v0 = {V0:.2f} m/s; alpha = {ALPHA_DEG:.2f} deg")
    print(f"vx0 = {vx0:.4f} m/s; vy0 = {vy0:.4f} m/s")
    print(
        "No drag: "
        f"flight time = {czas_lotu_idealny:.4f} s, "
        f"range = {zasieg_idealny:.4f} m, "
        f"maximum height = {wysokosc_idealna:.4f} m"
    )
    print(
        "Quadratic drag: "
        f"flight time = {czas_opor[-1]:.4f} s, "
        f"range = {x_opor[-1]:.4f} m, "
        f"maximum height = {np.max(y_opor):.4f} m"
    )
    print(
        f"Air model: A = {AREA:.6f} m^2, "
        f"rho*Cd*A/(2m) = {DRAG_FACTOR:.6f} 1/m"
    )

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    ax.plot(x0, y0, color=BLUE, label="Bez oporu")
    ax.plot(x_opor, y_opor, color=ORANGE, label="Opór kwadratowy")
    ax.scatter([0, zasieg_idealny, x_opor[-1]], [0, 0, 0],
               color=[INK, BLUE, ORANGE], zorder=3)
    ax.set(
        xlabel="Położenie poziome x [m]",
        ylabel="Wysokość y [m]",
        title=f"Tor lotu: v₀ = {V0:g} m/s, α = {ALPHA_DEG:g}°",
    )
    ax.legend()
    zapisz(fig, "tory_lotu")

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                             layout="constrained")
    axes[0].plot(czas_idealny, x0, color=BLUE, label="x bez oporu")
    axes[0].plot(czas_opor, x_opor, color=ORANGE, label="x z oporem")
    axes[0].plot(czas_idealny, y0, color=GREEN, label="y bez oporu")
    axes[0].plot(czas_opor, y_opor, color=PURPLE, label="y z oporem")
    axes[0].set_ylabel("Położenie [m]")
    axes[0].set_title("Położenie w funkcji czasu")
    axes[0].legend(ncol=2, fontsize=9)
    axes[1].plot(czas_idealny, vx0_t, color=BLUE, label="vₓ bez oporu")
    axes[1].plot(czas_opor, vx_opor, color=ORANGE, label="vₓ z oporem")
    axes[1].plot(czas_idealny, vy0_t, color=GREEN, label="vᵧ bez oporu")
    axes[1].plot(czas_opor, vy_opor, color=PURPLE, label="vᵧ z oporem")
    axes[1].set(xlabel="Czas t [s]", ylabel="Prędkość [m/s]")
    axes[1].set_title("Składowe prędkości")
    axes[1].legend(ncol=2, fontsize=9)
    zapisz(fig, "polozenie_predkosc")

    energia_idealna = 0.5 * MASS * (vx0_t**2 + vy0_t**2) + MASS * G * y0
    energia_opor = (
        0.5 * MASS * predkosc_opor**2 + MASS * G * y_opor
    )
    moc_tracona = (
        0.5 * RHO * CD * AREA * predkosc_opor**3
    )
    praca_oporu = np.concatenate((
        [0.0],
        np.cumsum(
            0.5 * (moc_tracona[1:] + moc_tracona[:-1])
            * np.diff(czas_opor)
        ),
    ))
    energia_poczatkowa = 0.5 * MASS * V0**2
    blad_bilansu = np.max(np.abs(energia_opor + praca_oporu
                                 - energia_poczatkowa))

    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    ax.plot(czas_idealny, energia_idealna, color=BLUE,
            label="Energia mechaniczna bez oporu")
    ax.plot(czas_opor, energia_opor, color=ORANGE,
            label="Energia mechaniczna z oporem")
    ax.plot(czas_opor, praca_oporu, color=GREEN, label="Energia rozproszona")
    ax.plot(czas_opor, energia_opor + praca_oporu, "--", color=PURPLE,
            label="Energia mechaniczna + rozproszona")
    ax.set(
        xlabel="Czas t [s]",
        ylabel="Energia [J]",
        title="Zmiana energii wskutek oporu powietrza",
    )
    ax.legend(fontsize=8, ncol=2)
    zapisz(fig, "energia")

    katy = np.arange(1.0, 90.0, 1.0)
    zasiegi_opor = np.array([
        z_oporem(kat, dt=0.02)[1][-1, 0] for kat in katy
    ])
    zasiegi_idealne = V0**2 * np.sin(2 * np.deg2rad(katy)) / G
    kat_opt = katy[np.argmax(zasiegi_opor)]
    zasieg_opt = np.max(zasiegi_opor)

    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    ax.plot(katy, zasiegi_idealne, color=BLUE, label="Bez oporu")
    ax.plot(katy, zasiegi_opor, color=ORANGE, label="Opór kwadratowy")
    ax.scatter([ALPHA_DEG], [x_opor[-1]], color=PURPLE, zorder=3,
               label=f"α = {ALPHA_DEG:g}°")
    ax.axvline(kat_opt, color=INK, linestyle=":", linewidth=1)
    ax.set(
        xlabel="Kąt wyrzutu α [°]",
        ylabel="Zasięg [m]",
        title="Zasięg w zależności od kąta wyrzutu",
        xlim=(0, 90),
    )
    ax.legend()
    zapisz(fig, "zasieg_kat")

    print(
        f"Maximum range with drag in 1-degree scan: "
        f"alpha = {kat_opt:.1f} deg, range = {zasieg_opt:.4f} m"
    )
    print(f"Maximum energy-balance residual: {blad_bilansu:.4e} J")
    if not np.all(np.isfinite(stan_opor)) or blad_bilansu > 0.05:
        raise RuntimeError("Simulation validation failed.")
    print(f"Saved 4 figures to {OUT}")


if __name__ == "__main__":
    main()
