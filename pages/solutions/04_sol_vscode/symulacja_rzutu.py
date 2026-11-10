from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# Dane indeksu: 36 to dwie pierwsze cyfry, a 70 to dwie ostatnie.
V0 = 36.0
ALPHA_DEG = 70.0
G = 9.81
MASS = 1.0
DRAG_COEFFICIENT = 0.05

ALPHA = np.deg2rad(ALPHA_DEG)
VX0 = V0 * np.cos(ALPHA)
VY0 = V0 * np.sin(ALPHA)
OUTPUT_DIR = Path(__file__).parent / "files"


def vacuum_motion(time):
    """Return position and velocity for a projectile without air resistance."""
    x = VX0 * time
    y = VY0 * time - 0.5 * G * time**2
    vx = np.full_like(time, VX0)
    vy = VY0 - G * time
    return x, y, vx, vy


def drag_equations(_time, state):
    x, y, vx, vy = state
    speed = np.hypot(vx, vy)
    ax = -(DRAG_COEFFICIENT / MASS) * speed * vx
    ay = -G - (DRAG_COEFFICIENT / MASS) * speed * vy
    return vx, vy, ax, ay


def ground_contact(_time, state):
    return state[1]


ground_contact.terminal = True
ground_contact.direction = -1


def save_figure(figure, filename):
    figure.tight_layout()
    figure.savefig(OUTPUT_DIR / filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    flight_time_vacuum = 2 * VY0 / G
    range_vacuum = V0**2 * np.sin(2 * ALPHA) / G
    height_vacuum = VY0**2 / (2 * G)
    time_vacuum = np.linspace(0, flight_time_vacuum, 1200)
    x_vacuum, y_vacuum, vx_vacuum, vy_vacuum = vacuum_motion(time_vacuum)

    solution = solve_ivp(
        drag_equations,
        (0, 100),
        (0.0, 0.0, VX0, VY0),
        events=ground_contact,
        dense_output=True,
        rtol=1e-9,
        atol=1e-11,
        max_step=0.02,
    )
    if not solution.success or len(solution.t_events[0]) == 0:
        raise RuntimeError("Symulacja z oporem nie zakończyła się uderzeniem w ziemię.")

    flight_time_drag = solution.t_events[0][0]
    time_drag = np.linspace(0, flight_time_drag, 1200)
    x_drag, y_drag, vx_drag, vy_drag = solution.sol(time_drag)
    time_apex_drag = brentq(lambda time: solution.sol(time)[3], 0, flight_time_drag)
    height_drag = solution.sol(time_apex_drag)[1]
    range_drag = solution.sol(flight_time_drag)[0]

    figure, axis = plt.subplots(figsize=(9, 5.5))
    axis.plot(x_vacuum, y_vacuum, label="Bez oporu", linewidth=2.2)
    axis.plot(x_drag, y_drag, label="Z oporem kwadratowym", linewidth=2.2)
    axis.scatter([0, range_vacuum], [0, 0], color="C0", s=28, zorder=3)
    axis.scatter([range_drag], [0], color="C1", s=28, zorder=3)
    axis.set(title="Tor lotu", xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]")
    axis.grid(True, alpha=0.3)
    axis.legend()
    save_figure(figure, "tor_lotu.png")

    figure, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(time_vacuum, x_vacuum, label="Bez oporu")
    axes[0].plot(time_drag, x_drag, label="Z oporem")
    axes[0].set_ylabel("Położenie x [m]")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].plot(time_vacuum, y_vacuum, label="Bez oporu")
    axes[1].plot(time_drag, y_drag, label="Z oporem")
    axes[1].set(xlabel="Czas t [s]", ylabel="Wysokość y [m]")
    axes[1].grid(True, alpha=0.3)
    save_figure(figure, "polozenie_w_czasie.png")

    figure, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(time_vacuum, vx_vacuum, label="Bez oporu")
    axes[0].plot(time_drag, vx_drag, label="Z oporem")
    axes[0].set_ylabel("Prędkość vx [m/s]")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].plot(time_vacuum, vy_vacuum, label="Bez oporu")
    axes[1].plot(time_drag, vy_drag, label="Z oporem")
    axes[1].set(xlabel="Czas t [s]", ylabel="Prędkość vy [m/s]")
    axes[1].grid(True, alpha=0.3)
    save_figure(figure, "predkosci_w_czasie.png")

    speed_vacuum = np.hypot(vx_vacuum, vy_vacuum)
    speed_drag = np.hypot(vx_drag, vy_drag)
    figure, axis = plt.subplots(figsize=(9, 5.5))
    axis.plot(time_vacuum, speed_vacuum, label="Bez oporu", linewidth=2.2)
    axis.plot(time_drag, speed_drag, label="Z oporem", linewidth=2.2)
    axis.set(title="Wartość prędkości", xlabel="Czas t [s]", ylabel="Szybkość v [m/s]")
    axis.grid(True, alpha=0.3)
    axis.legend()
    save_figure(figure, "szybkosc_w_czasie.png")

    print(f"v0 = {V0:.2f} m/s, alpha = {ALPHA_DEG:.2f} deg ({ALPHA:.6f} rad)")
    print(f"g = {G:.2f} m/s^2, m = {MASS:.2f} kg, k = {DRAG_COEFFICIENT:.3f} kg/m")
    print(f"Zasięg bez oporu: {range_vacuum:.3f} m")
    print(f"Zasięg z oporem:  {range_drag:.3f} m")
    print(f"Wysokość bez oporu: {height_vacuum:.3f} m")
    print(f"Wysokość z oporem:  {height_drag:.3f} m")
    print(f"Czas lotu bez oporu: {flight_time_vacuum:.3f} s")
    print(f"Czas lotu z oporem:  {flight_time_drag:.3f} s")
    print(f"Wykresy zapisano w: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()