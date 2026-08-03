#!/usr/bin/env python3
"""Exact dimensionless steady-uniform solution for two power-law layers.

This is the dimensionless counterpart of
``two_layer_powerlaw_steady_uniform_dimensional.py``. It uses the report's
lower-normal-flow scaling

    H = H_l,
    U = Q_l^0/H_l,
    L = H/S_0,
    T = L/U,
    tau_ref = rho_l S_0 U^2.

The independent dimensionless input set is

    (Fr_l, h_r, r_rho, n_l, n_u, kappa_K),

where

    h_r       = H_u/H_l,
    r_rho     = rho_u/rho_l,
    Fr_l      = U/sqrt(g cos(theta) H_l),
    kappa_K   = (K_u/K_l) (U/H_l)^(n_u-n_l)
              = Lambda_u/Lambda_l.

The lower normal-flow depth and mean velocity are both one, so q_l^0=1.
The code evaluates all remaining tractions, velocities, discharges, and
momentum moments analytically and writes dimensionless profile files.

Only NumPy is required for the calculation. Matplotlib is optional.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import math

import numpy as np


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

@dataclass(frozen=True)
class DimensionlessParameters:
    """Independent dimensionless inputs of the report."""

    froude_lower: float
    depth_ratio: float
    density_ratio: float
    lower_index: float
    upper_index: float
    scaled_consistency_ratio: float


# These defaults are exactly the dimensionless form of the default physical
# case in the dimensional companion program.
PARAMETERS = DimensionlessParameters(
    froude_lower=0.13155093503115464,
    depth_ratio=0.8,
    density_ratio=0.875,
    lower_index=0.60,
    upper_index=0.80,
    scaled_consistency_ratio=0.6634223762441973,
)

POINTS_PER_LAYER = 501
MAKE_PLOTS = True
OUTPUT_DIRECTORY = (
    Path(__file__).resolve().parent / "dimensionless_steady_uniform_output"
)


# =============================================================================
# 2. DATA STRUCTURE
# =============================================================================

@dataclass(frozen=True)
class DimensionlessSolution:
    lambda_0: float
    m_lower: float
    m_upper: float
    C_lower: float
    A_lower: float
    B_lower: float
    A_upper: float
    B_upper: float
    Lambda_lower: float
    Lambda_upper: float
    tau_b: float
    tau_interface: float
    interface_velocity: float
    upper_velocity_increment: float
    surface_velocity: float
    lower_mean_velocity: float
    upper_mean_velocity: float
    lower_discharge: float
    upper_discharge: float
    lower_momentum: float
    upper_momentum: float


# =============================================================================
# 3. ANALYTICAL SHAPE FUNCTIONS
# =============================================================================

def validate_parameters(p: DimensionlessParameters) -> None:
    positive_values = {
        "froude_lower": p.froude_lower,
        "depth_ratio": p.depth_ratio,
        "density_ratio": p.density_ratio,
        "lower_index": p.lower_index,
        "upper_index": p.upper_index,
        "scaled_consistency_ratio": p.scaled_consistency_ratio,
    }
    invalid = [
        name
        for name, value in positive_values.items()
        if (not math.isfinite(value)) or value <= 0.0
    ]
    if invalid:
        raise ValueError(
            "These parameters must be finite and positive: " + ", ".join(invalid)
        )
    if not (0.0 < p.lower_index <= 1.0):
        raise ValueError("lower_index must satisfy 0 < n_l <= 1.")
    if not (0.0 < p.upper_index <= 1.0):
        raise ValueError("upper_index must satisfy 0 < n_u <= 1.")
    if POINTS_PER_LAYER < 11:
        raise ValueError("POINTS_PER_LAYER must be at least 11.")


def lower_shape_coefficients(index: float, lam: float) -> Tuple[float, float, float]:
    """Return the analytical lower-profile coefficients C_l, A_l, and B_l."""

    if not (0.0 <= lam < 1.0):
        raise ValueError("The traction ratio lambda must satisfy 0 <= lambda < 1.")

    m = (index + 1.0) / index
    if lam <= 1.0e-14:
        c = 1.0 / m
        a = m / (m + 1.0)
        b = 1.0 - 2.0 / (m + 1.0) + 1.0 / (2.0 * m + 1.0)
        return c, a, b

    eps = 1.0 - lam
    if eps < 1.0e-4:
        c = (
            1.0
            - 0.5 * (m - 1.0) * eps
            + (m - 1.0) * (m - 2.0) * eps**2 / 6.0
            - (m - 1.0) * (m - 2.0) * (m - 3.0) * eps**3 / 24.0
            + (m - 1.0)
            * (m - 2.0)
            * (m - 3.0)
            * (m - 4.0)
            * eps**4
            / 120.0
        )
        a = (
            0.5
            + (m - 1.0) * eps / 12.0
            + (m - 1.0) * eps**2 / 24.0
            - (m - 1.0) * (m**2 - 19.0) * eps**3 / 720.0
            - (m - 3.0) * (m - 1.0) * (m + 3.0) * eps**4 / 480.0
        )
        b = (
            1.0 / 3.0
            + (m - 1.0) * eps / 12.0
            + (m - 1.0) * (m + 7.0) * eps**2 / 180.0
            - (m - 1.0) * (m**2 - 4.0 * m - 17.0) * eps**3 / 720.0
            - (m - 1.0)
            * (m + 2.0)
            * (m**2 + 8.0 * m - 41.0)
            * eps**4
            / 5040.0
        )
        return c, a, b

    log_lam = math.log(lam)
    one_minus_lam_m = -math.expm1(m * log_lam)
    j_m = -math.expm1((m + 1.0) * log_lam) / ((m + 1.0) * eps)
    j_2m = -math.expm1((2.0 * m + 1.0) * log_lam) / (
        (2.0 * m + 1.0) * eps
    )
    c = one_minus_lam_m / (m * eps)
    a = (1.0 - j_m) / one_minus_lam_m
    b = (1.0 - 2.0 * j_m + j_2m) / one_minus_lam_m**2
    return c, a, b


def lower_shape(s: np.ndarray, index: float, lam: float) -> np.ndarray:
    s = np.asarray(s, dtype=float)
    if np.any((s < 0.0) | (s > 1.0)):
        raise ValueError("The lower local coordinate s must lie in [0, 1].")
    m = (index + 1.0) / index
    if lam <= 1.0e-14:
        return 1.0 - (1.0 - s) ** m
    eps = 1.0 - lam
    numerator = -np.expm1(m * np.log1p(-eps * s))
    denominator = -math.expm1(m * math.log(lam))
    return numerator / denominator


def upper_shape(r: np.ndarray, index: float) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if np.any((r < 0.0) | (r > 1.0)):
        raise ValueError("The upper local coordinate r must lie in [0, 1].")
    m = (index + 1.0) / index
    return 1.0 - (1.0 - r) ** m


# =============================================================================
# 4. EXACT DIMENSIONLESS NORMAL FLOW
# =============================================================================

def solve_dimensionless(p: DimensionlessParameters) -> DimensionlessSolution:
    """Evaluate the exact dimensionless steady-uniform state."""

    validate_parameters(p)

    h_r = p.depth_ratio
    r_rho = p.density_ratio
    Fr = p.froude_lower

    lam = r_rho * h_r / (1.0 + r_rho * h_r)
    m_l = (p.lower_index + 1.0) / p.lower_index
    m_u = (p.upper_index + 1.0) / p.upper_index
    C_l, A_l, B_l = lower_shape_coefficients(p.lower_index, lam)

    A_u = m_u / (m_u + 1.0)
    B_u = 1.0 - 2.0 / (m_u + 1.0) + 1.0 / (2.0 * m_u + 1.0)

    # Gravity-based normal-flow tractions in tau_ref=rho_l S_0 U_l^2 units.
    tau_b = (1.0 + r_rho * h_r) / Fr**2
    tau_I = r_rho * h_r / Fr**2

    # Lower compatibility condition after h_l=1 and mean u_l=1 are imposed.
    Lambda_l = tau_b * (A_l * C_l) ** p.lower_index
    Lambda_u = p.scaled_consistency_ratio * Lambda_l

    U_I = 1.0 / A_l
    W_u = h_r / m_u * (tau_I / Lambda_u) ** (1.0 / p.upper_index)
    Ubar_l = 1.0
    Ubar_u = U_I + A_u * W_u
    q_l = 1.0
    q_u = h_r * Ubar_u
    M_l = U_I**2 * B_l
    M_u = h_r * (U_I**2 + 2.0 * A_u * U_I * W_u + B_u * W_u**2)

    calculated_values = {
        "lambda_0": lam,
        "Lambda_lower": Lambda_l,
        "Lambda_upper": Lambda_u,
        "tau_b": tau_b,
        "tau_interface": tau_I,
        "interface_velocity": U_I,
        "upper_velocity_increment": W_u,
        "upper_mean_velocity": Ubar_u,
        "upper_discharge": q_u,
        "lower_momentum": M_l,
        "upper_momentum": M_u,
    }
    nonfinite = [
        name for name, value in calculated_values.items() if not math.isfinite(value)
    ]
    if nonfinite:
        raise FloatingPointError(
            "The analytical evaluation produced non-finite values: "
            + ", ".join(nonfinite)
            + ". Check the parameter magnitudes; very small power-law indices can "
            "overflow double precision."
        )

    return DimensionlessSolution(
        lambda_0=lam,
        m_lower=m_l,
        m_upper=m_u,
        C_lower=C_l,
        A_lower=A_l,
        B_lower=B_l,
        A_upper=A_u,
        B_upper=B_u,
        Lambda_lower=Lambda_l,
        Lambda_upper=Lambda_u,
        tau_b=tau_b,
        tau_interface=tau_I,
        interface_velocity=U_I,
        upper_velocity_increment=W_u,
        surface_velocity=U_I + W_u,
        lower_mean_velocity=Ubar_l,
        upper_mean_velocity=Ubar_u,
        lower_discharge=q_l,
        upper_discharge=q_u,
        lower_momentum=M_l,
        upper_momentum=M_u,
    )


def dimensionless_profiles(
    p: DimensionlessParameters,
    s: DimensionlessSolution,
) -> Dict[str, np.ndarray]:
    """Return dimensionless velocity, traction, and shear-rate profiles."""

    z_l = np.linspace(0.0, 1.0, POINTS_PER_LAYER)
    local_l = z_l.copy()
    velocity_l = s.interface_velocity * lower_shape(
        local_l, p.lower_index, s.lambda_0
    )
    traction_l = s.tau_b - (s.tau_b - s.tau_interface) * local_l
    shear_l = (
        traction_l / s.Lambda_lower
    ) ** (1.0 / p.lower_index)

    y_u = np.linspace(0.0, p.depth_ratio, POINTS_PER_LAYER)
    local_u = y_u / p.depth_ratio
    z_u = 1.0 + y_u
    velocity_u = s.interface_velocity + s.upper_velocity_increment * upper_shape(
        local_u, p.upper_index
    )
    traction_u = s.tau_interface * (1.0 - local_u)
    shear_u = (
        traction_u / s.Lambda_upper
    ) ** (1.0 / p.upper_index)

    return {
        "z_lower": z_l,
        "z_upper": z_u,
        "local_lower": local_l,
        "local_upper": local_u,
        "velocity_lower": velocity_l,
        "velocity_upper": velocity_u,
        "traction_lower": traction_l,
        "traction_upper": traction_u,
        "shear_rate_lower": shear_l,
        "shear_rate_upper": shear_u,
    }


# =============================================================================
# 5. OUTPUT AND VERIFICATION
# =============================================================================

def verification_metrics(
    p: DimensionlessParameters,
    s: DimensionlessSolution,
    profiles: Dict[str, np.ndarray],
) -> Dict[str, float]:
    """Check analytical moments with high-order Gauss quadrature."""

    nodes, weights = np.polynomial.legendre.leggauss(200)
    xi = 0.5 * (nodes + 1.0)
    wi = 0.5 * weights

    velocity_l = s.interface_velocity * lower_shape(
        xi, p.lower_index, s.lambda_0
    )
    velocity_u = s.interface_velocity + s.upper_velocity_increment * upper_shape(
        xi, p.upper_index
    )

    q_l_num = float(np.sum(wi * velocity_l))
    q_u_num = p.depth_ratio * float(np.sum(wi * velocity_u))
    M_l_num = float(np.sum(wi * velocity_l**2))
    M_u_num = p.depth_ratio * float(np.sum(wi * velocity_u**2))

    source_l = 1.0 / p.froude_lower**2 + s.tau_interface - s.tau_b
    source_u = p.depth_ratio / p.froude_lower**2 - s.tau_interface / p.density_ratio

    lower_rheology_tau = s.Lambda_lower * (
        s.interface_velocity / s.C_lower
    ) ** p.lower_index
    upper_rheology_tau = s.Lambda_upper * (
        s.m_upper * s.upper_velocity_increment / p.depth_ratio
    ) ** p.upper_index

    return {
        "lower_discharge_quadrature_error": abs(q_l_num - s.lower_discharge),
        "upper_discharge_quadrature_error": abs(q_u_num - s.upper_discharge),
        "lower_momentum_quadrature_error": abs(M_l_num - s.lower_momentum),
        "upper_momentum_quadrature_error": abs(M_u_num - s.upper_momentum),
        "velocity_continuity_error": abs(
            profiles["velocity_lower"][-1] - profiles["velocity_upper"][0]
        ),
        "traction_continuity_error": abs(
            profiles["traction_lower"][-1] - profiles["traction_upper"][0]
        ),
        "free_surface_traction": abs(profiles["traction_upper"][-1]),
        "lower_source_residual": abs(source_l),
        "upper_source_residual": abs(source_u),
        "lower_rheology_traction_error": abs(lower_rheology_tau - s.tau_b),
        "upper_rheology_traction_error": abs(upper_rheology_tau - s.tau_interface),
        "Lambda_ratio_error": abs(
            s.Lambda_upper / s.Lambda_lower - p.scaled_consistency_ratio
        ),
    }


def write_profile_table(
    p: DimensionlessParameters,
    profiles: Dict[str, np.ndarray],
    filename: Path,
) -> None:
    """Write one Origin-compatible dimensionless profile file."""

    z = np.concatenate((profiles["z_lower"], profiles["z_upper"][1:]))
    layer = np.concatenate(
        (
            np.zeros_like(profiles["z_lower"]),
            np.ones_like(profiles["z_upper"][1:]),
        )
    )
    local = np.concatenate((profiles["local_lower"], profiles["local_upper"][1:]))
    velocity = np.concatenate(
        (profiles["velocity_lower"], profiles["velocity_upper"][1:])
    )
    traction = np.concatenate(
        (profiles["traction_lower"], profiles["traction_upper"][1:])
    )
    shear_rate = np.concatenate(
        (profiles["shear_rate_lower"], profiles["shear_rate_upper"][1:])
    )

    table = np.column_stack((z, layer, local, velocity, traction, shear_rate))
    np.savetxt(
        filename,
        table,
        delimiter="\t",
        comments="",
        fmt="%.12e",
        header=(
            "z_star\tlayer_0_lower_1_upper\tlocal_coordinate"
            "\tu_star\ttau_star\tdu_star_dz_star"
        ),
    )


def summary_lines(
    p: DimensionlessParameters,
    s: DimensionlessSolution,
    checks: Dict[str, float],
) -> List[str]:
    lines = [
        "EXACT DIMENSIONLESS STEADY-UNIFORM TWO-LAYER POWER-LAW FLOW",
        "=" * 72,
        "",
        "INDEPENDENT INPUTS",
        f"Fr_l                                   = {p.froude_lower:.12g}",
        f"h_r = H_u/H_l                          = {p.depth_ratio:.12g}",
        f"r_rho = rho_u/rho_l                    = {p.density_ratio:.12g}",
        f"n_l                                    = {p.lower_index:.12g}",
        f"n_u                                    = {p.upper_index:.12g}",
        f"kappa_K                                = {p.scaled_consistency_ratio:.12g}",
        "",
        "NORMAL-FLOW TRACTIONS AND RHEOLOGICAL COEFFICIENTS",
        f"lambda_0                               = {s.lambda_0:.12g}",
        f"tau_b^0                                = {s.tau_b:.12g}",
        f"tau_I^0                                = {s.tau_interface:.12g}",
        f"Lambda_l                               = {s.Lambda_lower:.12g}",
        f"Lambda_u                               = {s.Lambda_upper:.12g}",
        "",
        "SHAPE COEFFICIENTS",
        f"m_l                                    = {s.m_lower:.12g}",
        f"C_l(lambda_0)                          = {s.C_lower:.12g}",
        f"A_l(lambda_0)                          = {s.A_lower:.12g}",
        f"B_l(lambda_0)                          = {s.B_lower:.12g}",
        f"m_u                                    = {s.m_upper:.12g}",
        f"a_u                                    = {s.A_upper:.12g}",
        f"b_u                                    = {s.B_upper:.12g}",
        "",
        "DIMENSIONLESS VELOCITIES, DISCHARGES, AND MOMENTUM MOMENTS",
        f"h_l^0                                  = 1",
        f"h_u^0                                  = {p.depth_ratio:.12g}",
        f"U_I^0                                  = {s.interface_velocity:.12g}",
        f"W_u^0                                  = {s.upper_velocity_increment:.12g}",
        f"U_s^0                                  = {s.surface_velocity:.12g}",
        f"mean lower velocity                    = {s.lower_mean_velocity:.12g}",
        f"mean upper velocity                    = {s.upper_mean_velocity:.12g}",
        f"q_l^0                                  = {s.lower_discharge:.12g}",
        f"q_u^0                                  = {s.upper_discharge:.12g}",
        f"M_l^0                                  = {s.lower_momentum:.12g}",
        f"M_u^0                                  = {s.upper_momentum:.12g}",
        "",
        "VERIFICATION",
    ]
    for name, value in checks.items():
        lines.append(f"{name:<44s} = {value:.6e}")
    return lines


def make_plots(
    p: DimensionlessParameters,
    profiles: Dict[str, np.ndarray],
    output_directory: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is unavailable; profile tables were written without plots.")
        return

    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 11.5,
            "axes.labelsize": 12.0,
            "legend.fontsize": 10.5,
        }
    )

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.plot(profiles["velocity_lower"], profiles["z_lower"], label="lower layer")
    ax.plot(profiles["velocity_upper"], profiles["z_upper"], label="upper layer")
    ax.axhline(1.0, linestyle="--", linewidth=0.9, label="interface")
    ax.set_xlabel(r"dimensionless velocity $u/U_l$")
    ax.set_ylabel(r"dimensionless coordinate $z/H_l$")
    ax.legend(frameon=False)
    ax.grid(True, linewidth=0.45, alpha=0.35)
    fig.tight_layout()
    fig.savefig(output_directory / "velocity_profile_dimensionless.pdf")
    fig.savefig(output_directory / "velocity_profile_dimensionless.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.plot(profiles["traction_lower"], profiles["z_lower"], label="lower layer")
    ax.plot(profiles["traction_upper"], profiles["z_upper"], label="upper layer")
    ax.axhline(1.0, linestyle="--", linewidth=0.9, label="interface")
    ax.set_xlabel(r"dimensionless shear traction $\tau/\tau_{ref}$")
    ax.set_ylabel(r"dimensionless coordinate $z/H_l$")
    ax.legend(frameon=False)
    ax.grid(True, linewidth=0.45, alpha=0.35)
    fig.tight_layout()
    fig.savefig(output_directory / "traction_profile_dimensionless.pdf")
    fig.savefig(output_directory / "traction_profile_dimensionless.png", dpi=300)
    plt.close(fig)


def main() -> None:
    solution = solve_dimensionless(PARAMETERS)
    profiles = dimensionless_profiles(PARAMETERS, solution)
    checks = verification_metrics(PARAMETERS, solution, profiles)

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_profile_table(
        PARAMETERS,
        profiles,
        OUTPUT_DIRECTORY / "steady_uniform_dimensionless_profile.txt",
    )

    text = "\n".join(summary_lines(PARAMETERS, solution, checks)) + "\n"
    (OUTPUT_DIRECTORY / "steady_uniform_dimensionless_summary.txt").write_text(
        text, encoding="utf-8"
    )

    if MAKE_PLOTS:
        make_plots(PARAMETERS, profiles, OUTPUT_DIRECTORY)

    print(text)
    print(f"Output directory: {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
