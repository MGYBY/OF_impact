#!/usr/bin/env python3
"""Exact dimensional steady-uniform solution for two power-law layers.

The code implements the analytical normal-flow solution derived in the report
"A Two-Layer Karman-Pohlhausen Shallow-Layer Model for Shear-Thinning
Power-Law Fluids".

Geometry
--------
The lower layer occupies 0 <= z <= H_l and the upper layer occupies
H_l <= z <= H_l + H_u. The x direction is downslope along a rigid plane
inclined by ``theta``. Both layers obey

    tau_i = K_i |du_i/dz|^(n_i - 1) du_i/dz,

with 0 < n_i <= 1 on the positive monotone-shear branch.

The program calculates
----------------------
* exact basal and interfacial tractions;
* lower- and upper-layer velocity profiles;
* interfacial and free-surface velocities;
* layer mean velocities and discharges;
* depth-integrated momentum moments;
* the dimensionless parameters used by the companion dimensionless program;
* tab-delimited profiles and publication-friendly PDF/PNG plots.

Only NumPy is required for the calculation. Matplotlib is optional and is used
only when ``MAKE_PLOTS`` is True.
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
class DimensionalParameters:
    """Physical input data in SI units."""

    gravity: float                  # m/s^2
    slope_angle_degrees: float      # degrees above horizontal
    lower_depth: float              # H_l [m]
    upper_depth: float              # H_u [m]
    lower_density: float            # rho_l [kg/m^3]
    upper_density: float            # rho_u [kg/m^3]
    lower_consistency: float        # K_l [Pa s^n_l]
    upper_consistency: float        # K_u [Pa s^n_u]
    lower_index: float              # n_l
    upper_index: float              # n_u


# Edit this block for a new physical case.
PARAMETERS = DimensionalParameters(
    gravity=9.81,
    slope_angle_degrees=5.0,
    lower_depth=0.010,
    upper_depth=0.008,
    lower_density=1200.0,
    upper_density=1050.0,
    lower_consistency=4.0,
    upper_consistency=2.0,
    lower_index=0.60,
    upper_index=0.80,
)

POINTS_PER_LAYER = 501
MAKE_PLOTS = True
OUTPUT_DIRECTORY = (
    Path(__file__).resolve().parent / "dimensional_steady_uniform_output"
)


# =============================================================================
# 2. DATA STRUCTURE FOR THE ANALYTICAL SOLUTION
# =============================================================================

@dataclass(frozen=True)
class DimensionalSolution:
    theta_radians: float
    slope: float
    g_streamwise: float
    g_normal: float
    tau_b: float
    tau_interface: float
    lambda_0: float
    m_lower: float
    m_upper: float
    C_lower: float
    A_lower: float
    B_lower: float
    A_upper: float
    B_upper: float
    interface_velocity: float
    upper_velocity_increment: float
    surface_velocity: float
    lower_mean_velocity: float
    upper_mean_velocity: float
    lower_discharge: float
    upper_discharge: float
    lower_momentum: float
    upper_momentum: float
    length_scale: float
    time_scale: float
    traction_scale: float
    froude_lower: float
    depth_ratio: float
    density_ratio: float
    scaled_consistency_ratio: float
    Lambda_lower: float
    Lambda_upper: float


# =============================================================================
# 3. ANALYTICAL SHAPE FUNCTIONS
# =============================================================================

def validate_parameters(p: DimensionalParameters) -> None:
    """Check the domain of the analytical positive-shear solution."""

    positive_values = {
        "gravity": p.gravity,
        "lower_depth": p.lower_depth,
        "upper_depth": p.upper_depth,
        "lower_density": p.lower_density,
        "upper_density": p.upper_density,
        "lower_consistency": p.lower_consistency,
        "upper_consistency": p.upper_consistency,
        "lower_index": p.lower_index,
        "upper_index": p.upper_index,
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
    if not math.isfinite(p.slope_angle_degrees):
        raise ValueError("slope_angle_degrees must be finite.")
    if not (0.0 < p.slope_angle_degrees < 90.0):
        raise ValueError(
            "slope_angle_degrees must lie strictly between 0 and 90. "
            "The zero-slope state has zero velocity and makes the chosen normal-flow "
            "length/time scaling singular."
        )
    if POINTS_PER_LAYER < 11:
        raise ValueError("POINTS_PER_LAYER must be at least 11.")


def lower_shape_coefficients(index: float, lam: float) -> Tuple[float, float, float]:
    """Return C_l(lambda), A_l(lambda), and B_l(lambda).

    The direct analytical formulas suffer cancellation when lambda approaches
    one. A fourth-order series is therefore used in that limit.
    """

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
        # Series in eps = 1-lambda. The lambda -> 1 profile is linear.
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
    """Evaluate phi_l(s; lambda), 0 <= s <= 1, stably."""

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
    """Evaluate phi_u(r)=1-(1-r)^m, 0 <= r <= 1."""

    r = np.asarray(r, dtype=float)
    if np.any((r < 0.0) | (r > 1.0)):
        raise ValueError("The upper local coordinate r must lie in [0, 1].")
    m = (index + 1.0) / index
    return 1.0 - (1.0 - r) ** m


# =============================================================================
# 4. EXACT DIMENSIONAL SOLUTION
# =============================================================================

def solve_dimensional(p: DimensionalParameters) -> DimensionalSolution:
    """Evaluate the exact steady-uniform two-layer solution."""

    validate_parameters(p)

    theta = math.radians(p.slope_angle_degrees)
    slope = math.tan(theta)
    g_s = p.gravity * math.sin(theta)
    g_n = p.gravity * math.cos(theta)

    H_l = p.lower_depth
    H_u = p.upper_depth
    rho_l = p.lower_density
    rho_u = p.upper_density

    # Exact linear traction fields from the steady momentum balances.
    tau_I = rho_u * g_s * H_u
    tau_b = g_s * (rho_l * H_l + rho_u * H_u)
    lam = tau_I / tau_b

    m_l = (p.lower_index + 1.0) / p.lower_index
    m_u = (p.upper_index + 1.0) / p.upper_index
    C_l, A_l, B_l = lower_shape_coefficients(p.lower_index, lam)

    # Exact velocity amplitudes.
    U_I = H_l * (tau_b / p.lower_consistency) ** (1.0 / p.lower_index) * C_l
    W_u = H_u / m_u * (tau_I / p.upper_consistency) ** (1.0 / p.upper_index)

    A_u = m_u / (m_u + 1.0)
    B_u = 1.0 - 2.0 / (m_u + 1.0) + 1.0 / (2.0 * m_u + 1.0)

    Ubar_l = U_I * A_l
    Ubar_u = U_I + A_u * W_u
    Q_l = H_l * Ubar_l
    Q_u = H_u * Ubar_u
    M_l = H_l * U_I**2 * B_l
    M_u = H_u * (U_I**2 + 2.0 * A_u * U_I * W_u + B_u * W_u**2)

    # Scales and dimensionless parameters used in the report/PyPDE model.
    L = H_l / slope
    T = L / Ubar_l
    tau_ref = rho_l * slope * Ubar_l**2
    Fr_l = Ubar_l / math.sqrt(g_n * H_l)
    h_r = H_u / H_l
    r_rho = rho_u / rho_l
    kappa_K = (
        p.upper_consistency
        / p.lower_consistency
        * (Ubar_l / H_l) ** (p.upper_index - p.lower_index)
    )
    Lambda_l = (
        p.lower_consistency
        * Ubar_l ** (p.lower_index - 2.0)
        / (rho_l * slope * H_l**p.lower_index)
    )
    Lambda_u = (
        p.upper_consistency
        * Ubar_l ** (p.upper_index - 2.0)
        / (rho_l * slope * H_l**p.upper_index)
    )

    calculated_values = {
        "tau_b": tau_b,
        "tau_interface": tau_I,
        "lambda_0": lam,
        "interface_velocity": U_I,
        "upper_velocity_increment": W_u,
        "lower_mean_velocity": Ubar_l,
        "upper_mean_velocity": Ubar_u,
        "lower_discharge": Q_l,
        "upper_discharge": Q_u,
        "lower_momentum": M_l,
        "upper_momentum": M_u,
        "length_scale": L,
        "time_scale": T,
        "traction_scale": tau_ref,
        "froude_lower": Fr_l,
        "scaled_consistency_ratio": kappa_K,
        "Lambda_lower": Lambda_l,
        "Lambda_upper": Lambda_u,
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

    return DimensionalSolution(
        theta_radians=theta,
        slope=slope,
        g_streamwise=g_s,
        g_normal=g_n,
        tau_b=tau_b,
        tau_interface=tau_I,
        lambda_0=lam,
        m_lower=m_l,
        m_upper=m_u,
        C_lower=C_l,
        A_lower=A_l,
        B_lower=B_l,
        A_upper=A_u,
        B_upper=B_u,
        interface_velocity=U_I,
        upper_velocity_increment=W_u,
        surface_velocity=U_I + W_u,
        lower_mean_velocity=Ubar_l,
        upper_mean_velocity=Ubar_u,
        lower_discharge=Q_l,
        upper_discharge=Q_u,
        lower_momentum=M_l,
        upper_momentum=M_u,
        length_scale=L,
        time_scale=T,
        traction_scale=tau_ref,
        froude_lower=Fr_l,
        depth_ratio=h_r,
        density_ratio=r_rho,
        scaled_consistency_ratio=kappa_K,
        Lambda_lower=Lambda_l,
        Lambda_upper=Lambda_u,
    )


def dimensional_profiles(
    p: DimensionalParameters,
    s: DimensionalSolution,
) -> Dict[str, np.ndarray]:
    """Return exact profile arrays for both layers."""

    z_l = np.linspace(0.0, p.lower_depth, POINTS_PER_LAYER)
    local_l = z_l / p.lower_depth
    tau_l = s.tau_b - (s.tau_b - s.tau_interface) * local_l
    shear_l = (tau_l / p.lower_consistency) ** (1.0 / p.lower_index)
    velocity_l = s.interface_velocity * lower_shape(
        local_l, p.lower_index, s.lambda_0
    )

    y_u = np.linspace(0.0, p.upper_depth, POINTS_PER_LAYER)
    local_u = y_u / p.upper_depth
    z_u = p.lower_depth + y_u
    tau_u = s.tau_interface * (1.0 - local_u)
    shear_u = (tau_u / p.upper_consistency) ** (1.0 / p.upper_index)
    velocity_u = s.interface_velocity + s.upper_velocity_increment * upper_shape(
        local_u, p.upper_index
    )

    # Apparent viscosity K*|gamma_dot|^(n-1).  For a shear-thinning
    # power-law fluid (n < 1), it diverges at the zero-shear free surface.
    # For the Newtonian endpoint n = 1, its limiting value is exactly K.
    def apparent_viscosity(
        consistency: float,
        index: float,
        shear_rate: np.ndarray,
    ) -> np.ndarray:
        if index == 1.0:
            return np.full_like(shear_rate, consistency)
        viscosity = np.full_like(shear_rate, np.inf)
        positive = shear_rate > 0.0
        viscosity[positive] = (
            consistency * shear_rate[positive] ** (index - 1.0)
        )
        return viscosity

    mu_l = apparent_viscosity(
        p.lower_consistency, p.lower_index, shear_l
    )
    mu_u = apparent_viscosity(
        p.upper_consistency, p.upper_index, shear_u
    )

    return {
        "z_lower": z_l,
        "z_upper": z_u,
        "local_lower": local_l,
        "local_upper": local_u,
        "velocity_lower": velocity_l,
        "velocity_upper": velocity_u,
        "traction_lower": tau_l,
        "traction_upper": tau_u,
        "shear_rate_lower": shear_l,
        "shear_rate_upper": shear_u,
        "apparent_viscosity_lower": mu_l,
        "apparent_viscosity_upper": mu_u,
    }


# =============================================================================
# 5. OUTPUT AND VERIFICATION
# =============================================================================

def verification_metrics(
    p: DimensionalParameters,
    s: DimensionalSolution,
    profiles: Dict[str, np.ndarray],
) -> Dict[str, float]:
    """Compare the closed-form moments with high-order Gauss quadrature."""

    nodes, weights = np.polynomial.legendre.leggauss(200)
    xi = 0.5 * (nodes + 1.0)
    wi = 0.5 * weights

    velocity_l = s.interface_velocity * lower_shape(
        xi, p.lower_index, s.lambda_0
    )
    velocity_u = s.interface_velocity + s.upper_velocity_increment * upper_shape(
        xi, p.upper_index
    )

    q_l_num = p.lower_depth * float(np.sum(wi * velocity_l))
    q_u_num = p.upper_depth * float(np.sum(wi * velocity_u))
    M_l_num = p.lower_depth * float(np.sum(wi * velocity_l**2))
    M_u_num = p.upper_depth * float(np.sum(wi * velocity_u**2))

    compatibility = (
        (1.0 + s.density_ratio * s.depth_ratio)
        / s.froude_lower**2
        * (s.A_lower * s.C_lower) ** p.lower_index
    )

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
        "lower_gravity_traction_balance": abs(
            s.tau_b - s.tau_interface - p.lower_density * s.g_streamwise * p.lower_depth
        ),
        "upper_gravity_traction_balance": abs(
            s.tau_interface - p.upper_density * s.g_streamwise * p.upper_depth
        ),
        "Lambda_lower_compatibility_error": abs(s.Lambda_lower - compatibility),
        "Lambda_ratio_error": abs(
            s.Lambda_upper / s.Lambda_lower - s.scaled_consistency_ratio
        ),
    }


def write_profile_table(
    p: DimensionalParameters,
    s: DimensionalSolution,
    profiles: Dict[str, np.ndarray],
    filename: Path,
) -> None:
    """Write one Origin-compatible profile file."""

    # Keep the lower interface point and omit the duplicated upper interface point.
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
    apparent_viscosity = np.concatenate(
        (
            profiles["apparent_viscosity_lower"],
            profiles["apparent_viscosity_upper"][1:],
        )
    )

    table = np.column_stack(
        (
            z,
            z / p.lower_depth,
            layer,
            local,
            velocity,
            velocity / s.lower_mean_velocity,
            traction,
            traction / s.traction_scale,
            shear_rate,
            apparent_viscosity,
        )
    )
    np.savetxt(
        filename,
        table,
        delimiter="\t",
        comments="",
        fmt="%.12e",
        header=(
            "z_m\tz_over_H_lower\tlayer_0_lower_1_upper\tlocal_coordinate"
            "\tu_m_per_s\tu_over_U_lower\ttau_Pa\ttau_over_tau_ref"
            "\tshear_rate_per_s\tapparent_viscosity_Pa_s"
        ),
    )


def summary_lines(
    p: DimensionalParameters,
    s: DimensionalSolution,
    checks: Dict[str, float],
) -> List[str]:
    """Create the human-readable summary."""

    lines = [
        "EXACT DIMENSIONAL STEADY-UNIFORM TWO-LAYER POWER-LAW FLOW",
        "=" * 72,
        "",
        "INPUTS (SI units)",
        f"gravity g                              = {p.gravity:.12g} m/s^2",
        f"slope angle theta                      = {p.slope_angle_degrees:.12g} deg",
        f"lower depth H_l                        = {p.lower_depth:.12g} m",
        f"upper depth H_u                        = {p.upper_depth:.12g} m",
        f"lower density rho_l                    = {p.lower_density:.12g} kg/m^3",
        f"upper density rho_u                    = {p.upper_density:.12g} kg/m^3",
        f"lower consistency K_l                  = {p.lower_consistency:.12g} Pa s^n_l",
        f"upper consistency K_u                  = {p.upper_consistency:.12g} Pa s^n_u",
        f"lower power-law index n_l              = {p.lower_index:.12g}",
        f"upper power-law index n_u              = {p.upper_index:.12g}",
        "",
        "GRAVITY AND TRACTIONS",
        f"S_0 = tan(theta)                       = {s.slope:.12g}",
        f"g_s = g sin(theta)                     = {s.g_streamwise:.12g} m/s^2",
        f"g_n = g cos(theta)                     = {s.g_normal:.12g} m/s^2",
        f"basal traction tau_b                   = {s.tau_b:.12g} Pa",
        f"interfacial traction tau_I             = {s.tau_interface:.12g} Pa",
        f"lambda_0 = tau_I/tau_b                 = {s.lambda_0:.12g}",
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
        "VELOCITIES, DISCHARGES, AND MOMENTUM MOMENTS",
        f"interface velocity U_I                 = {s.interface_velocity:.12g} m/s",
        f"upper increment W_u                    = {s.upper_velocity_increment:.12g} m/s",
        f"free-surface velocity U_s              = {s.surface_velocity:.12g} m/s",
        f"lower mean velocity U_l                = {s.lower_mean_velocity:.12g} m/s",
        f"upper mean velocity U_u                = {s.upper_mean_velocity:.12g} m/s",
        f"lower discharge Q_l                    = {s.lower_discharge:.12g} m^2/s",
        f"upper discharge Q_u                    = {s.upper_discharge:.12g} m^2/s",
        f"lower momentum moment M_l              = {s.lower_momentum:.12g} m^3/s^2",
        f"upper momentum moment M_u              = {s.upper_momentum:.12g} m^3/s^2",
        "",
        "REPORT SCALES AND DIMENSIONLESS PARAMETERS",
        f"H = H_l                                = {p.lower_depth:.12g} m",
        f"U = U_l                                = {s.lower_mean_velocity:.12g} m/s",
        f"L = H_l/S_0                            = {s.length_scale:.12g} m",
        f"T = L/U_l                              = {s.time_scale:.12g} s",
        f"tau_ref = rho_l S_0 U_l^2              = {s.traction_scale:.12g} Pa",
        f"Fr_l                                   = {s.froude_lower:.12g}",
        f"h_r = H_u/H_l                          = {s.depth_ratio:.12g}",
        f"r_rho = rho_u/rho_l                    = {s.density_ratio:.12g}",
        f"kappa_K                                = {s.scaled_consistency_ratio:.12g}",
        f"Lambda_l                               = {s.Lambda_lower:.12g}",
        f"Lambda_u                               = {s.Lambda_upper:.12g}",
        f"dimensionless U_I                      = {s.interface_velocity/s.lower_mean_velocity:.12g}",
        f"dimensionless W_u                      = {s.upper_velocity_increment/s.lower_mean_velocity:.12g}",
        f"dimensionless U_u                      = {s.upper_mean_velocity/s.lower_mean_velocity:.12g}",
        f"dimensionless q_l                      = {s.lower_discharge/(p.lower_depth*s.lower_mean_velocity):.12g}",
        f"dimensionless q_u                      = {s.upper_discharge/(p.lower_depth*s.lower_mean_velocity):.12g}",
        "",
        "VALUES TO COPY INTO THE DIMENSIONLESS PROGRAM",
        f"FR_LOWER = {s.froude_lower:.16g}",
        f"DEPTH_RATIO = {s.depth_ratio:.16g}",
        f"DENSITY_RATIO = {s.density_ratio:.16g}",
        f"N_LOWER = {p.lower_index:.16g}",
        f"N_UPPER = {p.upper_index:.16g}",
        f"SCALED_CONSISTENCY_RATIO = {s.scaled_consistency_ratio:.16g}",
        "",
        "VERIFICATION",
    ]
    for name, value in checks.items():
        lines.append(f"{name:<44s} = {value:.6e}")
    return lines


def make_plots(
    p: DimensionalParameters,
    profiles: Dict[str, np.ndarray],
    output_directory: Path,
) -> None:
    """Write separate velocity and traction figures."""

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
    ax.axhline(p.lower_depth, linestyle="--", linewidth=0.9, label="interface")
    ax.set_xlabel(r"velocity $u$ (m s$^{-1}$)")
    ax.set_ylabel(r"normal coordinate $z$ (m)")
    ax.legend(frameon=False)
    ax.grid(True, linewidth=0.45, alpha=0.35)
    fig.tight_layout()
    fig.savefig(output_directory / "velocity_profile_dimensional.pdf")
    fig.savefig(output_directory / "velocity_profile_dimensional.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.plot(profiles["traction_lower"], profiles["z_lower"], label="lower layer")
    ax.plot(profiles["traction_upper"], profiles["z_upper"], label="upper layer")
    ax.axhline(p.lower_depth, linestyle="--", linewidth=0.9, label="interface")
    ax.set_xlabel(r"shear traction $\tau$ (Pa)")
    ax.set_ylabel(r"normal coordinate $z$ (m)")
    ax.legend(frameon=False)
    ax.grid(True, linewidth=0.45, alpha=0.35)
    fig.tight_layout()
    fig.savefig(output_directory / "traction_profile_dimensional.pdf")
    fig.savefig(output_directory / "traction_profile_dimensional.png", dpi=300)
    plt.close(fig)


def main() -> None:
    """Calculate, verify, and write the exact dimensional solution."""

    solution = solve_dimensional(PARAMETERS)
    profiles = dimensional_profiles(PARAMETERS, solution)
    checks = verification_metrics(PARAMETERS, solution, profiles)

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_profile_table(
        PARAMETERS,
        solution,
        profiles,
        OUTPUT_DIRECTORY / "steady_uniform_dimensional_profile.txt",
    )

    text = "\n".join(summary_lines(PARAMETERS, solution, checks)) + "\n"
    (OUTPUT_DIRECTORY / "steady_uniform_dimensional_summary.txt").write_text(
        text, encoding="utf-8"
    )

    if MAKE_PLOTS:
        make_plots(PARAMETERS, profiles, OUTPUT_DIRECTORY)

    print(text)
    print(f"Output directory: {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
