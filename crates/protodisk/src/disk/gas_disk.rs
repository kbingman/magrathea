//! Gas disk structure and evolution
//!
//! Provides a protoplanetary gas disk with power-law surface density and
//! temperature profiles. All other quantities (scale height, midplane density,
//! pressure, etc.) are derived from these.
//!
//! # Physics
//!
//! The disk is assumed to be vertically isothermal and in hydrostatic equilibrium.
//! Key relations:
//!
//! - Scale height: h = c_s / Ω_K
//! - Midplane density: ρ = Σ / (√(2π) h)
//! - Pressure: P = ρ c_s²
//! - Pressure gradient parameter: η = -(h/r)² × (1/2) × d ln P / d ln r
//!
//! The pressure gradient parameter η controls the sub-Keplerian rotation of
//! gas and drives radial drift of solid particles.

use serde::{Deserialize, Serialize};
use units::{Length, Mass, SurfaceDensity, Temperature};

use crate::{MainSequenceStar, disk::constants::PI, solar_analog};

use super::disk_model::{DiskMass, DiskModel};

/// A protoplanetary gas disk with power-law profiles.
///
/// Surface density: Σ(r) = Σ_0 × (r / r_0)^(-p)
/// Temperature: T(r) = T_0 × (r / r_0)^(-q)
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GasDisk {
    /// Inner edge of the disk
    pub inner_radius: Length,
    /// Outer edge of the disk
    pub outer_radius: Length,

    /// Reference radius (typically 1 AU)
    pub r_0: Length,
    /// Surface density at r_0
    pub sigma_0: SurfaceDensity,
    /// Surface density power-law exponent (Σ ∝ r^(-p))
    /// Typical value: 1.0 (MMSN) to 1.5
    pub sigma_exponent: f64,

    /// Temperature at r_0
    pub temperature_0: Temperature,
    /// Temperature power-law exponent (T ∝ r^(-q))
    /// Typical value: 0.5 (flared disk) to 0.75 (flat disk)
    pub temp_exponent: f64,

    /// Central stellar mass
    pub stellar_mass: Mass,

    /// Shakura-Sunyaev viscosity parameter
    /// Typical value: 1e-3 to 1e-2
    pub alpha: f64,
}

impl GasDisk {
    /// Create a disk scaled to stellar properties.
    ///
    /// Temperature scaling for a passively irradiated disk:
    /// T(r) ∝ L*^(1/4) × r^(-1/2)
    ///
    /// At 1 AU around the Sun, T ≈ 280K. For other stars:
    /// T_0 = 280K × (L/L☉)^(1/4)
    ///
    /// Surface density is scaled to maintain similar disk-to-star mass ratio:
    /// Σ_0 = 1700 g/cm² × (M/M☉)
    pub fn for_star(star: &MainSequenceStar) -> Self {
        let mass_solar = star.mass.to_solar_masses();
        let luminosity_solar = star.luminosity; // already in solar luminosities

        // Temperature at 1 AU scales as L^(1/4)
        // Solar value: 280K at 1 AU
        let temp_at_1au = 280.0 * luminosity_solar.powf(0.25);

        // Surface density scales with stellar mass to maintain
        // similar disk-to-star mass ratio
        // Solar value: ~1700 g/cm² at 1 AU
        let sigma_at_1au = 1700.0 * mass_solar;

        // Disk extent scales with stellar mass (approximation)
        // More massive stars have more extended disks
        let inner = 0.1 * mass_solar.powf(0.5);
        let outer = 100.0 * mass_solar.powf(0.5);

        Self {
            inner_radius: Length::from_au(inner),
            outer_radius: Length::from_au(outer),
            r_0: Length::from_au(1.0),
            sigma_0: SurfaceDensity::from_grams_per_cm2(sigma_at_1au),
            sigma_exponent: 1.0,
            temperature_0: Temperature::from_kelvin(temp_at_1au),
            temp_exponent: 0.5,
            stellar_mass: star.mass,
            alpha: 1e-3,
        }
    }

    /// Create a Minimum Mass Solar Nebula disk around a solar-mass star.
    ///
    /// MMSN parameters from Hayashi (1981):
    /// - Σ = 1700 × (r/AU)^(-1.5) g/cm² (gas)
    /// - T = 280 × (r/AU)^(-0.5) K
    ///
    /// Note: We use a slightly shallower Σ exponent (1.0) which is more
    /// consistent with observed disks.
    pub fn mmsn() -> Self {
        Self::for_star(&solar_analog())
    }

    /// Create a disk with custom parameters.
    pub fn new(
        stellar_mass: Mass,
        sigma_0: SurfaceDensity,
        sigma_exponent: f64,
        temperature_0: Temperature,
        temp_exponent: f64,
        alpha: f64,
    ) -> Self {
        Self {
            inner_radius: Length::from_au(0.1),
            outer_radius: Length::from_au(100.0),
            r_0: Length::from_au(1.0),
            sigma_0,
            sigma_exponent,
            temperature_0,
            temp_exponent,
            stellar_mass,
            alpha,
        }
    }
}

// =============================================================================
// DiskModel trait implementation
// =============================================================================

impl DiskModel for GasDisk {
    fn surface_density(&self, r: Length) -> SurfaceDensity {
        let ratio = r.to_cm() / self.r_0.to_cm();
        SurfaceDensity::from_grams_per_cm2(
            self.sigma_0.to_grams_per_cm2() * ratio.powf(-self.sigma_exponent),
        )
    }

    fn temperature(&self, r: Length) -> Temperature {
        let ratio = r.to_cm() / self.r_0.to_cm();
        Temperature::from_kelvin(self.temperature_0.to_kelvin() * ratio.powf(-self.temp_exponent))
    }

    fn stellar_mass(&self) -> Mass {
        self.stellar_mass
    }

    fn alpha(&self) -> f64 {
        self.alpha
    }

    fn inner_radius(&self) -> Length {
        self.inner_radius
    }

    fn outer_radius(&self) -> Length {
        self.outer_radius
    }

    /// Analytical pressure gradient for power-law disk.
    ///
    /// For Σ ∝ r^(-p) and T ∝ r^(-q):
    /// d ln P / d ln r = -(p + (3+q)/2)
    fn pressure_gradient_log(&self, _r: Length) -> f64 {
        -(self.sigma_exponent + (3.0 + self.temp_exponent) / 2.0)
    }
}

impl DiskMass for GasDisk {
    /// Total disk mass using analytical integration.
    ///
    /// M_disk = ∫ 2πr Σ(r) dr
    ///
    /// For Σ ∝ r^(-p):
    /// - p ≠ 2: M = 2π Σ_0 r_0^p × (r_out^(2-p) - r_in^(2-p)) / (2-p)
    /// - p = 2: M = 2π Σ_0 r_0² × ln(r_out/r_in)
    fn total_mass(&self) -> Mass {
        let sigma_0 = self.sigma_0.to_grams_per_cm2();
        let r_0 = self.r_0.to_cm();
        let r_in = self.inner_radius.to_cm();
        let r_out = self.outer_radius.to_cm();
        let p = self.sigma_exponent;

        let mass = if (p - 2.0).abs() < 1e-10 {
            // p ≈ 2: logarithmic integral
            2.0 * PI * sigma_0 * r_0.powi(2) * (r_out / r_in).ln()
        } else {
            // General case
            let factor = 2.0 * PI * sigma_0 * r_0.powf(p) / (2.0 - p);
            factor * (r_out.powf(2.0 - p) - r_in.powf(2.0 - p))
        };

        Mass::from_grams(mass)
    }
}
