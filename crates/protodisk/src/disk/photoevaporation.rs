//! Photoevaporation models for disk dispersal.
//!
//! High-energy photons (EUV and X-ray) from the central star heat the disk
//! surface, driving mass loss through photoevaporative winds. This is one of
//! the primary mechanisms for disk dispersal.
//!
//! # Physics
//!
//! Two main regimes:
//!
//! - **EUV (Extreme Ultraviolet)**: Photons ionize H atoms, heating gas to
//!   ~10^4 K. Mass loss occurs beyond the gravitational radius r_g where the
//!   thermal velocity exceeds escape velocity.
//!
//! - **X-ray**: Higher energy photons penetrate deeper, heating gas over a
//!   larger radial range. Generally more effective at dispersing disks.
//!
//! # Disk Dispersal Timescale
//!
//! Typical disk lifetimes: 1-10 Myr
//!
//! The process is "two-timescale":
//! 1. Slow viscous decline (several Myr)
//! 2. Rapid inside-out clearing when photoevaporation dominates (~0.1 Myr)
//!
//! # References
//! - Clarke, Gendrin & Sotomayor (2001) - "The dispersal of circumstellar discs"
//! - Owen, Ercolano & Clarke (2011) - "X-ray Photoevaporation"
//! - Alexander, Clarke & Pringle (2006) - "Photoevaporation of protoplanetary discs"

use serde::{Deserialize, Serialize};
use units::{Length, Mass, MassRate};

use crate::disk::constants::{G, K_B, M_PROTON, MU, PI};

/// Photoevaporation model.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub enum PhotoevaporationModel {
    /// No photoevaporation
    #[default]
    None,

    /// EUV photoevaporation only.
    ///
    /// Mass loss beyond gravitational radius r_g where thermal velocity
    /// exceeds escape velocity.
    ///
    /// Typical EUV photon flux: Φ_EUV ~ 10^41 photons/s (T Tauri stars)
    Euv {
        /// EUV photon flux (photons/s)
        phi_euv: f64,
    },

    /// X-ray photoevaporation only.
    ///
    /// More penetrating than EUV, affects larger radial range.
    ///
    /// Typical X-ray luminosity: L_X ~ 10^30 erg/s
    Xray {
        /// X-ray luminosity (erg/s)
        luminosity_xray: f64,
    },

    /// Combined EUV + X-ray photoevaporation.
    Combined {
        /// EUV photon flux (photons/s)
        phi_euv: f64,
        /// X-ray luminosity (erg/s)
        luminosity_xray: f64,
    },
}

impl PhotoevaporationModel {
    /// Create a standard T Tauri star photoevaporation model.
    ///
    /// Typical values:
    /// - EUV flux: 10^41 photons/s
    /// - X-ray luminosity: 10^30 erg/s
    pub fn t_tauri() -> Self {
        Self::Combined {
            phi_euv: 1e41,
            luminosity_xray: 1e30,
        }
    }

    /// Create a weak photoevaporation model (low-mass stars).
    pub fn weak() -> Self {
        Self::Combined {
            phi_euv: 1e40,
            luminosity_xray: 1e29,
        }
    }

    /// Create a strong photoevaporation model (intermediate-mass stars).
    pub fn strong() -> Self {
        Self::Combined {
            phi_euv: 1e42,
            luminosity_xray: 1e31,
        }
    }

    /// Calculate mass loss rate at radius r.
    ///
    /// # Arguments
    /// * `r` - Radius (cm)
    /// * `stellar_mass` - Central star mass (g)
    /// * `surface_density` - Local surface density (g/cm²)
    ///
    /// # Returns
    /// Mass loss rate per unit area (g/cm²/s)
    pub fn mass_loss_rate_per_area(&self, r: f64, stellar_mass: f64, surface_density: f64) -> f64 {
        match self {
            Self::None => 0.0,
            Self::Euv { phi_euv } => self.euv_mass_loss_rate(r, stellar_mass, *phi_euv),
            Self::Xray { luminosity_xray } => {
                self.xray_mass_loss_rate(r, stellar_mass, surface_density, *luminosity_xray)
            }
            Self::Combined {
                phi_euv,
                luminosity_xray,
            } => {
                let euv = self.euv_mass_loss_rate(r, stellar_mass, *phi_euv);
                let xray =
                    self.xray_mass_loss_rate(r, stellar_mass, surface_density, *luminosity_xray);
                euv + xray
            }
        }
    }

    /// EUV photoevaporation mass loss rate.
    ///
    /// Based on Clarke et al. (2001). Mass loss occurs beyond the
    /// gravitational radius r_g where the sound speed at ~10^4 K exceeds
    /// the escape velocity.
    ///
    /// Mdot ~ 4π × r_g² × ρ(r_g) × c_s
    ///
    /// Distributed over the disk beyond r_g.
    fn euv_mass_loss_rate(&self, r: f64, stellar_mass: f64, phi_euv: f64) -> f64 {
        // Gravitational radius: r_g = G M_* / (2 c_s²)
        // For T ~ 10^4 K: c_s ~ 10 km/s = 10^6 cm/s
        let t_euv = 1e4; // K
        let c_s = (K_B * t_euv / (MU * M_PROTON)).sqrt();
        let r_g = G * stellar_mass / (2.0 * c_s * c_s);

        if r < r_g {
            // No EUV mass loss inside gravitational radius
            return 0.0;
        }

        // Mass loss rate scales with EUV flux
        // Normalization: Φ_EUV ~ 10^41 photons/s gives Ṁ ~ 10^-10 M_sun/yr
        // But distributed over outer disk (r > r_g, typically r_g ~ 5 AU)
        let mdot_total = 1e-10 * 1.989e33 / 3.156e7 * (phi_euv / 1e41); // g/s

        // Distribute over disk beyond r_g with 1/r falloff
        // Σ̇(r) = Ṁ × f(r) where ∫ f(r) 2πr dr = 1
        // For f(r) = A/r × exp(-(r-r_g)/r_scale), need normalization

        let r_scale = 10.0 * r_g; // Scale length for EUV wind
        let profile = (-(r - r_g) / r_scale).exp() / r;

        // Normalization: ∫_{r_g}^{∞} (e^(-(r-r_g)/r_scale) / r) 2πr dr ≈ 2π r_scale
        let normalization = mdot_total / (2.0 * PI * r_scale);

        normalization * profile
    }

    /// X-ray photoevaporation mass loss rate.
    ///
    /// Based on Owen et al. (2011). X-rays penetrate deeper and heat
    /// gas over a wider range of radii than EUV.
    ///
    /// Mass loss rate depends on X-ray luminosity and disk structure.
    fn xray_mass_loss_rate(
        &self,
        r: f64,
        stellar_mass: f64,
        surface_density: f64,
        l_x: f64,
    ) -> f64 {
        // X-ray heating temperature ~few × 10^3 K
        let t_xray = 5e3; // K
        let c_s = (K_B * t_xray / (MU * M_PROTON)).sqrt();

        // Critical radius where thermal velocity ~ escape velocity
        // For solar mass: r_crit ~ 1-2 AU
        let r_crit = 0.5 * G * stellar_mass / (c_s * c_s);

        // Mass loss rate scales with X-ray luminosity
        // Normalization: L_X ~ 10^30 erg/s gives Ṁ ~ 10^-8 M_sun/yr
        let mdot_total = 1e-8 * 1.989e33 / 3.156e7 * (l_x / 1e30); // g/s

        // X-ray mass loss peaks at intermediate radii (1-10 AU typically)
        // Use Gaussian profile in log(r) centered at r_crit
        let log_r = r.ln();
        let log_r_peak = r_crit.ln();
        let log_width = 0.8; // Width in log space (covers ~factor of 3)

        let profile = (-(log_r - log_r_peak).powi(2) / (2.0 * log_width * log_width)).exp();

        // Mass loss also depends on column density - need sufficient
        // surface density to absorb X-rays but not too much
        let sigma_optimal = 10.0; // g/cm² - optimal column for X-ray heating
        let column_factor = (surface_density / sigma_optimal).min(1.0);

        // Normalize: ∫ profile × (1/r) × 2πr dr over disk
        // For Gaussian in log-space: integral ~ (2π)^(1/2) × r_peak × log_width × mdot_total
        let normalization = mdot_total / (2.0 * PI * r_crit * log_width * (2.0 * PI).sqrt());

        normalization * profile * column_factor / r
    }
}

/// Photoevaporation-enhanced disk model.
///
/// Wraps a disk model and adds photoevaporation mass loss.
pub struct PhotoevaporatingDisk<D> {
    /// Underlying disk
    pub disk: D,

    /// Photoevaporation model
    pub photoevaporation: PhotoevaporationModel,

    /// Total mass lost to photoevaporation
    pub mass_lost: Mass,
}

impl<D> PhotoevaporatingDisk<D> {
    /// Create a new photoevaporating disk.
    pub fn new(disk: D, photoevaporation: PhotoevaporationModel) -> Self {
        Self {
            disk,
            photoevaporation,
            mass_lost: Mass::zero(),
        }
    }

    /// Mass loss rate at a given radius.
    ///
    /// Returns the rate at which surface density decreases due to
    /// photoevaporation: dΣ/dt (g/cm²/s).
    pub fn mass_loss_rate_at(
        &self,
        r: Length,
        stellar_mass: Mass,
        surface_density: units::SurfaceDensity,
    ) -> f64 {
        self.photoevaporation.mass_loss_rate_per_area(
            r.to_cm(),
            stellar_mass.to_grams(),
            surface_density.to_grams_per_cm2(),
        )
    }

    /// Total photoevaporative mass loss rate from entire disk.
    ///
    /// Integrated over all radii: Ṁ = ∫ Σ̇ 2πr dr
    pub fn total_mass_loss_rate(
        &self,
        stellar_mass: Mass,
        radii: &[Length],
        sigma: &[f64],
    ) -> MassRate {
        let mut mdot = 0.0;

        for i in 0..radii.len() - 1 {
            let r1 = radii[i].to_cm();
            let r2 = radii[i + 1].to_cm();
            let s1 = sigma[i];
            let s2 = sigma[i + 1];

            let r_mid = (r1 * r2).sqrt();
            let s_mid = (s1 * s2).sqrt();
            let dr = r2 - r1;

            let sigma_dot = self.photoevaporation.mass_loss_rate_per_area(
                r_mid,
                stellar_mass.to_grams(),
                s_mid,
            );

            // dM/dt = Σ̇ × 2πr dr
            mdot += sigma_dot * 2.0 * PI * r_mid * dr;
        }

        MassRate::from_grams_per_year(mdot * 3.156e7) // Convert /s to /year
    }
}
