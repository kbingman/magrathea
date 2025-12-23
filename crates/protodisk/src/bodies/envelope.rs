//! Gas envelope accretion physics for protoplanets.
//!
//! Implements hydrostatic envelope solutions and runaway gas accretion based
//! on core accretion theory.
//!
//! # Physics
//!
//! A protoplanet's core can capture a gas envelope from the surrounding disk.
//! The envelope evolves through phases:
//!
//! 1. **Attached phase**: Small envelope in hydrostatic equilibrium
//! 2. **Detached phase**: Envelope extends to Hill radius or Bondi radius
//! 3. **Critical mass**: When M_env ≈ M_core, hydrostatic solution breaks down
//! 4. **Runaway accretion**: Rapid gas capture (10^5-10^6 years)
//!
//! # Critical Core Mass
//!
//! The critical core mass M_crit depends on:
//! - Planetesimal accretion rate (heating from accretion)
//! - Grain opacity (cooling efficiency)
//! - Distance from star (background temperature)
//!
//! Typical values: 5-20 M⊕
//!
//! # References
//! - Pollack et al. (1996) - "Formation of the Giant Planets"
//! - Ikoma et al. (2000) - "Formation of Giant Planets"
//! - Piso & Youdin (2014) - "On the minimum core mass for giant planets"

use units::{Length, Mass, MassRate, SurfaceDensity, Temperature, Time};

use crate::disk::DiskModel;

/// Gravitational constant (cgs units)
const G_CGS: f64 = 6.674e-8; // cm³ g⁻¹ s⁻²

/// Boltzmann constant (erg/K)
const K_B: f64 = 1.380649e-16;

/// Mean molecular weight for H₂-dominated gas
const MU: f64 = 2.3;

/// Proton mass (g)
const M_PROTON: f64 = 1.6726e-24;

/// Calculate critical core mass for runaway gas accretion.
///
/// The critical core mass is where the envelope mass equals the core mass
/// and hydrostatic equilibrium breaks down. Above this mass, runaway
/// gas accretion begins.
///
/// This uses the Piso & Youdin (2014) analytic approximation:
///
/// ```text
/// M_crit ≈ 10 M⊕ × (Ṁ_core / 10^-6 M⊕/yr)^0.25 × (κ / 1 cm²/g)^0.25
/// ```
///
/// where:
/// - Ṁ_core is the planetesimal accretion rate (heating source)
/// - κ is the grain opacity (cooling efficiency)
///
/// # Arguments
/// * `core_accretion_rate` - Rate of planetesimal accretion onto core
/// * `opacity` - Grain opacity in cm²/g (typically 0.01-10)
///
/// # Returns
/// Critical core mass in Earth masses
///
/// # References
/// - Piso & Youdin (2014) - "On the minimum core mass for giant planets"
/// - Rafikov (2006) - "Atmospheres of Protoplanetary Cores"
pub fn critical_core_mass(core_accretion_rate: MassRate, opacity: f64) -> Mass {
    // Normalization values
    let mdot_norm = 1e-6; // M⊕/yr
    let kappa_norm = 1.0; // cm²/g
    let m_crit_base = 10.0; // M⊕

    // Convert accretion rate to M⊕/yr
    let mdot_earth_per_year = core_accretion_rate.to_earth_masses_per_myr() / 1e6;

    // Scaling relations
    let mdot_factor = (mdot_earth_per_year / mdot_norm).powf(0.25);
    let kappa_factor = (opacity / kappa_norm).powf(0.25);

    let m_crit = m_crit_base * mdot_factor * kappa_factor;

    Mass::from_earth_masses(m_crit)
}

/// Calculate Kelvin-Helmholtz contraction timescale for envelope growth.
///
/// Before reaching the critical mass, the envelope grows slowly via
/// Kelvin-Helmholtz (KH) contraction - the envelope cools and contracts,
/// allowing more gas to be captured.
///
/// The KH timescale is:
///
/// ```text
/// τ_KH = G M_core M_env / (R_core L_core)
/// ```
///
/// where L_core is the core's luminosity from planetesimal accretion.
///
/// The envelope growth rate is approximately:
///
/// ```text
/// dM_env/dt ≈ M_env / τ_KH
/// ```
///
/// # Arguments
/// * `core_mass` - Mass of solid core
/// * `envelope_mass` - Current envelope mass
/// * `core_radius` - Physical radius of core
/// * `core_luminosity` - Accretion luminosity (from planetesimal bombardment)
///
/// # Returns
/// KH contraction timescale in years
///
/// # References
/// - Bodenheimer & Pollack (1986) - "Calculations of the accretion and evolution"
/// - Ikoma et al. (2000) - "Formation of Giant Planets"
pub fn kelvin_helmholtz_timescale(
    core_mass: Mass,
    envelope_mass: Mass,
    core_radius: Length,
    core_luminosity: f64, // erg/s
) -> Time {
    let m_core_g = core_mass.to_grams();
    let m_env_g = envelope_mass.to_grams();
    let r_core_cm = core_radius.to_cm();

    // τ_KH = G M_core M_env / (R_core L_core)
    // Result in seconds
    let tau_kh_seconds = (G_CGS * m_core_g * m_env_g) / (r_core_cm * core_luminosity);

    // Convert to years
    Time::from_years(tau_kh_seconds / 3.156e7)
}

/// Calculate envelope growth rate from Kelvin-Helmholtz contraction.
///
/// The envelope mass grows as it cools and contracts:
///
/// ```text
/// dM_env/dt = M_env / τ_KH
/// ```
///
/// This growth is slow (Myr timescales) compared to runaway accretion.
///
/// # Arguments
/// * `core_mass` - Mass of solid core
/// * `envelope_mass` - Current envelope mass
/// * `core_radius` - Physical radius of core
/// * `core_accretion_rate` - Planetesimal accretion rate (for luminosity)
///
/// # Returns
/// Envelope growth rate in solar masses per year
///
/// # Physics
///
/// The core luminosity from planetesimal accretion is:
/// ```text
/// L = G M_core Ṁ_core / R_core
/// ```
pub fn kelvin_helmholtz_growth_rate(
    core_mass: Mass,
    envelope_mass: Mass,
    core_radius: Length,
    core_accretion_rate: MassRate,
) -> MassRate {
    // Calculate core luminosity from planetesimal bombardment
    let m_core_g = core_mass.to_grams();
    let r_core_cm = core_radius.to_cm();
    let mdot_core_g_per_s = core_accretion_rate.to_grams_per_year() / 3.156e7;

    // L_core = G M_core Ṁ_core / R_core (erg/s)
    let l_core = (G_CGS * m_core_g * mdot_core_g_per_s) / r_core_cm;

    // Get KH timescale
    let tau_kh = kelvin_helmholtz_timescale(core_mass, envelope_mass, core_radius, l_core);

    // Growth rate: dM_env/dt = M_env / τ_KH
    let growth_rate = envelope_mass.to_solar_masses() / tau_kh.to_years();

    MassRate::from_solar_masses_per_year(growth_rate)
}

/// Calculate Bondi radius for gas capture.
///
/// The Bondi radius is the radius at which gas is gravitationally bound
/// to the protoplanet. It depends on the sound speed of the surrounding gas:
///
/// ```text
/// R_B = G M / c_s²
/// ```
///
/// The envelope cannot extend beyond the minimum of the Bondi radius
/// and the Hill radius.
///
/// # Arguments
/// * `mass` - Total mass of protoplanet (core + envelope)
/// * `sound_speed` - Sound speed of disk gas (cm/s)
///
/// # Returns
/// Bondi radius in AU
///
/// # References
/// - Bondi (1952) - "On spherically symmetrical accretion"
pub fn bondi_radius(mass: Mass, sound_speed: f64) -> Length {
    let m_g = mass.to_grams();
    let r_bondi_cm = (G_CGS * m_g) / (sound_speed * sound_speed);

    Length::from_cm(r_bondi_cm)
}

/// Calculate gas sound speed from temperature.
///
/// For an ideal gas:
/// ```text
/// c_s = sqrt(k_B T / (μ m_H))
/// ```
///
/// # Arguments
/// * `temperature` - Gas temperature in Kelvin
///
/// # Returns
/// Sound speed in cm/s
pub fn gas_sound_speed(temperature: Temperature) -> f64 {
    let t_k = temperature.to_kelvin();
    ((K_B * t_k) / (MU * M_PROTON)).sqrt()
}

/// Check if envelope can form at this location.
///
/// Requires minimum core mass (~0.1 M⊕) and sufficient gas density.
///
/// # Arguments
/// * `core_mass` - Mass of solid core
/// * `gas_surface_density` - Surface density of gas disk
///
/// # Returns
/// true if envelope capture is possible
pub fn can_capture_envelope(core_mass: Mass, gas_surface_density: SurfaceDensity) -> bool {
    let min_core = Mass::from_earth_masses(0.1);
    let min_sigma = SurfaceDensity::from_grams_per_cm2(1.0); // Very low threshold

    core_mass > min_core && gas_surface_density > min_sigma
}

/// Calculate supply-limited gas accretion rate.
///
/// When in runaway accretion, the rate is limited by how fast gas can
/// flow through the disk to reach the planet. This is approximately:
///
/// ```text
/// Ṁ_gas ≈ 3π × ν × Σ_gas
/// ```
///
/// where ν is the disk viscosity.
///
/// # Arguments
/// * `disk` - Gas disk model providing viscosity and surface density
/// * `location` - Radial location of protoplanet
///
/// # Returns
/// Supply-limited accretion rate in solar masses per year
///
/// # References
/// - Lubow & D'Angelo (2006) - "Gas Flow Across Gaps"
pub fn supply_limited_accretion_rate<D: DiskModel>(disk: &D, location: Length) -> MassRate {
    let sigma = disk.surface_density(location);
    let nu = disk.viscosity(location);

    // Ṁ = 3π × ν × Σ
    // In cgs: g/year = (cm²/s) × (g/cm²) × (s/year)
    let mdot_cgs_per_sec = 3.0 * std::f64::consts::PI * nu * sigma.to_grams_per_cm2();
    let mdot_cgs_per_year = mdot_cgs_per_sec * 3.156e7;

    MassRate::from_grams_per_year(mdot_cgs_per_year)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn critical_mass_increases_with_accretion_rate() {
        let mdot_low = MassRate::from_earth_masses_per_myr(1e-7 * 1e6);
        let mdot_high = MassRate::from_earth_masses_per_myr(1e-5 * 1e6);
        let opacity = 1.0;

        let m_crit_low = critical_core_mass(mdot_low, opacity);
        let m_crit_high = critical_core_mass(mdot_high, opacity);

        assert!(m_crit_high > m_crit_low);
    }

    #[test]
    fn critical_mass_increases_with_opacity() {
        let mdot = MassRate::from_earth_masses_per_myr(1e-6 * 1e6);
        let kappa_low = 0.1;
        let kappa_high = 10.0;

        let m_crit_low = critical_core_mass(mdot, kappa_low);
        let m_crit_high = critical_core_mass(mdot, kappa_high);

        assert!(m_crit_high > m_crit_low);
    }

    #[test]
    fn critical_mass_typical_range() {
        // Typical conditions
        let mdot = MassRate::from_earth_masses_per_myr(1e-6 * 1e6);
        let opacity = 1.0;

        let m_crit = critical_core_mass(mdot, opacity);

        // Should be 5-20 M⊕ range
        assert!(m_crit.to_earth_masses() > 5.0);
        assert!(m_crit.to_earth_masses() < 20.0);
    }

    #[test]
    fn kelvin_helmholtz_timescale_positive() {
        let core_mass = Mass::from_earth_masses(5.0);
        let envelope_mass = Mass::from_earth_masses(1.0);
        let core_radius = Length::from_earth_radii(2.0);
        let core_luminosity = 1e27; // erg/s

        let tau_kh = kelvin_helmholtz_timescale(
            core_mass,
            envelope_mass,
            core_radius,
            core_luminosity,
        );

        assert!(tau_kh.to_years() > 0.0);
    }

    #[test]
    fn kh_timescale_is_long() {
        // Typical sub-Neptune conditions
        let core_mass = Mass::from_earth_masses(5.0);
        let envelope_mass = Mass::from_earth_masses(0.5);
        let core_radius = Length::from_earth_radii(2.0);
        let core_luminosity = 1e27; // erg/s

        let tau_kh = kelvin_helmholtz_timescale(
            core_mass,
            envelope_mass,
            core_radius,
            core_luminosity,
        );

        // KH timescale should be ~Myr
        let tau_myr = tau_kh.to_years() / 1e6;
        assert!(tau_myr > 0.1 && tau_myr < 100.0);
    }

    #[test]
    fn kh_growth_rate_independent_of_envelope_mass() {
        let core_mass = Mass::from_earth_masses(5.0);
        let core_radius = Length::from_earth_radii(2.0);
        let core_accretion_rate = MassRate::from_earth_masses_per_myr(1e-6 * 1e6);

        let env_small = Mass::from_earth_masses(0.1);
        let env_large = Mass::from_earth_masses(1.0);

        let rate_small = kelvin_helmholtz_growth_rate(
            core_mass,
            env_small,
            core_radius,
            core_accretion_rate,
        );
        let rate_large = kelvin_helmholtz_growth_rate(
            core_mass,
            env_large,
            core_radius,
            core_accretion_rate,
        );

        // KH growth rate = M_env / τ_KH = (R_core L_core) / (G M_core)
        // It's actually independent of M_env! (τ_KH ∝ M_env)
        let ratio = rate_large.to_solar_masses_per_year() / rate_small.to_solar_masses_per_year();
        assert!((ratio - 1.0).abs() < 0.01);
    }

    #[test]
    fn bondi_radius_scales_with_mass() {
        let sound_speed = 1e5; // cm/s (typical disk)

        let m_small = Mass::from_earth_masses(1.0);
        let m_large = Mass::from_earth_masses(10.0);

        let r_b_small = bondi_radius(m_small, sound_speed);
        let r_b_large = bondi_radius(m_large, sound_speed);

        // Bondi radius ∝ M
        assert!(r_b_large.to_au() > r_b_small.to_au());
        assert!((r_b_large.to_au() / r_b_small.to_au() - 10.0).abs() < 0.1);
    }

    #[test]
    fn gas_sound_speed_scales_with_temperature() {
        let t_cold = Temperature::from_kelvin(150.0);
        let t_hot = Temperature::from_kelvin(600.0);

        let cs_cold = gas_sound_speed(t_cold);
        let cs_hot = gas_sound_speed(t_hot);

        // c_s ∝ √T
        assert!(cs_hot > cs_cold);
        let ratio = cs_hot / cs_cold;
        let expected_ratio = (600.0_f64 / 150.0_f64).sqrt();
        assert!((ratio - expected_ratio).abs() / expected_ratio < 0.01);
    }

    #[test]
    fn can_capture_requires_minimum_mass() {
        let sigma = SurfaceDensity::from_grams_per_cm2(100.0);

        let too_small = Mass::from_earth_masses(0.01);
        let sufficient = Mass::from_earth_masses(0.5);

        assert!(!can_capture_envelope(too_small, sigma));
        assert!(can_capture_envelope(sufficient, sigma));
    }

    #[test]
    fn can_capture_requires_gas_present() {
        let core_mass = Mass::from_earth_masses(1.0);

        let no_gas = SurfaceDensity::from_grams_per_cm2(0.1);
        let gas_present = SurfaceDensity::from_grams_per_cm2(100.0);

        assert!(!can_capture_envelope(core_mass, no_gas));
        assert!(can_capture_envelope(core_mass, gas_present));
    }
}
