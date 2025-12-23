//! Type I planetary migration physics.
//!
//! Implements the physics of Type I migration, where low-mass planets
//! embedded in gaseous disks experience torques from spiral density waves,
//! causing them to drift inward (or outward at specific locations).
//!
//! # Physical Model
//!
//! Type I migration timescale depends on:
//! - Planet mass (more massive → slower migration)
//! - Orbital radius (closer to star → faster migration)
//! - Disk surface density (denser disk → faster migration)
//! - Disk aspect ratio (thicker disk → slower migration)
//!
//! Migration velocity: v_mig = -r/τ_mig
//!
//! # References
//! - Ward (1997) - "Protoplanet Migration by Nebula Tides"
//! - Tanaka et al. (2002) - "Three-Dimensional Interaction between a Planet and Disk"
//! - Paardekooper et al. (2011) - "A torque formula for non-isothermal Type I migration"

use units::{Length, Mass, SurfaceDensity, Time};

/// Calculate Type I migration timescale.
///
/// The migration timescale τ_mig determines how quickly a planet drifts
/// through the disk due to gravitational torques from spiral density waves.
///
/// # Formula
///
/// ```text
/// τ_mig = (M_*/M_p) × (M_*/Σr²) × (h/r)² × Ω⁻¹
/// ```
///
/// Where:
/// - M_* is stellar mass
/// - M_p is planet mass
/// - Σ is disk surface density
/// - r is orbital radius
/// - h/r is disk aspect ratio
/// - Ω is Keplerian frequency
///
/// # Arguments
/// * `planet_mass` - Mass of the planet
/// * `stellar_mass` - Mass of the central star
/// * `orbital_radius` - Current orbital radius
/// * `surface_density` - Gas surface density at this location
/// * `aspect_ratio` - Disk aspect ratio (H/r) at this location
///
/// # Returns
/// Migration timescale in years
///
/// # Physics
///
/// Type I migration is driven by Lindblad and corotation torques from
/// spiral density waves launched by the planet. For typical disk parameters,
/// migration is inward with timescales of 10⁴-10⁶ years.
///
/// At special locations (e.g., ice lines, pressure maxima), corotation
/// torques can exceed Lindblad torques, causing outward migration and
/// creating "planet traps".
///
/// # References
/// - Tanaka et al. (2002) - Canonical Type I migration formula
/// - Ward (1997) - Original migration analysis
pub fn type_i_migration_timescale(
    planet_mass: Mass,
    stellar_mass: Mass,
    orbital_radius: Length,
    surface_density: SurfaceDensity,
    aspect_ratio: f64,
) -> Time {
    let m_p = planet_mass.to_solar_masses();
    let m_star = stellar_mass.to_solar_masses();
    let r_au = orbital_radius.to_au();
    let sigma_gcm2 = surface_density.to_grams_per_cm2();

    // Keplerian frequency: Ω = √(GM/r³)
    // In units where G = 1 (AU, M☉, years)
    let omega = (m_star / r_au.powi(3)).sqrt();

    // Surface density normalized to M☉/AU²
    // Conversion: 1 g/cm² ≈ 1.4e-7 M☉/AU²
    let sigma_normalized = sigma_gcm2 * 1.4e-7;

    // Type I timescale formula
    let tau_years =
        (m_star / m_p) * (m_star / (sigma_normalized * r_au.powi(2))) * aspect_ratio.powi(2)
            / omega;

    Time::from_years(tau_years)
}

/// Calculate migration velocity from migration timescale.
///
/// For Type I migration, the velocity is:
/// ```text
/// v_mig = dr/dt = -r/τ_mig
/// ```
///
/// The negative sign indicates inward migration (decreasing radius).
///
/// # Arguments
/// * `orbital_radius` - Current orbital radius
/// * `migration_timescale` - Migration timescale from `type_i_migration_timescale()`
///
/// # Returns
/// Migration velocity in AU/year (negative for inward migration)
///
/// # Physics
///
/// This exponential migration law means that closer planets migrate faster
/// in absolute terms, but the fractional migration rate (1/τ) is determined
/// by disk and planet properties.
///
/// # References
/// - Tanaka et al. (2002) - Migration rate formula
pub fn migration_velocity(orbital_radius: Length, migration_timescale: Time) -> f64 {
    // v = -r/τ (negative for inward migration)
    -orbital_radius.to_au() / migration_timescale.to_years()
}

/// Calculate migration distance over a timestep.
///
/// For small timesteps (Δt << τ), this is approximately:
/// ```text
/// Δr ≈ -r × (Δt/τ)
/// ```
///
/// For larger timesteps, we use the exact exponential solution:
/// ```text
/// Δr = -r × (1 - exp(-Δt/τ))
/// ```
///
/// # Arguments
/// * `orbital_radius` - Current orbital radius
/// * `migration_timescale` - Migration timescale
/// * `timestep` - Duration of timestep
///
/// # Returns
/// Change in orbital radius (negative for inward migration)
///
/// # Physics
///
/// The exponential form prevents unphysical jumps when Δt ≈ τ and
/// correctly captures that a planet can't migrate more than its current
/// orbital radius in a single step.
pub fn migration_distance(
    orbital_radius: Length,
    migration_timescale: Time,
    timestep: Time,
) -> f64 {
    let r = orbital_radius.to_au();
    let tau = migration_timescale.to_years();
    let dt = timestep.to_years();

    // Exponential decay: Δr = -r × (1 - exp(-Δt/τ))
    let dt_over_tau = dt / tau;
    let distance = -r * (1.0 - (-dt_over_tau).exp());

    // Cap at 0.5 AU per step to prevent extreme jumps
    distance.max(-0.5)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn migration_timescale_positive() {
        let tau = type_i_migration_timescale(
            Mass::from_earth_masses(5.0),
            Mass::from_solar_masses(1.0),
            Length::from_au(3.0),
            SurfaceDensity::from_grams_per_cm2(100.0),
            0.05,
        );

        assert!(tau.to_years() > 0.0);
        println!("Migration timescale: {:.2e} years", tau.to_years());
    }

    #[test]
    fn migration_timescale_decreases_with_planet_mass() {
        // In Type I migration, more massive planets migrate FASTER (shorter timescales)
        // This is because the torque scales as M_p² while angular momentum scales as M_p
        let stellar_mass = Mass::from_solar_masses(1.0);
        let radius = Length::from_au(5.0);
        let sigma = SurfaceDensity::from_grams_per_cm2(100.0);
        let h_over_r = 0.05;

        let tau_small = type_i_migration_timescale(
            Mass::from_earth_masses(1.0),
            stellar_mass,
            radius,
            sigma,
            h_over_r,
        );

        let tau_large = type_i_migration_timescale(
            Mass::from_earth_masses(10.0),
            stellar_mass,
            radius,
            sigma,
            h_over_r,
        );

        // More massive planet → shorter timescale (faster migration)
        assert!(tau_large.to_years() < tau_small.to_years());
    }

    #[test]
    fn migration_velocity_is_negative() {
        let r = Length::from_au(5.0);
        let tau = Time::from_years(100_000.0);

        let v = migration_velocity(r, tau);

        // Should be negative (inward migration)
        assert!(v < 0.0);
    }

    #[test]
    fn migration_distance_scales_with_timestep() {
        let r = Length::from_au(5.0);
        let tau = Time::from_years(100_000.0);

        let dr_small = migration_distance(r, tau, Time::from_years(1_000.0));
        let dr_large = migration_distance(r, tau, Time::from_years(10_000.0));

        // Larger timestep → larger migration distance
        assert!(dr_large.abs() > dr_small.abs());

        // But both should be negative (inward)
        assert!(dr_small < 0.0);
        assert!(dr_large < 0.0);
    }

    #[test]
    fn migration_distance_capped() {
        let r = Length::from_au(5.0);
        let tau = Time::from_years(1_000.0); // Very short timescale

        // Even with large timestep, distance capped at 0.5 AU
        let dr = migration_distance(r, tau, Time::from_years(10_000.0));

        assert!(dr >= -0.5, "Migration distance should be capped at 0.5 AU");
    }

    #[test]
    fn migration_timescale_decreases_with_surface_density() {
        let planet_mass = Mass::from_earth_masses(5.0);
        let stellar_mass = Mass::from_solar_masses(1.0);
        let radius = Length::from_au(5.0);
        let h_over_r = 0.05;

        let tau_low_sigma = type_i_migration_timescale(
            planet_mass,
            stellar_mass,
            radius,
            SurfaceDensity::from_grams_per_cm2(50.0),
            h_over_r,
        );

        let tau_high_sigma = type_i_migration_timescale(
            planet_mass,
            stellar_mass,
            radius,
            SurfaceDensity::from_grams_per_cm2(500.0),
            h_over_r,
        );

        // Denser disk → faster migration → shorter timescale
        assert!(tau_high_sigma.to_years() < tau_low_sigma.to_years());
    }
}
