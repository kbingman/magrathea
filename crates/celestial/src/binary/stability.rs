//! Planetary stability calculations for binary star systems.
//!
//! Provides functions for calculating stable orbital zones around binary stars,
//! based on empirical fits from numerical integrations.
//!
//! # Orbital Types
//!
//! - **S-type (circumstellar)**: Planet orbits one star, companion is external
//!   - Stable within ~1/3 to 1/2 of binary separation
//!   - Tighter constraint for eccentric binaries
//!
//! - **P-type (circumbinary)**: Planet orbits both stars
//!   - Stable beyond ~2-4× binary separation
//!   - Larger exclusion zone for eccentric binaries
//!
//! # References
//! - Holman & Wiegert (1999) - "Dynamical Stability of Planets in Binary Systems"
//! - Quarles et al. (2020) - "Stability Limits of Circumbinary Planets"

use super::orbital_params::{BinaryOrbitType, OrbitalParameters};
use stellar::StellarObject;
use units::Length;

/// Calculate the maximum stable orbital radius for an S-type orbit
///
/// Returns the critical semi-major axis beyond which planetary orbits
/// become unstable. Based on Holman & Wiegert (1999) empirical fits.
///
/// # Formula
/// For planets orbiting one star:
/// ```text
/// a_crit = a_bin * (0.464 - 0.38*μ - 0.631*e + 0.586*μ*e + 0.15*e² - 0.198*μ*e²)
/// ```
/// where:
/// - a_bin = binary semi-major axis
/// - μ = M₂/(M₁+M₂) (mass parameter)
/// - e = binary eccentricity
///
/// # Arguments
/// * `binary_orbit` - The binary orbital parameters
/// * `primary_mass` - Mass of the primary star in solar masses
/// * `secondary_mass` - Mass of the secondary star in solar masses
/// * `orbit_type` - Whether planet orbits primary or secondary
///
/// # Returns
/// Maximum stable semi-major axis in AU
///
/// # References
/// - Holman & Wiegert (1999) ApJ 117, 621
pub fn s_type_stability_limit(
    binary_orbit: &OrbitalParameters,
    primary_mass: f64,
    secondary_mass: f64,
    orbit_type: BinaryOrbitType,
) -> Length {
    let a_bin = binary_orbit.semi_major_axis.to_au();
    let e = binary_orbit.eccentricity;

    // Mass parameter μ = M₂/(M₁+M₂)
    let mu = match orbit_type {
        BinaryOrbitType::STypePrimary => secondary_mass / (primary_mass + secondary_mass),
        BinaryOrbitType::STypeSecondary => primary_mass / (primary_mass + secondary_mass), // Swap roles
        BinaryOrbitType::PType => {
            panic!("Use p_type_stability_limit for P-type orbits")
        }
    };

    // Holman & Wiegert (1999) empirical fit
    let factor =
        0.464 - 0.38 * mu - 0.631 * e + 0.586 * mu * e + 0.15 * e.powi(2) - 0.198 * mu * e.powi(2);

    let a_crit = a_bin * factor;

    Length::from_au(a_crit)
}

/// Calculate the minimum stable orbital radius for a P-type orbit
///
/// Returns the critical semi-major axis below which circumbinary planetary
/// orbits become unstable. Based on Holman & Wiegert (1999) empirical fits.
///
/// # Formula
/// ```text
/// a_crit = a_bin * (1.60 + 4.12*μ + 5.10*e - 4.27*μ*e - 2.22*e² + 4.61*μ²*e²)
/// ```
///
/// # Arguments
/// * `binary_orbit` - The binary orbital parameters
/// * `primary_mass` - Mass of the primary star in solar masses
/// * `secondary_mass` - Mass of the secondary star in solar masses
///
/// # Returns
/// Minimum stable semi-major axis in AU
///
/// # References
/// - Holman & Wiegert (1999) ApJ 117, 621
pub fn p_type_stability_limit(
    binary_orbit: &OrbitalParameters,
    primary_mass: f64,
    secondary_mass: f64,
) -> Length {
    let a_bin = binary_orbit.semi_major_axis.to_au();
    let e = binary_orbit.eccentricity;

    // Mass parameter μ = M₂/(M₁+M₂)
    let mu = secondary_mass / (primary_mass + secondary_mass);

    // Holman & Wiegert (1999) empirical fit
    let factor = 1.60 + 4.12 * mu + 5.10 * e - 4.27 * mu * e - 2.22 * e.powi(2)
        + 4.61 * mu.powi(2) * e.powi(2);

    let a_crit = a_bin * factor;

    Length::from_au(a_crit)
}

/// Calculate the stable orbital range for planets in a binary system
///
/// Returns the range of stable planetary orbits based on orbit type.
///
/// # Arguments
/// * `binary_orbit` - The binary orbital parameters
/// * `stars` - The stellar objects (must have exactly 2 stars)
/// * `orbit_type` - Type of planetary orbit
///
/// # Returns
/// A tuple `(min_au, max_au)` representing the stable orbital range.
/// For S-type: `(stellar_radius * 3, stability_limit)`
/// For P-type: `(stability_limit * 1.2, binary_separation * 100)`
///
/// # Example
/// ```
/// use celestial::binary::{OrbitalParameters, BinaryOrbitType, stable_orbital_range};
/// use stellar::{StellarObject, main_sequence_star};
/// use units::{Length, Time};
///
/// let orbit = OrbitalParameters {
///     semi_major_axis: Length::from_au(10.0),
///     eccentricity: 0.1,
///     inclination: 0.0,
///     longitude_of_ascending_node: 0.0,
///     argument_of_periapsis: 0.0,
///     mean_anomaly: 0.0,
///     period: Time::from_years(10.0),
/// };
///
/// let stars = vec![
///     StellarObject::MainSequence(main_sequence_star(1.0, 0.0, 5e9)),
///     StellarObject::MainSequence(main_sequence_star(0.5, 0.0, 5e9)),
/// ];
///
/// let (min_au, max_au) = stable_orbital_range(&orbit, &stars, BinaryOrbitType::STypePrimary);
/// println!("Planets can orbit between {} and {} AU", min_au, max_au);
/// ```
pub fn stable_orbital_range(
    binary_orbit: &OrbitalParameters,
    stars: &[StellarObject],
    orbit_type: BinaryOrbitType,
) -> (f64, f64) {
    assert_eq!(
        stars.len(),
        2,
        "stable_orbital_range requires exactly 2 stars"
    );

    let primary_mass = stars[0].mass().to_solar_masses();
    let secondary_mass = stars[1].mass().to_solar_masses();

    match orbit_type {
        BinaryOrbitType::STypePrimary => {
            let stellar_radius = stars[0].radius().to_solar_radii() * 0.00465; // Convert to AU
            let max_stable =
                s_type_stability_limit(binary_orbit, primary_mass, secondary_mass, orbit_type)
                    .to_au();

            // Conservative minimum: 3× stellar radius to avoid tidal disruption
            let min_orbit = (stellar_radius * 3.0).max(0.01);

            (min_orbit, max_stable)
        }
        BinaryOrbitType::STypeSecondary => {
            let stellar_radius = stars[1].radius().to_solar_radii() * 0.00465; // Convert to AU
            let max_stable =
                s_type_stability_limit(binary_orbit, primary_mass, secondary_mass, orbit_type)
                    .to_au();

            let min_orbit = (stellar_radius * 3.0).max(0.01);

            (min_orbit, max_stable)
        }
        BinaryOrbitType::PType => {
            let min_stable =
                p_type_stability_limit(binary_orbit, primary_mass, secondary_mass).to_au();

            // Conservative factor: add 20% margin for long-term stability
            let min_orbit = min_stable * 1.2;

            // Max distance limited by ~100× binary separation or 1000 AU
            let max_orbit = (binary_orbit.semi_major_axis.to_au() * 100.0).min(1000.0);

            (min_orbit, max_orbit)
        }
    }
}

/// Calculate the habitable zone in a binary system
///
/// Returns the conservative habitable zone boundaries accounting for
/// both stellar luminosities and stability constraints.
///
/// # Arguments
/// * `binary_orbit` - The binary orbital parameters
/// * `stars` - The stellar objects (must have exactly 2 stars)
/// * `orbit_type` - Type of planetary orbit
///
/// # Returns
/// A tuple `(inner_hz_au, outer_hz_au)` or `None` if no habitable zone exists
///
/// # Algorithm
/// 1. Calculate habitable zone based on total/individual luminosity
/// 2. Check if HZ overlaps with stable orbital range
/// 3. Return intersection if it exists
pub fn habitable_zone(
    binary_orbit: &OrbitalParameters,
    stars: &[StellarObject],
    orbit_type: BinaryOrbitType,
) -> Option<(f64, f64)> {
    assert_eq!(stars.len(), 2, "habitable_zone requires exactly 2 stars");

    // Calculate HZ boundaries based on orbit type
    let (hz_inner, hz_outer) = match orbit_type {
        BinaryOrbitType::STypePrimary => {
            let lum = stars[0].luminosity();
            calculate_hz_from_luminosity(lum)
        }
        BinaryOrbitType::STypeSecondary => {
            let lum = stars[1].luminosity();
            calculate_hz_from_luminosity(lum)
        }
        BinaryOrbitType::PType => {
            // For circumbinary, use total system luminosity
            let total_lum = stars[0].luminosity() + stars[1].luminosity();
            calculate_hz_from_luminosity(total_lum)
        }
    };

    // Get stable orbital range
    let (stable_min, stable_max) = stable_orbital_range(binary_orbit, stars, orbit_type);

    // Find intersection of HZ and stable range
    let effective_inner = hz_inner.max(stable_min);
    let effective_outer = hz_outer.min(stable_max);

    if effective_inner < effective_outer {
        Some((effective_inner, effective_outer))
    } else {
        None // No overlap - no habitable zone
    }
}

/// Calculate habitable zone from stellar luminosity
///
/// Uses simple scaling: HZ distance scales as sqrt(L/L☉)
/// Conservative limits:
/// - Inner edge: 0.75 AU at 1 L☉ (runaway greenhouse)
/// - Outer edge: 1.77 AU at 1 L☉ (maximum greenhouse)
fn calculate_hz_from_luminosity(luminosity: f64) -> (f64, f64) {
    let sqrt_l = luminosity.sqrt();

    let inner = 0.75 * sqrt_l; // Conservative inner edge
    let outer = 1.77 * sqrt_l; // Conservative outer edge

    (inner, outer)
}

#[cfg(test)]
mod tests {
    use super::*;
    use stellar::main_sequence_star;
    use units::Time;

    fn create_test_binary_orbit(semi_major_axis_au: f64, eccentricity: f64) -> OrbitalParameters {
        OrbitalParameters {
            semi_major_axis: Length::from_au(semi_major_axis_au),
            eccentricity,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0,
            period: Time::from_years(1.0), // Simplified for tests
        }
    }

    fn create_test_stars() -> Vec<StellarObject> {
        vec![
            StellarObject::MainSequence(main_sequence_star(1.0, 0.0, 5.0e9)),
            StellarObject::MainSequence(main_sequence_star(0.5, 0.0, 5.0e9)),
        ]
    }

    #[test]
    fn test_s_type_stability_limit() {
        let orbit = create_test_binary_orbit(10.0, 0.1);
        let stars = create_test_stars();

        let limit = s_type_stability_limit(
            &orbit,
            stars[0].mass().to_solar_masses(),
            stars[1].mass().to_solar_masses(),
            BinaryOrbitType::STypePrimary,
        );

        // Should be positive and less than binary separation
        assert!(limit.to_au() > 0.0);
        assert!(limit.to_au() < 10.0);
    }

    #[test]
    fn test_p_type_stability_limit() {
        let orbit = create_test_binary_orbit(1.0, 0.1);
        let stars = create_test_stars();

        let limit = p_type_stability_limit(
            &orbit,
            stars[0].mass().to_solar_masses(),
            stars[1].mass().to_solar_masses(),
        );

        // Inner exclusion zone should be larger than binary orbit
        assert!(limit.to_au() > 1.0);
    }

    #[test]
    fn test_stable_orbital_range_s_type() {
        let orbit = create_test_binary_orbit(20.0, 0.2);
        let stars = create_test_stars();

        let (min, max) = stable_orbital_range(&orbit, &stars, BinaryOrbitType::STypePrimary);

        assert!(min > 0.0);
        assert!(max > min);
        assert!(max < 20.0); // Should be less than binary separation
    }

    #[test]
    fn test_stable_orbital_range_p_type() {
        let orbit = create_test_binary_orbit(1.0, 0.1);
        let stars = create_test_stars();

        let (min, max) = stable_orbital_range(&orbit, &stars, BinaryOrbitType::PType);

        assert!(min > 1.0); // Should be beyond binary separation
        assert!(max > min);
    }

    #[test]
    fn test_habitable_zone_s_type() {
        let orbit = create_test_binary_orbit(50.0, 0.1);
        let stars = create_test_stars();

        let hz = habitable_zone(&orbit, &stars, BinaryOrbitType::STypePrimary);

        // Wide binary should have habitable zone around primary
        assert!(hz.is_some());

        if let Some((inner, outer)) = hz {
            assert!(inner > 0.0);
            assert!(outer > inner);
        }
    }

    #[test]
    fn test_habitable_zone_p_type() {
        let orbit = create_test_binary_orbit(0.5, 0.0);
        let stars = create_test_stars();

        let hz = habitable_zone(&orbit, &stars, BinaryOrbitType::PType);

        // Close binary might have circumbinary HZ
        if let Some((inner, outer)) = hz {
            assert!(outer > inner);
            assert!(inner > 0.5); // Should be beyond binary separation
        }
    }

    #[test]
    fn test_hz_scaling_with_luminosity() {
        let (inner1, outer1) = calculate_hz_from_luminosity(1.0);
        let (inner4, outer4) = calculate_hz_from_luminosity(4.0);

        // HZ should scale as sqrt(L)
        assert!((inner4 / inner1 - 2.0).abs() < 0.01);
        assert!((outer4 / outer1 - 2.0).abs() < 0.01);
    }
}
