//! Binary star system generation from stellar populations.
//!
//! Provides functions for generating binary star systems with physically realistic
//! mass ratios, orbital separations, and binary fractions based on observational data.
//!
//! # References
//! - Raghavan et al. (2010) - "A Survey of Stellar Families: Multiplicity of Solar-type Stars"
//! - Duchêne & Kraus (2013) - "Stellar Multiplicity"
//! - Moe & Di Stefano (2017) - "Mind Your Ps and Qs: The Interplay of Stellar and Planetary Multiplicity"

use rand::Rng;
use rand_chacha::ChaChaRng;
use star_system::binary::{BinaryConfiguration, BinaryOrbitType, OrbitalParameters};
use stellar::{StellarObject, sampling::sample_age_from_population, stellar_object};
use units::{Length, Time};

/// Determine if a system should be binary based on primary star mass
///
/// Binary fractions increase with stellar mass, ranging from ~25% for M-dwarfs
/// to ~80% for O-type stars.
///
/// # Mass-Dependent Binary Fractions
/// - M ≥ 16 M☉ (O stars): 80% binary
/// - M ≥ 3 M☉ (B/A stars): 60% binary
/// - M ≥ 0.8 M☉ (F/G stars): 44% binary
/// - M ≥ 0.45 M☉ (K stars): 35% binary
/// - M < 0.45 M☉ (M stars): 25% binary
///
/// # Arguments
/// * `primary_mass` - Mass of the primary star in solar masses
///
/// # Returns
/// Probability (0.0-1.0) that the system is binary
///
/// # References
/// - Raghavan et al. (2010) for FGK stars
/// - Duchêne & Kraus (2013) for M dwarfs
/// - Mason et al. (2009) for massive stars
pub fn binary_fraction(primary_mass: f64) -> f64 {
    match primary_mass {
        m if m >= 16.0 => 0.80, // O stars
        m if m >= 3.0 => 0.60,  // B & A stars
        m if m >= 0.8 => 0.44,  // F & G stars
        m if m >= 0.45 => 0.35, // K stars
        _ => 0.25,              // M stars
    }
}

/// Generate mass ratio q = M₂/M₁ based on observational distributions
///
/// The mass ratio distribution depends on primary mass:
/// - Massive stars (> 3 M☉): Roughly flat distribution, q ∈ [0.2, 1.0]
/// - Solar-type stars: Distribution favoring lower mass ratios
///
/// # Arguments
/// * `rng` - Random number generator
/// * `primary_mass` - Mass of the primary star in solar masses
///
/// # Returns
/// Mass ratio q = M₂/M₁, where 0.2 ≤ q ≤ 1.0
///
/// # References
/// - Moe & Di Stefano (2017) - mass ratio distributions
pub fn sample_mass_ratio(rng: &mut ChaChaRng, primary_mass: f64) -> f64 {
    match primary_mass {
        m if m > 3.0 => {
            // Massive stars: roughly flat distribution
            rng.random_range(0.2..1.0)
        }
        _ => {
            // Solar-type and lower: favor lower mass ratios
            let x: f64 = rng.random();
            0.2 + (x * x * 0.8) // Biased towards lower mass ratios
        }
    }
}

/// Generate orbital separation for a binary system
///
/// Returns semi-major axis in AU. Distribution is log-normal with mass-dependent
/// center:
/// - Massive stars (> 3 M☉): Median ~100 AU
/// - Solar-type (0.8-3 M☉): Median ~30 AU
/// - Low-mass (< 0.8 M☉): Median ~10 AU
///
/// # Arguments
/// * `rng` - Random number generator
/// * `primary_mass` - Mass of the primary star in solar masses
///
/// # Returns
/// Orbital semi-major axis in AU
///
/// # References
/// - Raghavan et al. (2010) - separation distribution for solar-type stars
pub fn sample_separation(rng: &mut ChaChaRng, primary_mass: f64) -> f64 {
    // Log-normal distribution with mass-dependent center
    let log_mean = match primary_mass {
        m if m > 3.0 => 2.0, // ~100 AU for massive stars
        m if m > 0.8 => 1.5, // ~30 AU for solar-type
        _ => 1.0,            // ~10 AU for M-dwarfs
    };

    let log_sep = rng.random_range(log_mean - 1.0..log_mean + 1.0);
    10.0_f64.powf(log_sep)
}

/// Generate orbital eccentricity for a binary system
///
/// Distribution depends on orbital period:
/// - Short period (< 10 days): Tidally circularized, e ~ 0
/// - Medium period (10-100 days): Low eccentricity, e < 0.3
/// - Long period (> 100 days): Thermal distribution, f(e) ∝ 2e
///
/// # Arguments
/// * `rng` - Random number generator
/// * `period_years` - Orbital period in years
///
/// # Returns
/// Orbital eccentricity (0.0-0.95)
pub fn sample_binary_eccentricity(rng: &mut ChaChaRng, period_years: f64) -> f64 {
    let period_days = period_years * 365.25;

    if period_days < 10.0 {
        // Tidally circularized
        rng.random_range(0.0..0.05)
    } else if period_days < 100.0 {
        // Partially circularized
        rng.random_range(0.0..0.3)
    } else {
        // Thermal distribution: f(e) ∝ 2e
        let u: f64 = rng.random();
        u.sqrt().min(0.95)
    }
}

/// Calculate orbital period from semi-major axis using Kepler's third law
///
/// P² = (4π²/G(M₁+M₂)) × a³
///
/// # Arguments
/// * `semi_major_axis_au` - Orbital separation in AU
/// * `total_mass_solar` - Combined mass of both stars in solar masses
///
/// # Returns
/// Orbital period in years
fn calculate_period(semi_major_axis_au: f64, total_mass_solar: f64) -> f64 {
    // Using Kepler's third law in convenient units:
    // P (years) = sqrt(a³ / M_total)
    // where a is in AU and M is in solar masses
    (semi_major_axis_au.powi(3) / total_mass_solar).sqrt()
}

/// Generate a binary companion for an existing primary star
///
/// Given a primary star, generates a companion with appropriate mass ratio,
/// orbital parameters, and age from the same stellar population.
///
/// # Arguments
/// * `rng` - Random number generator
/// * `primary` - The primary stellar object
///
/// # Returns
/// A tuple containing:
/// - Secondary `StellarObject`
/// - `BinaryConfiguration` with orbital parameters
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar::stellar_object;
/// use planetary_generator::binary::generation::generate_companion;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let primary = stellar_object(&mut rng, 1.0, 5e9, 0.0);
/// let (secondary, config) = generate_companion(&mut rng, &primary);
/// ```
pub fn generate_companion(
    rng: &mut ChaChaRng,
    primary: &StellarObject,
) -> (StellarObject, BinaryConfiguration) {
    let primary_mass = primary.mass().to_solar_masses();
    let metallicity = primary.metallicity();

    // Get primary age in years, or use a default 5 Gyr if not available
    let primary_age_years = primary.age().map(|t| t.to_years()).unwrap_or(5e9);

    // Sample binary parameters
    let mass_ratio = sample_mass_ratio(rng, primary_mass);
    let secondary_mass = primary_mass * mass_ratio;
    let total_mass = primary_mass + secondary_mass;

    // Sample separation and calculate period
    let semi_major_axis_au = sample_separation(rng, primary_mass);
    let period_years = calculate_period(semi_major_axis_au, total_mass);

    // Sample eccentricity
    let eccentricity = sample_binary_eccentricity(rng, period_years);

    // Sample random orbital angles (isotropic distribution)
    let inclination = sample_inclination(rng);
    let longitude_of_ascending_node = rng.random_range(0.0..360.0);
    let argument_of_periapsis = rng.random_range(0.0..360.0);
    let mean_anomaly = rng.random_range(0.0..360.0);

    // Generate secondary with same age as primary (coeval formation)
    let secondary = stellar_object(rng, secondary_mass, primary_age_years, metallicity);

    // Determine default orbit type (S-type around primary for planet hosting)
    let orbit_type = BinaryOrbitType::STypePrimary;

    let orbital_params = OrbitalParameters {
        semi_major_axis: Length::from_au(semi_major_axis_au),
        eccentricity,
        inclination,
        longitude_of_ascending_node,
        argument_of_periapsis,
        mean_anomaly,
        period: Time::from_years(period_years),
    };

    let config = BinaryConfiguration {
        orbit_type,
        orbital_params,
    };

    (secondary, config)
}

/// Generate a complete binary star system
///
/// Creates two stars with orbital parameters sampled from observational distributions.
/// The system includes full Keplerian orbital elements and binary configuration.
///
/// # Arguments
/// * `rng` - Random number generator
/// * `primary_mass` - Mass of the primary star in solar masses
/// * `metallicity` - Stellar metallicity [Fe/H] in dex
/// * `mean_age_gyr` - Mean population age in billions of years
/// * `age_spread_gyr` - Age spread (standard deviation) in billions of years
///
/// # Returns
/// A tuple containing:
/// - Primary `StellarObject`
/// - Secondary `StellarObject`
/// - `BinaryConfiguration` with orbital parameters
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use planetary_generator::binary::generation::generate_binary_system;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let (primary, secondary, config) = generate_binary_system(
///     &mut rng,
///     1.0,   // Solar-mass primary
///     0.0,   // Solar metallicity
///     5.0,   // 5 Gyr mean age
///     2.0    // 2 Gyr age spread
/// );
/// ```
pub fn generate_binary_system(
    rng: &mut ChaChaRng,
    primary_mass: f64,
    metallicity: f64,
    mean_age_gyr: f64,
    age_spread_gyr: f64,
) -> (StellarObject, StellarObject, BinaryConfiguration) {
    // Sample binary parameters
    let mass_ratio = sample_mass_ratio(rng, primary_mass);
    let secondary_mass = primary_mass * mass_ratio;
    let total_mass = primary_mass + secondary_mass;

    // Sample separation and calculate period
    let semi_major_axis_au = sample_separation(rng, primary_mass);
    let period_years = calculate_period(semi_major_axis_au, total_mass);

    // Sample eccentricity
    let eccentricity = sample_binary_eccentricity(rng, period_years);

    // Sample random orbital angles (isotropic distribution)
    let inclination = sample_inclination(rng);
    let longitude_of_ascending_node = rng.random_range(0.0..360.0);
    let argument_of_periapsis = rng.random_range(0.0..360.0);
    let mean_anomaly = rng.random_range(0.0..360.0);

    // Sample ages for both stars (from same population)
    let primary_age_gyr =
        sample_age_from_population(rng, primary_mass, mean_age_gyr, age_spread_gyr);
    let secondary_age_gyr =
        sample_age_from_population(rng, secondary_mass, mean_age_gyr, age_spread_gyr);

    // Generate stellar objects
    let primary = stellar_object(rng, primary_mass, primary_age_gyr * 1e9, metallicity);
    let secondary = stellar_object(rng, secondary_mass, secondary_age_gyr * 1e9, metallicity);

    // Determine default orbit type (S-type around primary for planet hosting)
    let orbit_type = BinaryOrbitType::STypePrimary;

    let orbital_params = OrbitalParameters {
        semi_major_axis: Length::from_au(semi_major_axis_au),
        eccentricity,
        inclination,
        longitude_of_ascending_node,
        argument_of_periapsis,
        mean_anomaly,
        period: Time::from_years(period_years),
    };

    let config = BinaryConfiguration {
        orbit_type,
        orbital_params,
    };

    (primary, secondary, config)
}

/// Sample orbital inclination from isotropic distribution
///
/// For randomly oriented orbits, P(i) ∝ sin(i), which gives
/// a uniform distribution in cos(i).
///
/// # Arguments
/// * `rng` - Random number generator
///
/// # Returns
/// Inclination in degrees (0-180°)
fn sample_inclination(rng: &mut ChaChaRng) -> f64 {
    let cos_i: f64 = rng.random_range(-1.0..1.0);
    cos_i.acos().to_degrees()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn test_binary_fraction_increases_with_mass() {
        assert!(binary_fraction(0.3) < binary_fraction(1.0));
        assert!(binary_fraction(1.0) < binary_fraction(5.0));
        assert!(binary_fraction(5.0) < binary_fraction(20.0));
    }

    #[test]
    fn test_mass_ratio_range() {
        let mut rng = ChaChaRng::seed_from_u64(42);

        for _ in 0..100 {
            let q = sample_mass_ratio(&mut rng, 1.0);
            assert!(q >= 0.2);
            assert!(q <= 1.0);
        }
    }

    #[test]
    fn test_separation_positive() {
        let mut rng = ChaChaRng::seed_from_u64(42);

        for mass in [0.3, 1.0, 5.0, 20.0] {
            let sep = sample_separation(&mut rng, mass);
            assert!(sep > 0.0);
            assert!(sep < 10000.0); // Reasonable upper limit
        }
    }

    #[test]
    fn test_eccentricity_depends_on_period() {
        let mut rng = ChaChaRng::seed_from_u64(42);

        // Short period should give low eccentricity
        let e_short = sample_binary_eccentricity(&mut rng, 0.01); // ~3.6 days
        assert!(e_short < 0.1);

        // Long period can give higher eccentricity
        let e_long = sample_binary_eccentricity(&mut rng, 100.0);
        // Just check it's in valid range
        assert!(e_long >= 0.0);
        assert!(e_long < 1.0);
    }

    #[test]
    fn test_generate_binary_system() {
        let mut rng = ChaChaRng::seed_from_u64(42);

        let (primary, secondary, config) = generate_binary_system(&mut rng, 1.0, 0.0, 5.0, 2.0);

        // Check masses are reasonable
        assert!(primary.mass().to_solar_masses() > 0.0);
        assert!(secondary.mass().to_solar_masses() > 0.0);
        assert!(secondary.mass().to_solar_masses() <= primary.mass().to_solar_masses());

        // Check orbital parameters are valid
        assert!(config.orbital_params.semi_major_axis.to_au() > 0.0);
        assert!(config.orbital_params.eccentricity >= 0.0);
        assert!(config.orbital_params.eccentricity < 1.0);
        assert!(config.orbital_params.period.to_years() > 0.0);
    }

    #[test]
    fn test_keplers_third_law() {
        // Test with Earth-Sun system
        let period = calculate_period(1.0, 1.0);
        assert!((period - 1.0).abs() < 0.01);

        // Binary with 2 M☉ total, same separation
        let period2 = calculate_period(1.0, 2.0);
        // Period should be shorter for more massive system
        assert!(period2 < period);
    }
}
