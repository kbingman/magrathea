use std::f64::consts::PI;

use rand::Rng;
use rand_chacha::ChaChaRng;

/// Sample from a Gaussian (normal) distribution using Box-Muller transform
///
/// # Arguments
/// * `rng` - Random number generator
/// * `mean` - Mean of the distribution
/// * `std_dev` - Standard deviation
///
/// # Returns
/// A sample from the normal distribution N(mean, std_dev²)
pub fn sample_gaussian(rng: &mut ChaChaRng, mean: f64, std_dev: f64) -> f64 {
    let u1: f64 = rng.random();
    let u2: f64 = rng.random();
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
    mean + std_dev * z
}

/// Sample from a power-law distribution
///
/// Samples from p(x) ∝ x^α between x_min and x_max using inverse transform sampling.
///
/// # Arguments
/// * `x_min` - Minimum value
/// * `x_max` - Maximum value
/// * `alpha` - Power-law exponent (e.g., -2.3 for Salpeter IMF)
/// * `rng` - Random number generator
///
/// # Returns
/// A sample from the power-law distribution
pub fn sample_power_law(x_min: f64, x_max: f64, alpha: f64, rng: &mut ChaChaRng) -> f64 {
    let u: f64 = rng.random();
    let alpha1 = alpha + 1.0;
    (u * (x_max.powf(alpha1) - x_min.powf(alpha1)) + x_min.powf(alpha1)).powf(1.0 / alpha1)
}

/// Sample metallicity from local galactic distribution
///
/// Returns [Fe/H] in dex, centered on solar (0.0) with σ ≈ 0.2 dex.
/// Values are clamped to the range [-0.5, 0.4] to avoid extreme outliers.
///
/// # Arguments
/// * `rng` - Random number generator
///
/// # Returns
/// Metallicity [Fe/H] in dex
pub fn sample_metallicity(rng: &mut ChaChaRng) -> f64 {
    sample_gaussian(rng, 0.0, 0.2).clamp(-0.5, 0.4)
}

/// Sample stellar mass from the Kroupa (2001) Initial Mass Function
///
/// The Kroupa IMF is a broken power law:
/// - 0.08 ≤ M < 0.5 M☉: α = -1.3 (brown dwarf/M dwarf transition)
/// - 0.5 ≤ M < 1.0 M☉: α = -2.3
/// - M ≥ 1.0 M☉: α = -2.3
///
/// # Arguments
/// * `rng` - Random number generator
/// * `max_mass` - Upper mass limit in solar masses (default behavior uses 150 M☉)
///
/// # Returns
/// Stellar mass in solar masses
pub fn sample_mass_kroupa(rng: &mut ChaChaRng, max_mass: f64) -> f64 {
    // Segment weights from integrating each segment of the IMF
    // These are approximate and assume the full mass range
    let segment_weights = [0.80, 0.15, 0.05];
    let rand: f64 = rng.random();

    if rand < segment_weights[0] {
        // Low mass: 0.08-0.5 M☉, slope -1.3
        sample_power_law(0.08, max_mass.min(0.5), -1.3, rng)
    } else if rand < segment_weights[0] + segment_weights[1] {
        // Intermediate: 0.5-1.0 M☉, slope -2.3
        sample_power_law(0.5, max_mass.min(1.0), -2.3, rng)
    } else {
        // High mass: 1.0-max_mass M☉, slope -2.3
        sample_power_law(1.0, max_mass, -2.3, rng)
    }
}

/// Sample stellar age from a population with given mean age and spread
///
/// Ages are drawn from a normal distribution, then clamped to physical limits:
/// - Minimum: 0.1 Gyr (very young populations)
/// - Maximum: min(stellar_lifetime, universe_age)
///
/// This allows modeling different stellar populations:
/// - Young populations (mean ~0.5 Gyr, spread ~0.3 Gyr): Star-forming regions
/// - Middle-aged (mean ~5 Gyr, spread ~2 Gyr): Thin disk like the Sun
/// - Old populations (mean ~10 Gyr, spread ~2 Gyr): Thick disk, halo
///
/// # Arguments
/// * `rng` - Random number generator
/// * `stellar_mass` - Mass of the star in solar masses (used to cap age at lifetime)
/// * `mean_age_gyr` - Mean population age in billions of years
/// * `age_spread_gyr` - Standard deviation of age distribution in billions of years
///
/// # Returns
/// Stellar age in billions of years
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar::sampling::sample_age_from_population;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
///
/// // Young star-forming region
/// let young_age = sample_age_from_population(&mut rng, 1.0, 0.5, 0.3);
///
/// // Solar neighborhood (thin disk)
/// let thin_disk_age = sample_age_from_population(&mut rng, 1.0, 5.0, 2.0);
/// ```
pub fn sample_age_from_population(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    mean_age_gyr: f64,
    age_spread_gyr: f64,
) -> f64 {
    const UNIVERSE_AGE_GYR: f64 = 13.8;
    const MIN_AGE_GYR: f64 = 0.1;

    // Calculate stellar lifetime
    let lifetime_gyr = crate::generation::estimate_lifetime(stellar_mass);

    // Maximum age is the minimum of stellar lifetime and universe age
    let max_age = lifetime_gyr.min(UNIVERSE_AGE_GYR);

    // Sample from normal distribution
    let age = sample_gaussian(rng, mean_age_gyr, age_spread_gyr);

    // Clamp to physical limits
    age.clamp(MIN_AGE_GYR, max_age)
}
