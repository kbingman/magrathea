//! Planet occurrence rate sampling
//!
//! Kepler/TESS-derived occurrence rates for statistical planet generation.

use rand::Rng;
use rand_chacha::ChaChaRng;

use planetary::planet_class::PlanetClass;

/// Sample planet mass using occurrence-weighted distribution
///
/// This function is used for inner system planets (within snow line) where
/// planet masses are sampled from Kepler-derived occurrence rates.
///
/// For outer system planets (giants, ice giants), use the dedicated functions
/// in generation.rs which properly account for stellar mass limits.
///
/// Note: Jupiter-mass planets (>160 M⊕) are NOT sampled here. Gas giants in the
/// inner system (Hot Jupiters) form via migration and are handled separately
/// in generate_giant_system.
pub fn sample_planet_mass(rng: &mut ChaChaRng, stellar_metallicity: f64) -> f64 {
    // Metallicity affects larger planet probability: P(sub-Saturn) ∝ 10^(2×[Fe/H])
    let metallicity_boost = 10.0_f64.powf(2.0 * stellar_metallicity);

    // Inner system mass distribution (Kepler-derived)
    // No Jupiter bin - those form via core accretion beyond snow line and migrate
    let base_weights = [
        0.15,                     // Sub-Earth (0.01-0.5 M⊕)
        0.30,                     // Earth (0.5-2 M⊕)
        0.30,                     // Super-Earth/Mini-Neptune (2-10 M⊕)
        0.20,                     // Neptune (10-50 M⊕)
        0.05 * metallicity_boost, // Sub-Saturn (50-160 M⊕) - rare, metallicity-dependent
    ];

    let total: f64 = base_weights.iter().sum();
    let roll: f64 = rng.random::<f64>() * total;
    let mut cumulative = 0.0;

    for (i, &weight) in base_weights.iter().enumerate() {
        cumulative += weight;
        if roll < cumulative {
            return sample_mass_in_bin(rng, i);
        }
    }

    1.0 // Fallback
}

fn sample_mass_in_bin(rng: &mut ChaChaRng, bin: usize) -> f64 {
    let (min, max): (f64, f64) = match bin {
        0 => (0.01, 0.5),
        1 => (0.5, 2.0),
        2 => (2.0, 10.0),
        3 => (10.0, 50.0),
        4 => (50.0, 160.0),
        5 => (160.0, 1000.0),
        _ => (0.5, 2.0),
    };

    // Log-uniform within bin
    let log_min = min.ln();
    let log_max = max.ln();
    (log_min + rng.random::<f64>() * (log_max - log_min)).exp()
}

/// Sample orbital period in days
pub fn sample_orbital_period(rng: &mut ChaChaRng, planet_class: &PlanetClass) -> f64 {
    match planet_class {
        PlanetClass::Compact | PlanetClass::Transitional => {
            let weights = [0.15, 0.30, 0.30, 0.15, 0.10];
            let bins = [
                (1.0, 10.0),
                (10.0, 30.0),
                (30.0, 100.0),
                (100.0, 300.0),
                (300.0, 1000.0),
            ];
            sample_from_bins(rng, &weights, &bins)
        }
        PlanetClass::Volatile => {
            let weights = [0.05, 0.15, 0.30, 0.30, 0.20];
            let bins = [
                (10.0, 30.0),
                (30.0, 100.0),
                (100.0, 365.0),
                (365.0, 1000.0),
                (1000.0, 5000.0),
            ];
            sample_from_bins(rng, &weights, &bins)
        }
        PlanetClass::Giant => {
            if rng.random::<f64>() < 0.15 {
                1.0 + rng.random::<f64>() * 9.0 // Hot Jupiter
            } else {
                let weights = [0.20, 0.40, 0.30, 0.10];
                let bins = [
                    (365.0, 1000.0),
                    (1000.0, 3000.0),
                    (3000.0, 5000.0),
                    (5000.0, 10000.0),
                ];
                sample_from_bins(rng, &weights, &bins)
            }
        }
    }
}

fn sample_from_bins(rng: &mut ChaChaRng, weights: &[f64], bins: &[(f64, f64)]) -> f64 {
    let total: f64 = weights.iter().sum();
    let roll: f64 = rng.random::<f64>() * total;
    let mut cumulative = 0.0;

    for (i, &weight) in weights.iter().enumerate() {
        cumulative += weight;
        if roll < cumulative {
            let (min, max) = bins[i];
            let log_min = min.ln();
            let log_max = max.ln();
            return (log_min + rng.random::<f64>() * (log_max - log_min)).exp();
        }
    }

    bins[0].0
}

/// Convert orbital period to semi-major axis (Kepler's third law)
pub fn period_to_semi_major_axis(period_days: f64, stellar_mass_solar: f64) -> f64 {
    let period_years = period_days / 365.25;
    (period_years.powi(2) * stellar_mass_solar).powf(1.0 / 3.0)
}

/// Sample eccentricity
pub fn sample_eccentricity(
    rng: &mut ChaChaRng,
    planet_class: &PlanetClass,
    period_days: f64,
) -> f64 {
    let tidal_damping = match period_days {
        p if p < 10.0 => 0.1,
        p if p < 30.0 => 0.5,
        _ => 1.0,
    };

    let (mean, sigma) = match planet_class {
        PlanetClass::Compact => (0.05, 0.03),
        PlanetClass::Transitional => (0.08, 0.05),
        PlanetClass::Volatile => (0.10, 0.08),
        PlanetClass::Giant => (0.15, 0.12),
    };

    sample_rayleigh(rng, mean * tidal_damping, sigma * tidal_damping).clamp(0.0, 0.8)
}

/// Sample inclination (radians)
pub fn sample_inclination(rng: &mut ChaChaRng, is_multiplanet: bool) -> f64 {
    let sigma_deg: f64 = if is_multiplanet { 1.5 } else { 5.0 };
    sample_rayleigh(rng, 0.0, sigma_deg.to_radians()).abs()
}

fn sample_rayleigh(rng: &mut ChaChaRng, mode: f64, sigma: f64) -> f64 {
    let u: f64 = rng.random();
    mode + sigma * (-2.0 * (1.0 - u).ln()).sqrt()
}

/// Sample number of planets for a system
pub fn sample_planet_count(rng: &mut ChaChaRng, spectral_type: &str, metallicity: f64) -> usize {
    let (mean, std) = match spectral_type {
        "M" => (2.5, 1.5),
        "K" => (2.0, 1.2),
        "G" => (1.5, 1.0),
        "F" => (1.2, 0.8),
        _ => (1.0, 1.0),
    };

    let adjusted_mean = mean * (1.0 + metallicity);
    let n = (adjusted_mean + sample_gaussian(rng, 0.0, std)).round() as i32;
    n.max(0) as usize
}

fn sample_gaussian(rng: &mut ChaChaRng, mean: f64, std_dev: f64) -> f64 {
    let u1: f64 = rng.random();
    let u2: f64 = rng.random();
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
    mean + std_dev * z
}
