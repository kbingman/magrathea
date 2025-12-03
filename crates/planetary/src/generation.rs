//! Planetary system generation pipeline
//!
//! Generates statistically realistic planetary systems based on occurrence rates from:
//! - Kepler/TESS surveys (inner system, P < 400 days)
//! - Radial velocity surveys (cold giants, 1-10 AU)
//! - Microlensing (ice giants near snow line)
//! - Direct imaging (wide companions, 50-1000 AU)
//!
//! See `docs/OCCURRENCE_RATES.md` for detailed calibration.

use rand::Rng;
use rand_chacha::ChaChaRng;
use stellar::MainSequenceStar;
use units::{EARTH_MASS_G, Length, Mass, SOLAR_MASS_G};

/// Earth masses per solar mass (M☉/M⊕)
const EARTH_MASSES_PER_SOLAR: f64 = SOLAR_MASS_G / EARTH_MASS_G;

use crate::composition::Composition;
use crate::planet::Planet;
use crate::planet_class::PlanetClass;
use crate::sampling::{
    period_to_semi_major_axis, sample_eccentricity, sample_inclination, sample_orbital_period,
    sample_planet_mass,
};
use crate::system::{HabitableZone, PlanetarySystem, SystemArchitecture, snow_line};

// =============================================================================
// Outer System Occurrence Rates (from RV, microlensing, direct imaging)
// =============================================================================

/// Base probability for cold giants (1.5-6× snow line) at solar metallicity
/// From Mayor+ 2011, Cumming+ 2008: ~10-15% of FGK stars
const COLD_GIANT_BASE_RATE: f64 = 0.15;

/// Base probability for ice giants (5-15× snow line)
/// From Cassan+ 2012, Suzuki+ 2016: ~30-40% of systems
const ICE_GIANT_BASE_RATE: f64 = 0.35;

/// Probability of wide companions (50-300 AU)
/// From Nielsen+ 2019, Vigan+ 2021: ~1-3% of systems
const WIDE_COMPANION_RATE: f64 = 0.02;

/// Maximum planet-to-star mass ratio
/// Prevents unrealistic massive planets around low-mass stars
/// ~1% is roughly the brown dwarf boundary for M-dwarfs
const MAX_PLANET_STAR_MASS_RATIO: f64 = 0.01;

// =============================================================================
// Stellar Context
// =============================================================================

/// Stellar properties needed for planet generation
///
/// Groups stellar attributes to reduce function parameter counts.
#[derive(Debug, Clone, Copy)]
struct StellarContext {
    /// Stellar mass in solar masses (M☉)
    mass: f64,
    /// Stellar luminosity in solar luminosities (L☉)
    luminosity: f64,
    /// Stellar metallicity [Fe/H] in dex
    metallicity: f64,
}

impl StellarContext {
    fn new(mass: f64, luminosity: f64, metallicity: f64) -> Self {
        Self {
            mass,
            luminosity,
            metallicity,
        }
    }

    /// Snow line distance in AU
    fn snow_line(&self) -> f64 {
        snow_line(self.luminosity)
    }

    /// Maximum planet mass in Earth masses based on stellar mass
    ///
    /// Caps planet mass at ~1% of stellar mass to prevent unrealistic
    /// massive planets around low-mass stars.
    fn max_planet_mass(&self) -> f64 {
        self.mass * EARTH_MASSES_PER_SOLAR * MAX_PLANET_STAR_MASS_RATIO
    }

    /// Giant planet occurrence scaling factor based on stellar mass
    ///
    /// Johnson et al. (2010) found giant occurrence scales roughly as M_star^1.0-1.5
    /// This means M-dwarfs should have ~10× fewer giants than Sun-like stars.
    fn giant_occurrence_scaling(&self) -> f64 {
        // Scale relative to solar mass with exponent ~1.5
        // For 0.1 M☉ star: 0.1^1.5 = 0.03 (3% of solar rate)
        // For 0.5 M☉ star: 0.5^1.5 = 0.35 (35% of solar rate)
        self.mass.powf(1.5)
    }
}

/// Generate a complete planetary system
pub fn generate_planetary_system(
    rng: &mut ChaChaRng,
    stellar_mass: f64,
    stellar_luminosity: f64,
    stellar_temperature: f64,
    stellar_metallicity: f64,
    spectral_type: &str,
) -> PlanetarySystem {
    let star = StellarContext::new(stellar_mass, stellar_luminosity, stellar_metallicity);
    let architecture = SystemArchitecture::sample(rng, spectral_type, stellar_metallicity);

    let planets = match architecture {
        SystemArchitecture::CompactMulti => generate_compact_system(rng, &star),
        SystemArchitecture::Mixed => generate_mixed_system(rng, &star),
        SystemArchitecture::GiantDominated => generate_giant_system(rng, &star),
        SystemArchitecture::Sparse => generate_sparse_system(rng, &star),
    };

    PlanetarySystem::new(
        stellar_mass,
        stellar_luminosity,
        stellar_temperature,
        stellar_metallicity,
        spectral_type.to_string(),
        planets,
        architecture,
    )
}

/// Generate a planetary system from a `MainSequenceStar`
///
/// This is a convenience function that extracts stellar properties from a
/// `MainSequenceStar` and delegates to `generate_planetary_system`.
///
/// # Example
/// ```ignore
/// use rand::SeedableRng;
/// use rand_chacha::ChaChaRng;
/// use stellar::sample_main_sequence_star;
/// use planetary::from_star;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let star = sample_main_sequence_star(&mut rng);
/// let system = from_star(&mut rng, &star);
/// ```
pub fn from_star(rng: &mut ChaChaRng, star: &MainSequenceStar) -> PlanetarySystem {
    let stellar_mass = star.mass.to_solar_masses();
    let stellar_luminosity = star.luminosity;
    let stellar_temperature = star.temperature.to_kelvin();
    let stellar_metallicity = star.metallicity;
    let spectral_type = format!("{}", star.spectral_type);

    generate_planetary_system(
        rng,
        stellar_mass,
        stellar_luminosity,
        stellar_temperature,
        stellar_metallicity,
        &spectral_type,
    )
}

/// Generate a compact multi-planet system (Kepler-like)
///
/// Compact systems have 4-7 tightly-packed inner planets with low
/// mutual inclinations. They can still have outer system companions,
/// but at reduced rates compared to mixed systems.
fn generate_compact_system(rng: &mut ChaChaRng, star: &StellarContext) -> Vec<Planet> {
    let n_planets: usize = rng.random_range(4..=7);
    let hz = HabitableZone::from_luminosity(star.luminosity);

    let inner_au = 0.01 * star.luminosity.sqrt();
    let outer_au = hz.outer_edge * 1.5;

    let mut planets =
        generate_n_spaced_planets(rng, star, n_planets, inner_au, outer_au, 0.01..5.0);

    // Compact systems have reduced outer planet rates
    // Cold giants might destabilize inner compact systems
    // Modifiers: 0.5× for cold giants, 0.7× for ice giants
    let outer = generate_outer_system(rng, star, 0.5, 0.7);
    planets.extend(outer);

    planets
}

/// Generate a mixed architecture system (inner terrestrials + outer giants)
///
/// These are "Solar System-like" architectures with small rocky/icy planets
/// in the inner system and gas/ice giants in the outer system. This is
/// likely the most common architecture for FGK stars with giant planets.
fn generate_mixed_system(rng: &mut ChaChaRng, star: &StellarContext) -> Vec<Planet> {
    let sl = star.snow_line();

    // Inner terrestrial zone: inside snow line
    let n_inner: usize = rng.random_range(1..=4);
    let inner_au = 0.3 * star.luminosity.sqrt();
    let mut planets = generate_n_spaced_planets(rng, star, n_inner, inner_au, sl * 0.8, 0.05..5.0);

    // Mixed systems have enhanced outer system rates
    // These are the "full" planetary systems
    // Modifiers: 1.5× for cold giants (they define this architecture), 1.2× for ice giants
    let outer = generate_outer_system(rng, star, 1.5, 1.2);
    planets.extend(outer);

    planets
}

/// Generate a giant-dominated system
///
/// These systems are defined by having one or more massive gas giants.
/// Two sub-types:
/// - Hot Jupiter systems (~20%): Migrated giant close-in, usually "lonely"
/// - Cold Jupiter systems (~80%): Giant(s) near snow line, may have ice giants
fn generate_giant_system(rng: &mut ChaChaRng, star: &StellarContext) -> Vec<Planet> {
    let sl = star.snow_line();
    let max_mass = star.max_planet_mass();

    // For low-mass stars, cap giant mass appropriately
    // Minimum mass for a "giant" is ~50 M⊕ (Saturn-class)
    let min_giant_mass = 50.0_f64;
    if max_mass < min_giant_mass {
        // Star too small for true giants - generate sub-giants instead
        return generate_n_spaced_planets(rng, star, 1, sl * 1.0, sl * 4.0, 10.0..max_mass);
    }

    let hot_jupiter = rng.random::<f64>() < 0.2;

    // Scale giant mass range by stellar mass
    let giant_max = max_mass.min(1000.0);
    let giant_min = min_giant_mass.min(giant_max * 0.5);

    let mut planets = if hot_jupiter {
        // Hot Jupiter: P < 10 days, migrated inward
        // Hot Jupiters are "lonely" - migration clears inner system
        let mass = sample_giant_mass_for_star(rng, star);
        let period = rng.random_range(2.0..10.0);
        let sma = period_to_semi_major_axis(period, star.mass);
        vec![create_planet(rng, star, mass, sma)]
    } else {
        // Cold giant(s) in Jupiter zone
        let n_giants: usize = rng.random_range(1..=2);
        generate_n_spaced_planets(
            rng,
            star,
            n_giants,
            sl * 1.5,
            sl * 6.0,
            giant_min..giant_max,
        )
    };

    // Cold giant systems often have ice giants further out
    // Enhanced rate (~50%) since giant formation was successful
    if !hot_jupiter {
        let ice_giants = generate_ice_giants(rng, star, ICE_GIANT_BASE_RATE * 1.4);
        planets.extend(ice_giants);

        // Rare wide companions (slightly enhanced for giant systems)
        let wide = generate_wide_companion(rng, star);
        planets.extend(wide);
    }

    // Surviving inner planets (rare with hot jupiter due to migration clearing)
    let inner_prob = if hot_jupiter { 0.1 } else { 0.3 };
    if rng.random::<f64>() < inner_prob {
        let n_inner = rng.random_range(1..=2);
        let inner = generate_n_spaced_planets(rng, star, n_inner, 0.3, sl * 0.5, 0.1..3.0);
        planets = [inner, planets].concat();
    }

    planets
}

/// Generate a sparse system (0-1 detected planets)
///
/// These systems either have no detectable planets or just one lonely planet.
/// ~50% are truly empty (or have only very small/distant planets).
/// May still have outer system planets at reduced rates.
fn generate_sparse_system(rng: &mut ChaChaRng, star: &StellarContext) -> Vec<Planet> {
    let mut planets = Vec::new();

    // 50% chance of having an inner planet
    if rng.random::<f64>() < 0.5 {
        let mass = sample_planet_mass(rng, star.metallicity);
        let class = PlanetClass::from_earth_masses(mass);
        let period = sample_orbital_period(rng, &class);
        let sma = period_to_semi_major_axis(period, star.mass);
        planets.push(create_planet(rng, star, mass, sma));
    }

    // Sparse systems can still have outer planets, but at reduced rates
    // Modifiers: 0.3× for cold giants, 0.4× for ice giants
    let outer = generate_outer_system(rng, star, 0.3, 0.4);
    planets.extend(outer);

    planets
}

/// Generate the outer system: cold giants + ice giants + rare wide companions
///
/// Outer system zones (scaled by snow line, ~2.7 AU for Sun):
/// - Jupiter zone: 1.5-4× SL (~4-11 AU) - gas giants peak here
/// - Saturn zone: 3-7× SL (~8-19 AU) - secondary giants, large ice giants
/// - Ice giant zone: 5-15× SL (~14-40 AU) - Uranus/Neptune analogs
/// - Wide companion zone: 20-110× SL (~50-300 AU) - rare massive planets
///
/// # Arguments
/// * `cold_giant_modifier` - Multiplier on base cold giant rate (architecture-dependent)
/// * `ice_giant_modifier` - Multiplier on base ice giant rate
fn generate_outer_system(
    rng: &mut ChaChaRng,
    star: &StellarContext,
    cold_giant_modifier: f64,
    ice_giant_modifier: f64,
) -> Vec<Planet> {
    let mut planets = Vec::new();
    let sl = star.snow_line();

    // Cold giant(s) in the Jupiter zone (1.5-6× snow line)
    // Strong metallicity dependence: P(giant) ∝ 10^(2×[Fe/H])
    // Also scales with stellar mass: Johnson+ 2010 found P(giant) ∝ M_star^1-1.5
    let metallicity_boost = 10.0_f64.powf(2.0 * star.metallicity);
    let stellar_mass_scaling = star.giant_occurrence_scaling();
    let cold_giant_prob =
        (COLD_GIANT_BASE_RATE * cold_giant_modifier * metallicity_boost * stellar_mass_scaling)
            .min(0.9);

    let giant_roll = rng.random::<f64>();
    if giant_roll < cold_giant_prob {
        // Of systems with cold giants: ~70% have 1, ~30% have 2
        let n_giants = if giant_roll < cold_giant_prob * 0.3 {
            2
        } else {
            1
        };

        let giants =
            generate_n_spaced_planets(rng, star, n_giants, sl * 1.5, sl * 6.0, 80.0..800.0);
        planets.extend(giants);
    }

    // Ice giants in the outer zone (5-15× snow line)
    // Weak metallicity dependence
    let ice_giant_prob = ICE_GIANT_BASE_RATE * ice_giant_modifier;
    let ice = generate_ice_giants(rng, star, ice_giant_prob);
    planets.extend(ice);

    // Wide companions (50-300 AU) - rare, massive planets
    // May form via disk instability or capture; no metallicity dependence
    let wide = generate_wide_companion(rng, star);
    planets.extend(wide);

    planets
}

/// Generate ice giants in the outer system
///
/// Ice giants form beyond the snow line where slower accretion prevented
/// runaway gas accretion. Microlensing surveys (Cassan+ 2012, Suzuki+ 2016)
/// suggest cold Neptunes are very common (~50% of stars).
///
/// Zone: 5-15× snow line (~14-40 AU for Sun-like star)
/// Mass range: 8-30 M⊕ (sub-Neptune to super-Neptune)
///
/// Multiplicity distribution (of systems that have ice giants):
/// - 50% have 1 ice giant
/// - 40% have 2 ice giants (like Solar System)
/// - 10% have 3 ice giants
fn generate_ice_giants(
    rng: &mut ChaChaRng,
    star: &StellarContext,
    probability: f64,
) -> Vec<Planet> {
    // Ice giants have weaker stellar mass dependence than gas giants
    // Use sqrt scaling (exponent 0.5) instead of 1.5
    let stellar_scaling = star.mass.powf(0.5);
    let scaled_prob = probability * stellar_scaling;

    if rng.random::<f64>() > scaled_prob {
        return vec![];
    }

    let sl = star.snow_line();

    // Ice giant zone: 5-15× snow line
    // Uranus ~19 AU (~7× SL), Neptune ~30 AU (~11× SL) for our solar system
    let inner_ice = sl * 5.0;
    let outer_ice = sl * 15.0;

    // Multiplicity: 50% have 1, 40% have 2, 10% have 3
    let n_ice_giants = match rng.random::<f64>() {
        x if x < 0.50 => 1,
        x if x < 0.90 => 2,
        _ => 3,
    };

    // Use a zero-metallicity context for ice giants (weak metallicity dependence)
    let ice_star = StellarContext::new(star.mass, star.luminosity, 0.0);
    generate_n_spaced_planets(
        rng,
        &ice_star,
        n_ice_giants,
        inner_ice,
        outer_ice,
        8.0..30.0,
    )
}

/// Generate rare wide-separation companions (50-300 AU)
///
/// Direct imaging surveys (GPIES, SPHERE/SHINE) find that ~1-3% of FGK stars
/// have massive planetary companions at wide separations. These may form via:
/// - Gravitational instability in massive disks
/// - Scattering from inner system
/// - Capture during cluster dissolution
///
/// Zone: 50-300 AU (~20-110× snow line for Sun)
/// Mass: 1-13 M_J (above detection threshold, below brown dwarf)
///
/// # References
/// - Nielsen+ 2019 (GPIES): 0.8 ± 0.5% for 5-13 M_J at 10-100 AU
/// - Vigan+ 2021 (SHINE): 1-3% for 1-20 M_J at 10-100 AU
fn generate_wide_companion(rng: &mut ChaChaRng, star: &StellarContext) -> Vec<Planet> {
    // Wide companions also scale with stellar mass (disk instability needs massive disk)
    let scaled_rate = WIDE_COMPANION_RATE * star.giant_occurrence_scaling();
    if rng.random::<f64>() > scaled_rate {
        return vec![];
    }

    // Wide companion zone: 50-300 AU
    // Semi-major axis follows roughly log-uniform distribution
    let inner_wide: f64 = 50.0;
    let outer_wide: f64 = 300.0;
    let log_inner = inner_wide.ln();
    let log_outer = outer_wide.ln();
    let sma = (log_inner + rng.random::<f64>() * (log_outer - log_inner)).exp();

    // Mass distribution: 1-13 M_J, weighted toward lower masses
    // Direct imaging is biased toward detecting more massive planets
    // True distribution likely falls off steeply with mass
    // Cap at stellar mass limit
    let min_mass: f64 = 318.0; // 1 M_J in Earth masses
    let max_mass: f64 = (318.0_f64 * 13.0).min(star.max_planet_mass()); // 13 M_J or stellar limit

    // If stellar mass limit is below 1 M_J, skip wide companion
    if max_mass < min_mass {
        return vec![];
    }

    let log_min = min_mass.ln();
    let log_max = max_mass.ln();
    // Power law favoring lower masses: dN/dM ∝ M^(-1.3)
    let u: f64 = rng.random();
    let mass = (log_min + u.powf(1.3) * (log_max - log_min)).exp();

    vec![create_planet(rng, star, mass, sma)]
}

fn generate_n_spaced_planets(
    rng: &mut ChaChaRng,
    star: &StellarContext,
    n_planets: usize,
    inner_au: f64,
    outer_au: f64,
    mass_range: std::ops::Range<f64>,
) -> Vec<Planet> {
    if n_planets == 0 {
        return vec![];
    }

    for _ in 0..50 {
        let mut planets = Vec::with_capacity(n_planets);

        let log_inner = inner_au.ln();
        let log_outer = outer_au.ln();
        let log_spacing = (log_outer - log_inner) / n_planets as f64;

        for i in 0..n_planets {
            let log_base = log_inner + log_spacing * (i as f64 + 0.5);
            let scatter = rng.random_range(-0.3..0.3) * log_spacing;
            let sma = (log_base + scatter).exp();

            let mass =
                sample_mass_in_range(rng, mass_range.start, mass_range.end, star.metallicity);
            planets.push(create_planet(rng, star, mass, sma));
        }

        planets.sort_by(|a, b| {
            a.semi_major_axis
                .to_au()
                .partial_cmp(&b.semi_major_axis.to_au())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        if is_stable(&planets, star.mass) {
            return planets;
        }
    }

    vec![] // Fallback
}

fn sample_mass_in_range(rng: &mut ChaChaRng, min: f64, max: f64, _metallicity: f64) -> f64 {
    let log_min = min.max(0.01).ln();
    let log_max = max.ln();
    let u: f64 = rng.random();
    (log_min + u.powf(0.7) * (log_max - log_min)).exp()
}

fn sample_giant_mass(rng: &mut ChaChaRng) -> f64 {
    let log_jupiter = 318.0_f64.ln();
    let z: f64 = rng.random_range(-0.5..0.5);
    (log_jupiter + z).exp().clamp(50.0, 2000.0)
}

/// Sample giant mass scaled to stellar mass limits
fn sample_giant_mass_for_star(rng: &mut ChaChaRng, star: &StellarContext) -> f64 {
    let max_mass = star.max_planet_mass();
    let base_mass = sample_giant_mass(rng);
    base_mass.min(max_mass)
}

fn create_planet(
    rng: &mut ChaChaRng,
    star: &StellarContext,
    mass_earth: f64,
    sma_au: f64,
) -> Planet {
    // Cap planet mass at ~1% of stellar mass to prevent unrealistic planets
    let mass_earth = mass_earth.min(star.max_planet_mass());
    let class = PlanetClass::from_earth_masses(mass_earth);
    let period = (sma_au.powi(3) / star.mass).sqrt() * 365.25;
    let eccentricity = sample_eccentricity(rng, &class, period);
    let inclination = sample_inclination(rng, true);
    let composition = sample_composition(rng, &class, sma_au, star.luminosity);

    Planet::from_mass(
        Mass::from_earth_masses(mass_earth),
        Length::from_au(sma_au),
        eccentricity,
        inclination,
        composition,
        star.luminosity,
        rng,
    )
}

fn sample_composition(
    rng: &mut ChaChaRng,
    class: &PlanetClass,
    sma_au: f64,
    stellar_luminosity: f64,
) -> Composition {
    let sl = snow_line(stellar_luminosity);
    let beyond_snow_line = sma_au > sl;

    match class {
        PlanetClass::Rocky => {
            if beyond_snow_line {
                Composition::new(
                    0.15 + rng.random::<f64>() * 0.10,
                    0.35 + rng.random::<f64>() * 0.10,
                    0.40 + rng.random::<f64>() * 0.20,
                    0.0,
                )
            } else {
                Composition::new(
                    0.25 + rng.random::<f64>() * 0.15,
                    0.60 + rng.random::<f64>() * 0.15,
                    rng.random::<f64>() * 0.05,
                    0.0,
                )
            }
        }
        PlanetClass::Transitional => {
            let water = if beyond_snow_line {
                0.20 + rng.random::<f64>() * 0.30
            } else {
                rng.random::<f64>() * 0.10
            };
            let envelope = rng.random::<f64>() * 0.15;
            Composition::new(
                0.15 + rng.random::<f64>() * 0.10,
                0.40 + rng.random::<f64>() * 0.10,
                water,
                envelope,
            )
        }
        PlanetClass::Volatile => Composition::ice_giant(),
        PlanetClass::Giant => Composition::gas_giant(),
    }
}

fn is_stable(planets: &[Planet], stellar_mass: f64) -> bool {
    if planets.len() < 2 {
        return true;
    }

    for window in planets.windows(2) {
        let m1 = window[0].mass.to_earth_masses() / EARTH_MASSES_PER_SOLAR;
        let m2 = window[1].mass.to_earth_masses() / EARTH_MASSES_PER_SOLAR;
        let a1 = window[0].semi_major_axis.to_au();
        let a2 = window[1].semi_major_axis.to_au();

        let mutual_hill = ((m1 + m2) / (3.0 * stellar_mass)).powf(1.0 / 3.0) * ((a1 + a2) / 2.0);
        if (a2 - a1) / mutual_hill < 8.0 {
            return false;
        }
    }

    true
}
