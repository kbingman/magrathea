//! Factory functions for generating stellar objects.
//!
//! This module contains functions for creating stars and stellar remnants
//! with physically realistic properties based on mass, age, and metallicity.

use rand::Rng;
use rand_chacha::ChaChaRng;
use units::time::Time;
use units::{Mass, Temperature};

use super::sampling::{sample_gaussian, sample_mass_kroupa, sample_metallicity};
use super::spectral::{LuminosityClass, SpectralType, VariabilityType};
use super::stellar_color::StellarColor;
use super::stellar_objects::{
    BlackHole, GiantStar, MainSequenceStar, NeutronStar, StellarObject, WhiteDwarf, WhiteDwarfType,
};
use super::stellar_radius::StellarRadius;

/// Solar temperature in Kelvin
const SOLAR_TEMP: f64 = 5778.0;

/// Temperature boundaries for spectral types (in Kelvin)
const TEMP_BOUNDS: [(SpectralType, f64); 10] = [
    (SpectralType::O, 30000.0),
    (SpectralType::B, 10000.0),
    (SpectralType::A, 7500.0),
    (SpectralType::F, 6000.0),
    (SpectralType::G, 5200.0),
    (SpectralType::K, 3700.0),
    (SpectralType::M, 2400.0),
    (SpectralType::L, 1300.0),
    (SpectralType::T, 550.0),
    (SpectralType::Y, 0.0),
];

// ============================================================================
// Main Sequence Stars
// ============================================================================

/// Creates a new main sequence star with physically realistic properties based on mass.
///
/// Uses mass-luminosity and mass-temperature relations to derive all properties.
///
/// # Mass Categories and Properties
/// * Very Massive (>30 M☉): O-type stars, L ∝ M^3.5, T ~ 38,000K
/// * Massive (8-30 M☉): B-type stars, L ∝ M^3.5, T ~ 22,000K
/// * Intermediate (2-8 M☉): A/F-type stars, L ∝ M^3.5, T ~ 9,000K
/// * Solar-type (0.8-2 M☉): G-type stars, L ∝ M^4, T ~ 5,800K
/// * Low Mass (<0.8 M☉): K/M-type stars, L ∝ M^2.3, T ~ 2,500-4,800K
///
/// # Arguments
/// * `mass_solar` - Stellar mass in solar masses
/// * `metallicity` - Metallicity [Fe/H] in dex (0.0 = solar)
/// * `age_myr` - Stellar age in millions of years
///
/// # Example
/// ```
/// use stellar_forge::stellar::generation::main_sequence_star;
///
/// let sun_like = main_sequence_star(1.0, 0.0, 4600.0);
/// ```
pub fn main_sequence_star(mass_solar: f64, metallicity: f64, age_myr: f64) -> MainSequenceStar {
    let (luminosity, temperature) = calculate_ms_properties(mass_solar);

    let radius_solar = calculate_radius(luminosity, temperature);
    let spectral_type = spectral_type_from_temp(temperature);
    let subtype = calculate_subtype(temperature);

    MainSequenceStar {
        mass: Mass::from_solar_masses(mass_solar),
        spectral_type,
        luminosity,
        temperature: Temperature::from_kelvin(temperature),
        luminosity_class: LuminosityClass::V,
        radius: StellarRadius::from_solar_radii(radius_solar),
        subtype,
        variability: VariabilityType::determine_variability(
            mass_solar,
            temperature,
            LuminosityClass::V,
        ),
        metallicity,
        age: Time::from_myr(age_myr),
        color: StellarColor::from_temperature(temperature),
    }
}

/// Create a solar analog (1 M☉, solar metallicity, 4.6 Gyr age)
pub fn solar_analog() -> MainSequenceStar {
    main_sequence_star(1.0, 0.0, 1.0)
}

/// Sample a random main sequence star from a realistic stellar population.
///
/// Uses the Kroupa (2001) IMF for mass distribution, limited to main sequence
/// stars capable of hosting planetary systems (0.08-3.0 M☉).
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar_forge::stellar::generation::sample_main_sequence_star;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let star = sample_main_sequence_star(&mut rng);
/// ```
pub fn sample_main_sequence_star(rng: &mut ChaChaRng) -> MainSequenceStar {
    // Sample mass from Kroupa IMF, limited to planet-hosting range (≤3 M☉)
    let mass = sample_mass_kroupa(rng, 3.0);
    let metallicity = sample_metallicity(rng);
    let age_myr = sample_young_star_age(rng);

    main_sequence_star(mass, metallicity, age_myr)
}

/// Calculate luminosity and temperature for a main sequence star of given mass.
fn calculate_ms_properties(mass_solar: f64) -> (f64, f64) {
    match mass_solar {
        // Very massive main sequence (> 30 M☉)
        m if m > 30.0 => {
            let luminosity = 3.0e4 * (m / 30.0).powf(3.5);
            let temperature = 38000.0 * (m / 30.0).powf(0.2);
            (luminosity, temperature)
        }
        // Massive main sequence (8-30 M☉)
        m if m > 8.0 => {
            let luminosity = 1.0e3 * (m / 8.0).powf(3.5);
            let temperature = 22000.0 * (m / 8.0).powf(0.2);
            (luminosity, temperature)
        }
        // Intermediate mass main sequence (2-8 M☉)
        m if m > 2.0 => {
            let luminosity = 25.0 * (m / 2.0).powf(3.5);
            let temperature = 9000.0 * (m / 2.0).powf(0.2);
            (luminosity, temperature)
        }
        // Solar-type main sequence (0.8-2 M☉)
        m if m > 0.8 => {
            let luminosity = m.powf(4.0);
            let temperature = 5800.0 * m.powf(0.1);
            (luminosity, temperature)
        }
        // Low mass main sequence (< 0.8 M☉)
        m => {
            let luminosity = m.powf(2.3);
            let temperature = if m < 0.45 {
                2500.0 * (m / 0.08).powf(0.23)
            } else {
                3700.0 * (m / 0.45).powf(0.45)
            };
            (luminosity, temperature)
        }
    }
}

/// Sample stellar age appropriate for planet formation studies (1-10 Myr).
fn sample_young_star_age(rng: &mut ChaChaRng) -> f64 {
    1.0 + rng.random::<f64>() * 9.0
}

// ============================================================================
// Giant Stars
// ============================================================================

/// Creates a giant star based on initial mass with physically realistic properties.
///
/// # Mass Categories
/// * Hypergiants (>30 M☉): Class Ia+/Ia, L ~ 5×10^5 L☉
/// * Supergiants (8-30 M☉): Class Ia/Ib, L ~ 5×10^4 L☉
/// * Regular Giants (2-8 M☉): Class II/III, L ~ 100 L☉
/// * Small Giants (0.8-2 M☉): Class III/IV, L ~ 10 L☉
///
/// # Arguments
/// * `rng` - Random number generator for property variation
/// * `initial_mass` - Initial mass in solar masses
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar_forge::stellar::generation::giant_star;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let supergiant = giant_star(&mut rng, 15.0);
/// ```
pub fn giant_star(rng: &mut ChaChaRng, initial_mass: f64) -> GiantStar {
    let (mass_solar, luminosity, temp_kelvin, luminosity_class) = match initial_mass {
        // Hypergiant (> 30 M☉)
        m if m > 30.0 => {
            let mass = m * rng.random_range(0.78..0.82);
            let luminosity = 5.0e5 * (m / 30.0).powf(2.0);
            let temperature = rng.random_range(4000.0..30000.0);
            let luminosity_class = if luminosity > 500000.0 {
                LuminosityClass::IAPLUS
            } else {
                LuminosityClass::IA
            };
            (mass, luminosity, temperature, luminosity_class)
        }
        // Supergiant (8-30 M☉)
        m if m > 8.0 => {
            let mass = m * 0.85;
            let luminosity = 5.0e4 * (m / 8.0).powf(2.0);
            let temperature = rng.random_range(3500.0..20000.0);
            let luminosity_class = if luminosity > 50000.0 {
                LuminosityClass::IA
            } else {
                LuminosityClass::IB
            };
            (mass, luminosity, temperature, luminosity_class)
        }
        // Regular giant (2-8 M☉)
        m if m > 2.0 => {
            let mass = m * 0.95;
            let luminosity = 100.0 * (m / 2.0).powf(2.0);
            let temperature = 4500.0;
            let luminosity_class = if luminosity > 1000.0 {
                LuminosityClass::II
            } else {
                LuminosityClass::III
            };
            (mass, luminosity, temperature, luminosity_class)
        }
        // Smaller giants (0.8-2 M☉)
        _ => {
            let mass = initial_mass * 0.98;
            let luminosity = 10.0 * initial_mass.powf(2.0);
            let temperature = 4800.0;
            let luminosity_class = if luminosity > 10.0 {
                LuminosityClass::III
            } else {
                LuminosityClass::IV
            };
            (mass, luminosity, temperature, luminosity_class)
        }
    };

    let radius_solar = calculate_radius(luminosity, temp_kelvin);
    let spectral_type = spectral_type_from_temp(temp_kelvin);
    let subtype = calculate_subtype(temp_kelvin);

    GiantStar {
        mass: Mass::from_solar_masses(mass_solar),
        radius: StellarRadius::from_solar_radii(radius_solar),
        luminosity,
        temperature: Temperature::from_kelvin(temp_kelvin),
        spectral_type,
        subtype,
        luminosity_class,
        variability: VariabilityType::determine_variability(
            mass_solar,
            temp_kelvin,
            luminosity_class,
        ),
        color: StellarColor::from_temperature(temp_kelvin),
    }
}

// ============================================================================
// White Dwarfs
// ============================================================================

/// Creates a white dwarf with physically accurate properties.
///
/// Mass is randomly sampled from 0.17-1.44 M☉ (Chandrasekhar limit).
/// Radius follows the mass-radius relation for degenerate matter.
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar_forge::stellar::generation::white_dwarf;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let wd = white_dwarf(&mut rng);
/// ```
pub fn white_dwarf(rng: &mut ChaChaRng) -> WhiteDwarf {
    let mass_solar: f64 = rng.random_range(0.17..=1.44);
    let radius_solar = 0.01 * (0.6 / mass_solar).powf(1.0 / 3.0);
    let base_luminosity = 0.001;
    let luminosity = base_luminosity * (mass_solar / 0.6).powf(1.0);
    let spectral_type = determine_wd_spectral_class(luminosity);

    // Calculate temperature from Stefan-Boltzmann: T/T_sun = (L/R²)^(1/4)
    let temp_kelvin = SOLAR_TEMP * (luminosity / (radius_solar * radius_solar)).powf(0.25);

    WhiteDwarf {
        mass: Mass::from_solar_masses(mass_solar),
        radius: StellarRadius::from_solar_radii(radius_solar),
        luminosity,
        temperature: Temperature::from_kelvin(temp_kelvin),
        spectral_type,
        color: StellarColor::from_temperature(temp_kelvin),
    }
}

fn determine_wd_spectral_class(luminosity: f64) -> WhiteDwarfType {
    match luminosity {
        l if l > 0.01 => WhiteDwarfType::DA,
        l if l > 0.001 => WhiteDwarfType::DB,
        l if l > 0.0001 => WhiteDwarfType::DC,
        l if l > 0.00001 => WhiteDwarfType::DQ,
        _ => WhiteDwarfType::DZ,
    }
}

// ============================================================================
// Neutron Stars
// ============================================================================

/// Creates a neutron star with physically realistic properties.
///
/// Mass follows the Özel & Freire (2016) distribution peaked at 1.4 M☉.
/// ~10% are magnetars with fields > 10^14 G.
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar_forge::stellar::generation::neutron_star;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let ns = neutron_star(&mut rng);
/// ```
pub fn neutron_star(rng: &mut ChaChaRng) -> NeutronStar {
    let mass_solar = sample_ns_mass(rng);
    let radius_km = calculate_ns_radius_km(mass_solar);

    // Magnetic field distribution is bimodal
    let magnetic_field = if rng.random_bool(0.1) {
        rng.random_range(14.0..=15.0) // Magnetar
    } else {
        rng.random_range(11.0..=13.0) // Regular NS
    };

    let magnetar = magnetic_field >= 14.0;
    let pulsar = if magnetar {
        rng.random_bool(0.3)
    } else {
        rng.random_bool(0.8)
    };

    let color = if magnetar {
        StellarColor::new(200, 220, 255) // Bright blue-white
    } else if pulsar {
        StellarColor::new(180, 200, 255) // Blue-white
    } else {
        StellarColor::new(160, 180, 220) // Dimmer blue
    };

    NeutronStar {
        mass: Mass::from_solar_masses(mass_solar),
        radius: StellarRadius::from_km(radius_km),
        magnetic_field,
        pulsar,
        magnetar,
        color,
    }
}

fn sample_ns_mass(rng: &mut ChaChaRng) -> f64 {
    // Skewed normal distribution (Özel & Freire 2016)
    let z = sample_gaussian(rng, 0.0, 1.0);
    let mean = 1.4;
    let std_dev = 0.2;
    let skew_factor = 0.3;
    let skewed = mean + std_dev * z + skew_factor * z.abs();
    skewed.clamp(1.1, 2.5)
}

fn calculate_ns_radius_km(mass_solar: f64) -> f64 {
    12.0 * (mass_solar / 1.4).powf(-0.3)
}

// ============================================================================
// Black Holes
// ============================================================================

/// Creates a black hole based on progenitor star parameters.
///
/// Final mass depends on initial mass, metallicity (affects wind loss),
/// and supernova mass ejection.
///
/// # Arguments
/// * `rng` - Random number generator
/// * `initial_mass` - Mass of the progenitor star in solar masses
/// * `metallicity` - Metallicity of the progenitor star
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar_forge::stellar::generation::black_hole;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let bh = black_hole(&mut rng, 25.0, 0.02);
/// ```
pub fn black_hole(rng: &mut ChaChaRng, initial_mass: f64, metallicity: f64) -> BlackHole {
    let mass_solar = calculate_bh_mass(rng, initial_mass, metallicity);
    let spin = calculate_bh_spin(rng, initial_mass);
    let has_accretion = rng.random_bool(0.1);

    let color = if has_accretion {
        StellarColor::new(255, 200, 150) // Orange-white accretion glow
    } else {
        StellarColor::new(20, 20, 30) // Nearly black with faint blue
    };

    BlackHole {
        mass: Mass::from_solar_masses(mass_solar),
        has_accretion,
        spin,
        color,
    }
}

fn calculate_bh_mass(rng: &mut ChaChaRng, initial_mass: f64, metallicity: f64) -> f64 {
    let wind_loss = initial_mass * metallicity * 0.4;
    let sn_loss = initial_mass * rng.random_range(0.2..0.5);
    let final_mass = initial_mass - wind_loss - sn_loss;
    let variation = rng.random_range(0.8..1.2);

    final_mass * variation
}

fn calculate_bh_spin(rng: &mut ChaChaRng, initial_mass: f64) -> f64 {
    let base_spin = (initial_mass - 20.0) / 100.0;
    let variation = rng.random_range(-0.2..0.2);

    (base_spin + variation).clamp(0.0, 1.0)
}

// ============================================================================
// Stellar Object (any type)
// ============================================================================

/// Creates a stellar object based on initial parameters.
///
/// Determines whether the star is still "alive" (main sequence or giant)
/// or has become a remnant (white dwarf, neutron star, black hole) based
/// on the age relative to estimated lifetime.
///
/// # Arguments
/// * `rng` - Random number generator
/// * `initial_mass` - The star's initial mass in solar masses
/// * `age_years` - The current age of the star in years
/// * `metallicity` - The star's metal content [Fe/H] in dex
///
/// # Example
/// ```
/// use rand_chacha::ChaChaRng;
/// use rand::SeedableRng;
/// use stellar_forge::stellar::generation::stellar_object;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let star = stellar_object(&mut rng, 1.0, 1.0e9, 0.0);
/// ```
pub fn stellar_object(
    rng: &mut ChaChaRng,
    initial_mass: f64,
    age_years: f64,
    metallicity: f64,
) -> StellarObject {
    let lifetime_gyr = estimate_lifetime(initial_mass);
    let age_gyr = age_years / 1.0e9;
    let life_fraction = age_gyr / lifetime_gyr;

    if life_fraction > 1.0 {
        return new_remnant(rng, initial_mass, metallicity);
    }

    let age_myr = Time::from_years(age_years).to_myr();

    new_living_star(rng, life_fraction, initial_mass, metallicity, age_myr)
}

/// Sample a stellar object from realistic galactic distributions.
///
/// # Example
/// ```
/// use rand::SeedableRng;
/// use rand_chacha::ChaChaRng;
/// use stellar_forge::stellar::generation::sample_stellar_object;
///
/// let mut rng = ChaChaRng::seed_from_u64(42);
/// let star = sample_stellar_object(&mut rng);
/// ```
pub fn sample_stellar_object(rng: &mut ChaChaRng) -> StellarObject {
    let mass = sample_mass_kroupa(rng, 150.0);
    let metallicity = sample_metallicity(rng);
    let max_lifetime_gyr = estimate_lifetime(mass);
    let age_gyr = sample_age_capped(rng, max_lifetime_gyr);
    let age_years = age_gyr * 1.0e9;

    stellar_object(rng, mass, age_years, metallicity)
}

fn new_living_star(
    rng: &mut ChaChaRng,
    life_fraction: f64,
    initial_mass: f64,
    metallicity: f64,
    age_myr: f64,
) -> StellarObject {
    match initial_mass {
        m if m > 30.0 => {
            if life_fraction < 0.9 {
                StellarObject::MainSequence(main_sequence_star(m, metallicity, age_myr))
            } else {
                StellarObject::Giant(giant_star(rng, m))
            }
        }
        m if m > 8.0 => {
            if life_fraction < 0.9 {
                StellarObject::MainSequence(main_sequence_star(m, metallicity, age_myr))
            } else {
                StellarObject::Giant(giant_star(rng, m))
            }
        }
        m if m > 2.0 => match life_fraction {
            x if x < 0.85 => {
                StellarObject::MainSequence(main_sequence_star(m, metallicity, age_myr))
            }
            x if x < 0.95 => StellarObject::Giant(giant_star(rng, m)),
            _ => StellarObject::WhiteDwarf(white_dwarf(rng)),
        },
        m if m > 0.8 => match life_fraction {
            x if x < 0.90 => {
                StellarObject::MainSequence(main_sequence_star(m, metallicity, age_myr))
            }
            x if x < 0.95 => StellarObject::Giant(giant_star(rng, m)),
            _ => StellarObject::WhiteDwarf(white_dwarf(rng)),
        },
        m => StellarObject::MainSequence(main_sequence_star(m, metallicity, age_myr)),
    }
}

fn new_remnant(rng: &mut ChaChaRng, initial_mass: f64, metallicity: f64) -> StellarObject {
    match initial_mass {
        m if m < 8.0 => StellarObject::WhiteDwarf(white_dwarf(rng)),
        m if m < 20.0 => StellarObject::NeutronStar(neutron_star(rng)),
        _ => StellarObject::BlackHole(black_hole(rng, initial_mass, metallicity)),
    }
}

/// Estimate stellar lifetime in billions of years.
pub fn estimate_lifetime(mass: f64) -> f64 {
    match mass {
        m if m > 10.0 => 0.0001 * (m / 10.0).powf(-2.0),
        m if m > 2.0 => 0.01 * (m / 2.0).powf(-2.5),
        m if m > 0.5 => 10.0 * m.powf(-3.0),
        _ => 100.0 * mass.powf(-2.5),
    }
}

fn sample_age_capped(rng: &mut ChaChaRng, max_lifetime_gyr: f64) -> f64 {
    const THIN_DISK_MAX_AGE: f64 = 10.0;
    const MIN_AGE: f64 = 0.1;

    let effective_max = max_lifetime_gyr.min(THIN_DISK_MAX_AGE);

    if effective_max <= MIN_AGE {
        return MIN_AGE * rng.random::<f64>();
    }

    MIN_AGE + rng.random::<f64>() * (effective_max - MIN_AGE)
}

// ============================================================================
// Shared Utilities
// ============================================================================

/// Calculate stellar radius from luminosity and temperature.
///
/// Uses Stefan-Boltzmann: L = 4πR²σT⁴
pub fn calculate_radius(luminosity: f64, temperature: f64) -> f64 {
    luminosity.powf(0.5) / (temperature / SOLAR_TEMP).powf(2.0)
}

/// Determine spectral type from temperature.
pub fn spectral_type_from_temp(temperature: f64) -> SpectralType {
    for (spec_type, temp_bound) in TEMP_BOUNDS.iter() {
        if temperature >= *temp_bound {
            return *spec_type;
        }
    }
    SpectralType::Y
}

/// Calculate spectral subtype (0-9) from temperature.
pub fn calculate_subtype(temperature: f64) -> u8 {
    let mut upper_bound = TEMP_BOUNDS[0].1;
    let mut lower_bound = TEMP_BOUNDS[1].1;

    for window in TEMP_BOUNDS.windows(2) {
        if temperature >= window[1].1 {
            upper_bound = window[0].1;
            lower_bound = window[1].1;
            break;
        }
    }

    let temp_range = upper_bound - lower_bound;
    let temp_position = upper_bound - temperature;
    let subtype = (9.0 * temp_position / temp_range).round() as u8;
    subtype.clamp(0, 9)
}
