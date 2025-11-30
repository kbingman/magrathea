use std::f64::consts::PI;

use rand::Rng;
use rand_chacha::ChaChaRng;
use serde::{Deserialize, Serialize};
use units::time::Time;
use units::{Length, Mass, Temperature};

use crate::StellarColor;

use super::spectral_class::{LuminosityClass, SpectralType, VariabilityType};

// Constants for stellar distribution
const SOLAR_TEMP: f64 = 5778.0; // Kelvin

/// Main sequence stars: core hydrogen burning stars
/// mass: typically 0.08-150 solar masses
/// temperature: 2400K (M) to 50,000K (O)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MainSequenceStar {
    pub mass: Mass,
    pub radius: Length,
    pub luminosity: f64,
    pub temperature: Temperature,
    pub spectral_type: SpectralType,
    pub luminosity_class: LuminosityClass,
    pub subtype: u8,
    pub variability: VariabilityType,
    pub metallicity: f64,
    pub age: Time,
}

impl MainSequenceStar {
    /// Creates a new main sequence star with physically realistic properties based on mass
    /// using mass-luminosity and mass-temperature relations
    ///
    /// # Mass Categories and Properties:
    /// * Very Massive (>30 M☉):
    ///   - Luminosity: 3.0e4 L☉ × (M/30M☉)³·⁵
    ///   - Temperature: 38,000K × (M/30M☉)⁰·²
    ///   - Examples: O-type stars
    ///
    /// * Massive (8-30 M☉):
    ///   - Luminosity: 1.0e3 L☉ × (M/8M☉)³·⁵
    ///   - Temperature: 22,000K × (M/8M☉)⁰·²
    ///   - Examples: B-type stars
    ///
    /// * Intermediate (2-8 M☉):
    ///   - Luminosity: 25 L☉ × (M/2M☉)³·⁵
    ///   - Temperature: 9,000K × (M/2M☉)⁰·²
    ///   - Examples: A and F-type stars
    ///
    /// * Solar-type (0.8-2 M☉):
    ///   - Luminosity: M⁴ L☉
    ///   - Temperature: 5,800K × M⁰·¹
    ///   - Examples: G-type stars like the Sun
    ///
    /// * Low Mass (<0.8 M☉):
    ///   - Luminosity: M²·³ L☉
    ///   - Temperature: 3,500K × (M/0.1M☉)⁰·²
    ///   - Examples: K and M-type stars
    ///
    /// # Arguments
    /// * `mass_solar` - Stellar mass in solar masses
    /// * `metallicity` - Metallicity [Fe/H] in dex (0.0 = solar)
    /// * `age` - Stellar age
    ///
    /// # Returns
    /// Returns a MainSequenceStar struct with:
    /// * Mass-derived luminosity and temperature
    /// * Calculated radius from L-T relation
    /// * Spectral type based on temperature
    /// * Luminosity class V (main sequence)
    ///
    /// # Example
    /// ```
    /// use stellar_forge::stellar::MainSequenceStar;
    /// use units::time::Time;
    ///
    /// let sun_like = MainSequenceStar::new(1.0, 0.0, 1.0);
    /// ```
    pub fn new(mass_solar: f64, metallicity: f64, age: f64) -> Self {
        let (luminosity, temperature) = match mass_solar {
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
            // K dwarfs: 0.45-0.8 M☉, 3700-5200K
            // M dwarfs: 0.08-0.45 M☉, 2400-3700K
            m => {
                let luminosity = m.powf(2.3);
                // Piecewise linear in log-log space to match spectral boundaries:
                // M dwarfs: 0.08 → 2500K, 0.45 → 3700K
                // K dwarfs: 0.45 → 3700K, 0.8 → 4800K
                let temperature = if m < 0.45 {
                    // M dwarf range: steep slope
                    // α = ln(3700/2500) / ln(0.45/0.08) ≈ 0.23
                    2500.0 * (m / 0.08).powf(0.23)
                } else {
                    // K dwarf range: gentler slope
                    // α = ln(4800/3700) / ln(0.8/0.45) ≈ 0.45
                    3700.0 * (m / 0.45).powf(0.45)
                };
                (luminosity, temperature)
            }
        };

        let radius_solar = StellarObject::calculate_radius(luminosity, temperature);
        let spectral_type = StellarObject::spectral_type_from_temp(temperature);
        let subtype = StellarObject::calculate_subtype(temperature);

        Self {
            mass: Mass::from_solar_masses(mass_solar),
            spectral_type,
            luminosity,
            temperature: Temperature::from_kelvin(temperature),
            luminosity_class: LuminosityClass::V,
            radius: Length::from_solar_radii(radius_solar),
            subtype,
            variability: VariabilityType::determine_variability(
                mass_solar,
                temperature,
                LuminosityClass::V,
            ),
            metallicity,
            age: Time::from_myr(age),
        }
    }

    /// Create a solar analog (1 M☉, solar metallicity, young age for disk studies)
    pub fn solar_analog() -> Self {
        Self::new(1.0, 0.0, 1.0)
    }

    /// Sample a random main sequence star from a realistic stellar population
    ///
    /// Uses the Kroupa (2001) IMF for mass distribution, limited to main sequence
    /// stars capable of hosting planetary systems (0.08-3.0 M☉). More massive stars
    /// have short lifetimes and are unlikely to form stable planetary systems.
    ///
    /// # Arguments
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// A MainSequenceStar with realistic mass, metallicity, and age distributions
    ///
    /// # Example
    /// ```
    /// use rand_chacha::ChaChaRng;
    /// use rand::SeedableRng;
    /// use stellar_forge::stellar::MainSequenceStar;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let star = MainSequenceStar::sample(&mut rng);
    /// ```
    pub fn sample(rng: &mut ChaChaRng) -> Self {
        // Sample mass from Kroupa IMF, limited to planet-hosting range
        let mass = Self::sample_mass_kroupa_planet_hosts(rng);

        // Sample metallicity from local galactic distribution
        let metallicity = Self::sample_metallicity(rng);

        // Sample age appropriate for disk formation studies (young to middle-aged)
        let age = Self::sample_age(rng, mass);

        Self::new(mass, metallicity, age.to_myr())
    }

    /// Sample stellar mass from Kroupa IMF, limited to planet-hosting range
    ///
    /// The Kroupa IMF is a broken power law. We limit to 0.08-3.0 M☉ because:
    /// - Below 0.08 M☉: brown dwarfs, not true stars
    /// - Above 3.0 M☉: short lifetimes (<500 Myr), unlikely stable planets
    ///
    /// Distribution heavily favors M dwarfs (~75%), then K dwarfs (~15%),
    /// then G/F stars (~10%).
    fn sample_mass_kroupa_planet_hosts(rng: &mut ChaChaRng) -> f64 {
        // Segment weights for planet-hosting mass range
        // Integrated from Kroupa IMF over each segment
        let segment_weights = [0.75, 0.15, 0.10]; // M, K, G/F stars
        let rand: f64 = rng.random();

        if rand < segment_weights[0] {
            // M dwarfs: 0.08-0.45 M☉, slope -1.3
            Self::sample_power_law_mass(0.08, 0.45, -1.3, rng)
        } else if rand < segment_weights[0] + segment_weights[1] {
            // K dwarfs: 0.45-0.8 M☉, slope -2.3
            Self::sample_power_law_mass(0.45, 0.8, -2.3, rng)
        } else {
            // G/F stars: 0.8-3.0 M☉, slope -2.3
            Self::sample_power_law_mass(0.8, 3.0, -2.3, rng)
        }
    }

    /// Sample from a power-law distribution
    fn sample_power_law_mass(m_min: f64, m_max: f64, alpha: f64, rng: &mut ChaChaRng) -> f64 {
        let x: f64 = rng.random();
        let alpha1 = alpha + 1.0;
        (x * (m_max.powf(alpha1) - m_min.powf(alpha1)) + m_min.powf(alpha1)).powf(1.0 / alpha1)
    }

    /// Sample metallicity from local galactic distribution
    ///
    /// Returns [Fe/H] in dex, centered on solar (0.0) with σ ≈ 0.2 dex
    fn sample_metallicity(rng: &mut ChaChaRng) -> f64 {
        // Box-Muller transform for Gaussian
        let u1: f64 = rng.random();
        let u2: f64 = rng.random();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();

        // Mean = 0.0 (solar), σ = 0.2 dex, clamped to reasonable range
        (z * 0.2).clamp(-0.5, 0.4)
    }

    /// Sample stellar age appropriate for planet formation studies
    ///
    /// Young stars (< 10 Myr) still have protoplanetary disks.
    /// We sample from 1-10 Myr to ensure disk presence for formation.
    fn sample_age(rng: &mut ChaChaRng, _mass_solar: f64) -> Time {
        // For planet formation studies, we need young stars with disks
        // Disks typically last 1-10 Myr, so sample in this range
        // This ensures the population can actually form planets
        let age_myr = 1.0 + rng.random::<f64>() * 9.0; // 1-10 Myr

        Time::from_myr(age_myr)
    }
}

/// Giant stars: evolved stars burning heavier elements
/// mass: typically 0.3-8 solar masses
/// radius: 10-1000 times main sequence radius
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GiantStar {
    pub mass: f64,
    pub radius: f64,
    pub luminosity: f64,
    pub temperature: f64,
    pub spectral_type: SpectralType,
    pub subtype: u8,
    pub luminosity_class: LuminosityClass,
    pub variability: VariabilityType,
}

impl GiantStar {
    /// Creates a new giant star based on initial mass, with physically realistic properties
    /// scaled according to different giant star categories
    ///
    /// # Mass Categories and Properties:
    /// * Hypergiants (>30 M☉):
    ///   - Mass: 78-82% of initial mass
    ///   - Luminosity: ~5e5 L☉ × (M/30M☉)²
    ///   - Temperature: 4,000-30,000K
    ///   - Class: Ia+ (>500,000 L☉) or Ia
    ///
    /// * Supergiants (8-30 M☉):
    ///   - Mass: 85% of initial mass
    ///   - Luminosity: ~5e4 L☉ × (M/8M☉)²
    ///   - Temperature: 3,500-20,000K
    ///   - Class: Ia (>50,000 L☉) or Ib
    ///
    /// * Regular Giants (2-8 M☉):
    ///   - Mass: 95% of initial mass
    ///   - Luminosity: 100 L☉ × (M/2M☉)²
    ///   - Temperature: ~4,500K
    ///   - Class: II (>1,000 L☉) or III
    ///
    /// * Small Giants (0.8-2 M☉):
    ///   - Mass: 98% of initial mass
    ///   - Luminosity: 10 L☉ × M²
    ///   - Temperature: ~4,800K
    ///   - Class: III (>10 L☉) or IV
    ///
    /// # Arguments
    /// * `rng` - Random number generator for stellar properties
    /// * `initial_mass` - Initial mass in solar masses
    ///
    /// # Returns
    /// Returns a GiantStar struct with mass-appropriate properties including:
    /// radius (derived from luminosity and temperature), spectral type,
    /// and luminosity class
    ///
    /// # Example
    /// ```
    /// use rand_chacha::ChaChaRng;
    /// use rand::SeedableRng;
    /// use stellar_forge::stellar::GiantStar;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let supergiant = GiantStar::new(&mut rng, 15.0); // Creates a supergiant
    /// ```
    pub fn new(rng: &mut ChaChaRng, initial_mass: f64) -> Self {
        let (mass, luminosity, temperature, luminosity_class) = match initial_mass {
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

        let radius = StellarObject::calculate_radius(luminosity, temperature);
        let spectral_type = StellarObject::spectral_type_from_temp(temperature);
        let subtype = StellarObject::calculate_subtype(temperature);

        GiantStar {
            mass,
            radius,
            luminosity,
            temperature,
            spectral_type,
            subtype,
            luminosity_class,
            variability: VariabilityType::determine_variability(
                mass,
                temperature,
                LuminosityClass::V,
            ),
        }
    }
}

/// White dwarfs: dense stellar remnants of low/medium mass stars
/// mass: typically 0.17-1.44 solar masses (Chandrasekhar limit)
/// temperature: 4,000K to 150,000K (initially very hot, cooling over time)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WhiteDwarf {
    pub mass: f64,
    pub radius: f64,
    pub luminosity: f64,
    pub spectral_type: WhiteDwarfType,
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
pub enum WhiteDwarfType {
    DA,
    DB,
    DC,
    DQ,
    DZ,
}

impl WhiteDwarf {
    /// Creates a new white dwarf with physically accurate properties based on mass-radius-luminosity relations
    ///
    /// # Properties Generated:
    /// * Mass: 0.17-1.44 solar masses (up to Chandrasekhar limit)
    /// * Radius: Calculated using mass-radius relation R ∝ M^(-1/3)
    /// * Luminosity: Scaled from typical middle-aged WD (0.001 solar luminosities)
    /// * Spectral Type: Determined from luminosity
    ///
    /// # Physical Relations:
    /// * Radius relation assumes electron-degenerate matter
    /// * Base luminosity of 0.001 L☉ represents typical cooling white dwarf
    /// * Upper mass limit is the Chandrasekhar limit (1.44 M☉)
    ///
    /// # Arguments
    /// * `rng` - Random number generator for white dwarf properties
    ///
    /// # Returns
    /// Returns a WhiteDwarf struct with correlated physical properties
    ///
    /// # Example
    /// ```
    /// use rand_chacha::ChaChaRng;
    /// use rand::SeedableRng;
    /// use stellar_forge::stellar::WhiteDwarf;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let white_dwarf = WhiteDwarf::new(&mut rng);
    /// ```
    pub fn new(rng: &mut ChaChaRng) -> Self {
        let mass: f64 = rng.random_range(0.17..=1.44);
        let radius = 0.01 * (0.6 / mass).powf(1.0 / 3.0);
        let base_luminosity = 0.001; // Typical middle-aged WD
        let luminosity = base_luminosity * (mass / 0.6).powf(1.0);
        let spectral_type = Self::determine_spectral_class(luminosity);

        Self {
            mass,
            radius,
            luminosity,
            spectral_type,
        }
    }

    fn determine_spectral_class(luminosity: f64) -> WhiteDwarfType {
        // White dwarf spectral classification based on temperature/luminosity
        match luminosity {
            l if l > 0.01 => WhiteDwarfType::DA,
            l if l > 0.001 => WhiteDwarfType::DB,
            l if l > 0.0001 => WhiteDwarfType::DC,
            l if l > 0.00001 => WhiteDwarfType::DQ,
            _ => WhiteDwarfType::DZ,
        }
    }
}

/// Neutron stars: ultra-dense stellar remnants
/// mass: 1.1-2.5 solar masses
/// radius: ~12km (remarkably consistent)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NeutronStar {
    pub mass: f64,           // Solar masses
    pub radius: f64,         // km
    pub magnetic_field: f64, // Log Gauss
    pub pulsar: bool,
    pub magnetar: bool,
}

impl NeutronStar {
    /// Creates a new neutron star with randomized but physically realistic properties
    ///
    /// # Properties Generated:
    /// * Mass: 1.4-2.0 solar masses (typical observed range)
    /// * Radius: 11-13 km (typical neutron star size)
    /// * Magnetic Field: Log(B) of 11-13 Gauss (typical field strength)
    /// * Type Classification:
    ///   - Magnetar: If magnetic field > 14 Log(B)
    ///   - Pulsar: 80% chance if not a magnetar
    ///
    /// # Arguments
    /// * `rng` - Random number generator for neutron star properties
    ///
    /// # Returns
    /// Returns a NeutronStar struct with randomized physical properties
    /// within observed astronomical ranges
    ///
    /// # Example
    /// ```
    /// use rand_chacha::ChaChaRng;
    /// use rand::SeedableRng;
    /// use stellar_forge::stellar::NeutronStar;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let neutron_star = NeutronStar::new(&mut rng);
    /// ```
    pub fn new(rng: &mut ChaChaRng) -> Self {
        // Mass distribution following Özel & Freire 2016
        // Using a peaked distribution around 1.4 M☉
        let mass = Self::sample_mass(rng);

        // Radius based on mass-radius relationship
        let radius = Self::calculate_radius(mass);

        // Magnetic field distribution is bimodal
        // Regular NS: 11-13, Magnetars: 14-15
        let magnetic_field = match rng.random_bool(0.1) {
            // ~10% magnetar fraction
            true => rng.random_range(14.0..=15.0),  // Magnetar
            false => rng.random_range(11.0..=13.0), // Regular NS
        };

        let magnetar = magnetic_field >= 14.0;

        // Pulsars are more common among regular NS
        let pulsar = match magnetar {
            true => rng.random_bool(0.3),  // 30% of magnetars are pulsars
            false => rng.random_bool(0.8), // 80% of regular NS are pulsars
        };

        Self {
            mass,
            radius,
            magnetic_field,
            pulsar,
            magnetar,
        }
    }

    fn sample_mass(rng: &mut ChaChaRng) -> f64 {
        // Box-Muller transform to generate normal distribution
        let u1: f64 = rng.random();
        let u2: f64 = rng.random();

        let z: f64 = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();

        // Parameters for the distribution
        let mean = 1.4;
        let std_dev = 0.2;
        let skew_factor = 0.3;

        // Apply skew transformation and constraints
        let skewed = mean + std_dev * z + skew_factor * z.abs();

        // Clamp to physical limits
        skewed.clamp(1.1, 2.5)
    }

    fn calculate_radius(mass: f64) -> f64 {
        // Approximation of mass-radius relationship
        // Based on various EOS models
        // Returns radius in km
        12.0 * (mass / 1.4).powf(-0.3)
    }
}

/// Black holes: objects with gravity so strong that nothing can escape
/// mass: typically 3-100 solar masses (stellar), millions to billions (supermassive)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BlackHole {
    pub mass: f64,
    pub spin: f64,
    pub has_accretion: bool,
}

impl BlackHole {
    /// Creates a new black hole based on progenitor star parameters
    ///
    /// # Arguments
    /// * `rng` - Random number generator for black hole properties
    /// * `initial_mass` - Mass of the progenitor star in solar masses
    /// * `metallicity` - Metallicity of the progenitor star
    ///
    /// # Returns
    /// Returns a BlackHole struct with:
    /// * Calculated final mass based on progenitor properties
    /// * Randomized spin parameter (0 to 1)
    /// * Random chance of having an accretion disk (10% probability)
    ///
    /// # Example
    /// ```
    /// use rand_chacha::ChaChaRng;
    /// use rand::SeedableRng;
    /// use stellar_forge::stellar::BlackHole;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let black_hole = BlackHole::new(&mut rng, 25.0, 0.02); // 25 solar mass progenitor, solar metallicity
    /// ```
    pub fn new(rng: &mut ChaChaRng, initial_mass: f64, metallicity: f64) -> Self {
        let mass = Self::calculate_mass(rng, initial_mass, metallicity);
        let spin = Self::calculate_spin(rng, initial_mass);
        let has_accretion = rng.random_bool(0.1); // 10% chance

        Self {
            mass,
            has_accretion,
            spin,
        }
    }

    fn calculate_mass(rng: &mut ChaChaRng, initial_mass: f64, metallicity: f64) -> f64 {
        // Mass loss due to winds scales with metallicity
        let wind_loss = initial_mass * metallicity * 0.4;

        // Additional mass loss in supernova
        let sn_loss = initial_mass * rng.random_range(0.2..0.5);

        // Final BH mass
        let final_mass = initial_mass - wind_loss - sn_loss;

        // Add some random variation
        let variation = rng.random_range(0.8..1.2);

        final_mass * variation
    }

    fn calculate_spin(rng: &mut ChaChaRng, initial_mass: f64) -> f64 {
        // More massive stars tend to produce higher spin BHs
        let base_spin = (initial_mass - 20.0) / 100.0;
        let variation = rng.random_range(-0.2..0.2);

        // Clamp to physical range [0,1]
        (base_spin + variation).clamp(0.0, 1.0)
    }

    // fn calculate_schwarzschild_radius(mass: f64) -> f64 {
    //     // R_s = 2GM/c^2 in solar radii
    //     2.95 * mass
    // }

    // fn calculate_isco_radius(mass: f64, spin: f64) -> f64 {
    //     // Approximate ISCO calculation including spin effects
    //     let z1 = 1.0
    //         + (1.0 - spin * spin).powf(1.0 / 3.0)
    //             * ((1.0 + spin).powf(1.0 / 3.0) + (1.0 - spin).powf(1.0 / 3.0));
    //     let z2 = (3.0 * spin * spin + z1 * z1).sqrt();

    //     mass * (3.0 + z2 - ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt())
    // }

    // fn calculate_hawking_temperature(mass: f64) -> f64 {
    //     // T = ℏc³/8πGMk_B in Kelvin
    //     6.169e-8 / mass
    // }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum StellarObject {
    MainSequence(MainSequenceStar),
    Giant(GiantStar),
    NeutronStar(NeutronStar),
    WhiteDwarf(WhiteDwarf),
    BlackHole(BlackHole),
}

impl StellarObject {
    /// Creates a new star based on its initial parameters
    ///
    /// # Arguments
    /// * `rng` - Random number generator for stellar calculations
    /// * `initial_mass` - The star's initial mass in solar masses
    /// * `age` - The current age of the star in years
    /// * `metallicity` - The star's metal content (Z) as a fraction
    ///
    /// # Returns
    /// Returns either a living star or stellar remnant depending on the age:
    /// * If age < lifetime: Returns a living star
    /// * If age > lifetime: Returns a remnant (white dwarf, neutron star, or black hole)
    ///
    /// # Example
    /// ```
    /// use rand_chacha::ChaChaRng;
    /// use rand::SeedableRng;
    /// use stellar_forge::stellar::StellarObject;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let star = StellarObject::new(&mut rng, 1.0, 1.0e9, 0.02); // 1 solar mass, 1 billion years old, solar metallicity
    /// ```
    pub fn new(rng: &mut ChaChaRng, initial_mass: f64, age_years: f64, metallicity: f64) -> Self {
        // estimate_lifetime returns Gyr, convert age_years to Gyr for comparison
        let lifetime_gyr = Self::estimate_lifetime(initial_mass);
        let age_gyr = age_years / 1.0e9;
        let life_fraction = age_gyr / lifetime_gyr;

        if life_fraction > 1.0 {
            return Self::new_remnant_star(rng, initial_mass, metallicity);
        }

        let age = Time::from_years(age_years);
        Self::new_living_star(rng, life_fraction, initial_mass, metallicity, age.to_myr())
    }

    fn new_living_star(
        rng: &mut ChaChaRng,
        life_fraction: f64,
        initial_mass: f64,
        metallicity: f64,
        age: f64,
    ) -> Self {
        match initial_mass {
            // Very massive stars (> 30 M☉)
            m if m > 30.0 => {
                if life_fraction < 0.9 {
                    StellarObject::MainSequence(MainSequenceStar::new(m, metallicity, age))
                } else {
                    StellarObject::Giant(GiantStar::new(rng, m))
                }
            }
            // Massive stars (8-30 M☉)
            m if m > 8.0 => {
                if life_fraction < 0.9 {
                    StellarObject::MainSequence(MainSequenceStar::new(m, metallicity, age))
                } else {
                    StellarObject::Giant(GiantStar::new(rng, m))
                }
            }
            // Intermediate mass stars (2-8 M☉)
            m if m > 2.0 => match life_fraction {
                x if x < 0.85 => {
                    StellarObject::MainSequence(MainSequenceStar::new(m, metallicity, age))
                }
                x if x < 0.95 => StellarObject::Giant(GiantStar::new(rng, m)),
                _ => StellarObject::WhiteDwarf(WhiteDwarf::new(rng)),
            },
            // Solar-type stars (0.8-2 M☉)
            m if m > 0.8 => match life_fraction {
                x if x < 0.90 => {
                    StellarObject::MainSequence(MainSequenceStar::new(m, metallicity, age))
                }
                x if x < 0.95 => StellarObject::Giant(GiantStar::new(rng, m)),
                _ => StellarObject::WhiteDwarf(WhiteDwarf::new(rng)),
            },
            // Low mass stars (< 0.8 M☉) - always main sequence
            m => StellarObject::MainSequence(MainSequenceStar::new(m, metallicity, age)),
        }
    }

    fn new_remnant_star(rng: &mut ChaChaRng, initial_mass: f64, metallicity: f64) -> Self {
        match initial_mass {
            // // Main sequence stars below 0.08 solar masses never form
            // m if m < 0.08 => StellarObject::BrownDwarf(BrownDwarf::new(initial_mass, temperature)),

            // // White dwarf range: 0.08-8 solar masses
            m if m < 8.0 => StellarObject::WhiteDwarf(WhiteDwarf::new(rng)),

            // Neutron star range: 8-20 solar masses
            m if m < 20.0 => StellarObject::NeutronStar(NeutronStar::new(rng)),

            // Black hole range: >20 solar masses
            _ => StellarObject::BlackHole(BlackHole::new(rng, initial_mass, metallicity)),
        }
    }

    pub fn estimate_lifetime(mass: f64) -> f64 {
        // Lifetime in billions of years
        match mass {
            m if m > 10.0 => 0.0001 * (m / 10.0).powf(-2.0),
            m if m > 2.0 => 0.01 * (m / 2.0).powf(-2.5),
            m if m > 0.5 => 10.0 * (m).powf(-3.0),
            _ => 100.0 * (mass).powf(-2.5), // Low mass stars live very long
        }
    }

    /// Sample a stellar object from realistic galactic distributions
    ///
    /// Draws mass from the Kroupa IMF, age from a thin disk distribution,
    /// and metallicity from a Gaussian centered on solar.
    ///
    /// # Arguments
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// A StellarObject appropriate for the sampled parameters
    ///
    /// # Example
    /// ```
    /// use rand::SeedableRng;
    /// use rand_chacha::ChaChaRng;
    /// use stellar_forge::stellar::StellarObject;
    ///
    /// let mut rng = ChaChaRng::seed_from_u64(42);
    /// let star = StellarObject::sample(&mut rng);
    /// ```
    pub fn sample(rng: &mut ChaChaRng) -> Self {
        let mass = Self::sample_mass_kroupa(rng);
        let metallicity = Self::sample_metallicity(rng);

        // Sample age, but cap at stellar lifetime to avoid all stars becoming remnants
        let max_lifetime_gyr = Self::estimate_lifetime(mass);
        let age_gyr = Self::sample_age_capped(rng, max_lifetime_gyr);

        let age_years = age_gyr * 1.0e9;
        Self::new(rng, mass, age_years, metallicity)
    }

    /// Sample stellar mass from Kroupa (2001) IMF
    ///
    /// The Kroupa IMF is a broken power law:
    /// - 0.08 ≤ M < 0.5: α = -1.3
    /// - 0.5 ≤ M < 1.0: α = -2.3
    /// - M ≥ 1.0: α = -2.3
    ///
    /// Returns mass in solar masses, range 0.08 - 150 M☉
    fn sample_mass_kroupa(rng: &mut ChaChaRng) -> f64 {
        // Segment weights from integrating each segment of the IMF
        let segment_weights = [0.80, 0.15, 0.05];
        let rand: f64 = rng.random();

        if rand < segment_weights[0] {
            // 0.08-0.5 M☉: slope -1.3
            Self::sample_power_law(0.08, 0.5, -1.3, rng)
        } else if rand < segment_weights[0] + segment_weights[1] {
            // 0.5-1.0 M☉: slope -2.3
            Self::sample_power_law(0.5, 1.0, -2.3, rng)
        } else {
            // 1.0-150.0 M☉: slope -2.3
            Self::sample_power_law(1.0, 150.0, -2.3, rng)
        }
    }

    /// Sample from a power-law distribution between m_min and m_max with given slope
    fn sample_power_law(m_min: f64, m_max: f64, alpha: f64, rng: &mut ChaChaRng) -> f64 {
        let x: f64 = rng.random();
        let alpha1 = alpha + 1.0;

        (x * (m_max.powf(alpha1) - m_min.powf(alpha1)) + m_min.powf(alpha1)).powf(1.0 / alpha1)
    }

    /// Sample stellar age from galactic thin disk distribution
    ///
    /// Samples age uniformly from the thin disk (0.1-10 Gyr), but caps
    /// at the stellar lifetime to ensure we get a mix of main sequence
    /// stars and evolved objects.
    ///
    /// # Arguments
    /// * `rng` - Random number generator
    /// * `max_lifetime_gyr` - Maximum stellar lifetime in Gyr (from estimate_lifetime)
    ///
    /// # Returns
    /// Age in Gyr, capped at min(10 Gyr, max_lifetime_gyr)
    fn sample_age_capped(rng: &mut ChaChaRng, max_lifetime_gyr: f64) -> f64 {
        const THIN_DISK_MAX_AGE: f64 = 10.0;
        const MIN_AGE: f64 = 0.1;

        // Cap the max age at either the disk age or stellar lifetime
        let effective_max = max_lifetime_gyr.min(THIN_DISK_MAX_AGE);

        // Ensure we have a valid range
        if effective_max <= MIN_AGE {
            // Star lives less than 100 Myr - give it a very young age
            return MIN_AGE * rng.random::<f64>();
        }

        // Uniform distribution from MIN_AGE to effective_max
        MIN_AGE + rng.random::<f64>() * (effective_max - MIN_AGE)
    }

    /// Sample metallicity from local galactic distribution
    ///
    /// Returns [Fe/H] in dex, centered on solar (0.0)
    /// Uses Gaussian with σ ≈ 0.2 dex
    fn sample_metallicity(rng: &mut ChaChaRng) -> f64 {
        // Box-Muller transform for Gaussian
        let u1: f64 = rng.random();
        let u2: f64 = rng.random();

        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();

        // Mean = 0.0 (solar), σ = 0.2 dex
        let metallicity = z * 0.2;

        // Clamp to reasonable range
        metallicity.clamp(-0.5, 0.4)
    }

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

    pub fn spectral_type_from_temp(temperature: f64) -> SpectralType {
        for (spec_type, temp_bound) in Self::TEMP_BOUNDS.iter() {
            if temperature >= *temp_bound {
                return *spec_type;
            }
        }
        SpectralType::Y
    }

    pub fn calculate_radius(luminosity: f64, temperature: f64) -> f64 {
        // Using L = 4πR²σT⁴, solved for R in solar units
        // Since L is in solar units and T is in Kelvin, we need to compare with solar temperature
        luminosity.powf(0.5) / (temperature / SOLAR_TEMP).powf(2.0)
    }

    pub fn calculate_subtype(temperature: f64) -> u8 {
        // Find the relevant temperature boundaries
        let mut upper_bound = Self::TEMP_BOUNDS[0].1;
        let mut lower_bound = Self::TEMP_BOUNDS[1].1;

        for window in Self::TEMP_BOUNDS.windows(2) {
            if temperature >= window[1].1 {
                upper_bound = window[0].1;
                lower_bound = window[1].1;
                break;
            }
        }

        // Calculate where this temperature falls within the range
        let temp_range = upper_bound - lower_bound;
        let temp_position = upper_bound - temperature;

        // Convert to 0-9 subtype scale
        let subtype = (9.0 * temp_position / temp_range).round() as u8;

        // Clamp to valid range
        subtype.clamp(0, 9)
    }

    /// Returns the apparent color of this stellar object as RGB values (0-255)
    ///
    /// Colors are based on blackbody radiation for temperature-based objects,
    /// with special handling for exotic objects like black holes.
    pub fn color(&self) -> StellarColor {
        match self {
            StellarObject::MainSequence(star) => {
                StellarColor::from_temperature(star.temperature.to_kelvin())
            }
            StellarObject::Giant(star) => StellarColor::from_temperature(star.temperature),
            StellarObject::WhiteDwarf(wd) => {
                // White dwarfs range from very hot (blue-white) to cool (red)
                // Use spectral type to estimate temperature
                let temp = match wd.spectral_type {
                    WhiteDwarfType::DA | WhiteDwarfType::DB => 15000.0, // Hot, blue-white
                    WhiteDwarfType::DC => 8000.0,                       // Cooler, white
                    WhiteDwarfType::DQ | WhiteDwarfType::DZ => 6000.0,  // Even cooler
                };
                StellarColor::from_temperature(temp)
            }
            StellarObject::NeutronStar(ns) => {
                // Neutron stars: young ones are hot (blue), old ones cool
                // Magnetars and pulsars tend to be more energetic
                if ns.magnetar {
                    StellarColor::new(200, 220, 255) // Bright blue-white
                } else if ns.pulsar {
                    StellarColor::new(180, 200, 255) // Blue-white
                } else {
                    StellarColor::new(160, 180, 220) // Dimmer blue
                }
            }
            StellarObject::BlackHole(bh) => {
                // Black holes: dark with possible accretion disk glow
                if bh.has_accretion {
                    StellarColor::new(255, 200, 150) // Orange-white accretion glow
                } else {
                    StellarColor::new(20, 20, 30) // Nearly black with faint blue
                }
            }
        }
    }
}
