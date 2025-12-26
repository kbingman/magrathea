//! Moon (satellite) system generation
//!
//! Generates moons for planets based on statistical occurrence rates and
//! physical constraints (Hill sphere, Roche limit, tidal locking).
//!
//! # Formation Mechanisms
//!
//! - **Giant Impact**: Large collision creates debris disk (Earth-Moon style)
//! - **Capture**: Gravitational capture of passing body (Triton style)
//! - **Co-accretion**: Formation in circumplanetary disk (Galilean style)
//!
//! # References
//!
//! - Canup (2004) - "Simulations of a late lunar-forming impact"
//! - Agnor & Hamilton (2006) - "Neptune's capture of its moon Triton"
//! - Sasaki et al. (2010) - "Origin of the regular satellites of giant planets"
//! - Heller & Barnes (2013) - "Exomoon habitability constrained by tidal heating"

use planetary::Planet;
use planetary::composition::Composition;
use planetary::moon::{Moon, MoonFormation, MoonSystem, moon_numeral};
use planetary::planet_class::PlanetClass;
use planetary::planet_type::PlanetType;
use rand::Rng;
use rand_chacha::ChaChaRng;
use units::{Length, Mass};

// =============================================================================
// Physical Constants
// =============================================================================

/// Seconds per year
const SECONDS_PER_YEAR: f64 = 3.154e7;

// =============================================================================
// Occurrence Rates by Planet Type
// =============================================================================

/// Moon count distribution for terrestrial planets (< 2 M⊕)
/// [none, one, two, three+]
const TERRESTRIAL_MOON_DIST: [f64; 4] = [0.70, 0.25, 0.05, 0.00];

/// Moon count distribution for super-Earths (2-10 M⊕)
const SUPER_EARTH_MOON_DIST: [f64; 4] = [0.60, 0.28, 0.10, 0.02];

/// Moon count distribution for ice giants (10-50 M⊕)
const ICE_GIANT_MOON_DIST: [f64; 4] = [0.15, 0.25, 0.35, 0.25];

/// Moon count distribution for gas giants (> 50 M⊕)
pub const GAS_GIANT_MOON_DIST: [f64; 4] = [0.03, 0.12, 0.35, 0.50];

/// Ring probability by planet class
const RING_PROBABILITY: [f64; 4] = [
    0.00, // Terrestrial
    0.01, // Super-Earth
    0.50, // Ice Giant
    0.80, // Gas Giant
];

// =============================================================================
// Hill Sphere and Roche Limit
// =============================================================================

/// Calculate Hill sphere radius in meters
///
/// The Hill sphere defines the region where the planet's gravity dominates
/// over the star's gravity. Moons must orbit well within this sphere
/// (~30% of Hill radius) for long-term stability.
///
/// R_Hill = a × (M_planet / (3 × M_star))^(1/3)
///
/// # Arguments
/// * `semi_major_axis_m` - Planet's orbital distance from star in meters
/// * `planet_mass_kg` - Planet mass in kg
/// * `stellar_mass_kg` - Star mass in kg
#[inline]
pub fn hill_radius(semi_major_axis_m: f64, planet_mass_kg: f64, stellar_mass_kg: f64) -> f64 {
    semi_major_axis_m * (planet_mass_kg / (3.0 * stellar_mass_kg)).powf(1.0 / 3.0)
}

/// Calculate Roche limit in meters
///
/// Inside the Roche limit, tidal forces exceed the moon's self-gravity,
/// causing disruption into a ring system.
///
/// R_Roche ≈ 2.44 × R_planet × (ρ_planet / ρ_moon)^(1/3)
///
/// # Arguments
/// * `planet_radius_m` - Planet radius in meters
/// * `planet_density` - Planet bulk density in kg/m³
/// * `moon_density` - Moon bulk density in kg/m³ (default ~2000 for icy, ~3000 for rocky)
pub fn roche_limit(planet_radius_m: f64, planet_density: f64, moon_density: f64) -> f64 {
    2.44 * planet_radius_m * (planet_density / moon_density).powf(1.0 / 3.0)
}

/// Calculate tidal locking timescale in years
///
/// Most close-in moons become tidally locked to their host planet.
///
/// τ_lock ∝ a^6 / (M_planet × R_moon^2 × k₂/Q)
///
/// Simplified formula assuming typical dissipation factor.
fn tidal_locking_timescale(
    semi_major_axis_m: f64,
    planet_mass_kg: f64,
    moon_radius_m: f64,
    moon_mass_kg: f64,
) -> f64 {
    // Simplified tidal locking formula
    // Using Q/k2 ~ 100 for rocky moons (Earth-like dissipation)
    let q_over_k2 = 100.0;
    let a6 = semi_major_axis_m.powi(6);
    let denominator = planet_mass_kg * moon_radius_m.powi(2) * moon_mass_kg;

    // Rough scaling to get years
    (a6 / denominator) * q_over_k2 * 1e-20 / SECONDS_PER_YEAR
}

/// Calculate tidal heat flux in W/m²
///
/// Tidal heating arises from gravitational flexing as the moon's orbit
/// is perturbed. Higher eccentricity and closer orbits produce more heating.
///
/// Q_tidal ∝ (M_planet × R_moon^5 × e^2) / (a^6 × Q)
fn tidal_heat_flux(
    planet_mass_kg: f64,
    moon_radius_m: f64,
    semi_major_axis_m: f64,
    eccentricity: f64,
) -> f64 {
    // Simplified tidal heating formula
    // Io's heat flux is ~2 W/m², Europa's is ~0.1 W/m²
    let q_factor = 100.0; // Tidal dissipation factor

    let numerator = planet_mass_kg * moon_radius_m.powi(5) * eccentricity.powi(2);
    let denominator = semi_major_axis_m.powi(6) * q_factor;

    // Scale factor calibrated to give ~2 W/m² for Io-like parameters
    (numerator / denominator) * 1e-10
}

// =============================================================================
// Moon Generation
// =============================================================================

/// Generate a moon system for a planet
///
/// This is the main entry point for moon generation. It:
/// 1. Determines if the planet should have moons based on its class
/// 2. Samples the number of moons from class-specific distributions
/// 3. Generates each moon with appropriate mass, orbit, and properties
/// 4. Determines if the planet has a ring system
///
/// # Arguments
/// * `rng` - Random number generator
/// * `planet` - The host planet
/// * `stellar_mass` - Host star mass in solar masses
/// * `planet_id` - Planet identifier for moon naming
/// * `catalog_name` - System catalog name for moon naming
pub fn generate_moon_system(
    rng: &mut ChaChaRng,
    planet: &Planet,
    stellar_mass: f64,
    planet_id: &str,
    catalog_name: &str,
) -> Option<MoonSystem> {
    let planet_mass_earth = planet.mass.to_earth_masses();

    // Very small bodies don't have major moons
    if planet_mass_earth < 0.05 {
        return None;
    }

    // Determine planet category for occurrence rates
    let category = planet_category(planet_mass_earth);
    let moon_dist = moon_distribution(category);
    let ring_prob = RING_PROBABILITY[category];

    // Sample moon count
    let n_moons = sample_moon_count(rng, moon_dist);

    // Check for rings (gas/ice giants only)
    let has_rings = category >= 2 && rng.random::<f64>() < ring_prob;

    // If no moons and no rings, return None
    if n_moons == 0 && !has_rings {
        return None;
    }

    // Generate moons
    let moons = if n_moons > 0 {
        generate_moons(
            rng,
            planet,
            stellar_mass,
            n_moons,
            category,
            planet_id,
            catalog_name,
        )
    } else {
        Vec::new()
    };

    Some(MoonSystem::with_moons(moons, has_rings))
}

/// Categorize planet by mass for moon occurrence rates
/// Returns: 0=Terrestrial, 1=Super-Earth, 2=Ice Giant, 3=Gas Giant
pub fn planet_category(mass_earth: f64) -> usize {
    match mass_earth {
        m if m < 2.0 => 0,  // Terrestrial
        m if m < 10.0 => 1, // Super-Earth
        m if m < 50.0 => 2, // Ice Giant
        _ => 3,             // Gas Giant
    }
}

/// Get moon count distribution for planet category
fn moon_distribution(category: usize) -> [f64; 4] {
    match category {
        0 => TERRESTRIAL_MOON_DIST,
        1 => SUPER_EARTH_MOON_DIST,
        2 => ICE_GIANT_MOON_DIST,
        _ => GAS_GIANT_MOON_DIST,
    }
}

/// Sample moon count from cumulative distribution
pub fn sample_moon_count(rng: &mut ChaChaRng, dist: [f64; 4]) -> usize {
    let roll: f64 = rng.random();
    let mut cumulative = 0.0;

    for (i, prob) in dist.iter().enumerate() {
        cumulative += prob;
        if roll < cumulative {
            return match i {
                0 => 0,
                1 => 1,
                2 => 2,
                3 => rng.random_range(3..=6), // 3+ moons: sample 3-6
                _ => 0,
            };
        }
    }
    0
}

/// Generate moons for a planet
fn generate_moons(
    rng: &mut ChaChaRng,
    planet: &Planet,
    stellar_mass: f64,
    n_moons: usize,
    category: usize,
    planet_id: &str,
    catalog_name: &str,
) -> Vec<Moon> {
    let planet_mass_kg = planet.mass.to_kg();
    let planet_radius_m = planet.radius.to_m();
    let stellar_mass_kg = Mass::from_solar_masses(stellar_mass).to_kg();
    let planet_sma_m = planet.semi_major_axis.to_m();

    // Calculate Hill sphere and Roche limit
    let r_hill = hill_radius(planet_sma_m, planet_mass_kg, stellar_mass_kg);
    let planet_density = planet.density() * 1000.0; // Convert g/cm³ to kg/m³

    // Stable moon zone: from Roche limit to ~30% of Hill sphere
    let r_roche = roche_limit(planet_radius_m, planet_density, 2500.0); // Assume rocky moon
    let r_stable = r_hill * 0.30;

    // If stable zone is too small, no moons possible
    if r_stable <= r_roche * 2.0 {
        return Vec::new();
    }

    let mut moons = Vec::with_capacity(n_moons);

    for i in 0..n_moons {
        let moon = generate_single_moon(
            rng,
            planet,
            stellar_mass,
            category,
            i,
            n_moons,
            r_roche,
            r_stable,
            planet_id,
            catalog_name,
        );
        moons.push(moon);
    }

    // Sort by semi-major axis
    moons.sort_by(|a, b| {
        a.semi_major_axis
            .to_m()
            .partial_cmp(&b.semi_major_axis.to_m())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Reassign names after sorting
    moons.iter_mut().enumerate().for_each(|(i, moon)| {
        moon.id = format!("{}-{}", planet_id, moon_numeral(i));
        moon.name = format!("{} {}", catalog_name, moon_numeral(i));
    });

    moons
}

/// Generate a single moon
#[allow(clippy::too_many_arguments)]
fn generate_single_moon(
    rng: &mut ChaChaRng,
    planet: &Planet,
    stellar_mass: f64,
    category: usize,
    index: usize,
    total_moons: usize,
    r_roche: f64,
    r_stable: f64,
    planet_id: &str,
    catalog_name: &str,
) -> Moon {
    let planet_mass_kg = planet.mass.to_kg();
    let planet_mass_earth = planet.mass.to_earth_masses();

    // Sample formation mechanism
    let formation = sample_formation(rng, category);

    // Sample mass ratio based on formation and planet type
    let mass_ratio = sample_mass_ratio(rng, category, &formation);
    let moon_mass_earth = planet_mass_earth * mass_ratio;
    let moon_mass = Mass::from_earth_masses(moon_mass_earth);

    // Calculate moon radius (simplified mass-radius relation for rocky/icy bodies)
    let moon_radius_earth = if moon_mass_earth > 0.01 {
        // Larger moons: use rocky planet scaling
        moon_mass_earth.powf(0.27)
    } else {
        // Small moons: approximate spherical
        moon_mass_earth.powf(0.33)
    };
    let moon_radius = Length::from_earth_radii(moon_radius_earth);

    // Sample orbital distance
    // Spread moons logarithmically between Roche limit and stable boundary
    let log_roche = (r_roche * 2.0).ln();
    let log_stable = r_stable.ln();
    let log_spacing = (log_stable - log_roche) / (total_moons as f64 + 1.0);
    let log_base = log_roche + log_spacing * (index as f64 + 1.0);
    let scatter = rng.random_range(-0.2..0.2) * log_spacing;
    let semi_major_axis_m = (log_base + scatter).exp();
    let semi_major_axis = Length::from_meters(semi_major_axis_m);

    // Sample eccentricity (lower for co-accretion, higher for capture)
    let eccentricity = match &formation {
        MoonFormation::CoAccretion { .. } => rng.random_range(0.001..0.02),
        MoonFormation::GiantImpact => rng.random_range(0.01..0.05),
        MoonFormation::Capture { .. } => rng.random_range(0.1..0.4),
    };

    // Calculate tidal heating
    let moon_radius_m = moon_radius.to_m();
    let heat_flux = tidal_heat_flux(
        planet_mass_kg,
        moon_radius_m,
        semi_major_axis_m,
        eccentricity,
    );

    // Determine if tidally locked
    let moon_mass_kg = moon_mass.to_kg();
    let lock_timescale = tidal_locking_timescale(
        semi_major_axis_m,
        planet_mass_kg,
        moon_radius_m,
        moon_mass_kg,
    );
    let tidally_locked = lock_timescale < 1.0e9; // Locked if timescale < 1 Gyr

    // Calculate surface temperature
    // Includes stellar heating (from planet's equilibrium temp) + tidal heating
    let stellar_temp = planet.equilibrium_temp;
    let tidal_temp_contribution = (heat_flux * 100.0).powf(0.25); // Rough conversion
    let surface_temp = (stellar_temp.powi(4) + tidal_temp_contribution.powi(4)).powf(0.25);

    // Determine moon type based on environment
    let class = PlanetClass::from_earth_masses(moon_mass_earth);
    let composition = sample_moon_composition(rng, category, heat_flux);

    // For moon classification, use planet's stellar context
    let incident_flux = 1.0 / planet.semi_major_axis.to_au().powi(2);

    let moon_type = if heat_flux > 1.0 {
        // Extreme tidal heating - lava world
        PlanetType::Lava {
            tidally_locked,
            surface_temp_k: surface_temp,
        }
    } else if heat_flux > 0.1 && composition.water > 0.2 {
        // Moderate tidal heating with ice - subsurface ocean likely
        PlanetType::Oceanic {
            ocean_fraction: 0.8 + rng.random::<f64>() * 0.2,
        }
    } else if surface_temp < 150.0 {
        // Cold moon
        PlanetType::Frozen {
            has_subsurface_ocean: heat_flux > 0.05,
        }
    } else {
        // Use standard classification
        PlanetType::from_environment(
            class,
            &composition,
            surface_temp,
            incident_flux,
            moon_mass_earth,
            stellar_mass,
            planet.semi_major_axis.to_au(),
        )
    };

    Moon::new(
        format!("{}-{}", planet_id, moon_numeral(index)),
        format!("{} {}", catalog_name, moon_numeral(index)),
        moon_mass,
        moon_radius,
        semi_major_axis,
        eccentricity,
        class,
        moon_type,
        formation,
        surface_temp,
        heat_flux,
        tidally_locked,
    )
}

/// Sample formation mechanism based on planet category
fn sample_formation(rng: &mut ChaChaRng, category: usize) -> MoonFormation {
    let roll: f64 = rng.random();

    match category {
        0 => {
            // Terrestrial: 60% impact, 30% capture, 10% co-accretion
            match roll {
                r if r < 0.60 => MoonFormation::GiantImpact,
                r if r < 0.90 => MoonFormation::Capture {
                    retrograde: rng.random::<f64>() < 0.3,
                },
                _ => MoonFormation::CoAccretion {
                    resonance_order: None,
                },
            }
        }
        1 => {
            // Super-Earth: 30% impact, 25% capture, 45% co-accretion
            match roll {
                r if r < 0.30 => MoonFormation::GiantImpact,
                r if r < 0.55 => MoonFormation::Capture {
                    retrograde: rng.random::<f64>() < 0.3,
                },
                _ => MoonFormation::CoAccretion {
                    resonance_order: sample_resonance(rng),
                },
            }
        }
        2 => {
            // Ice Giant: 10% impact, 40% capture, 50% co-accretion
            match roll {
                r if r < 0.10 => MoonFormation::GiantImpact,
                r if r < 0.50 => MoonFormation::Capture {
                    retrograde: rng.random::<f64>() < 0.4, // Higher retrograde chance (Triton)
                },
                _ => MoonFormation::CoAccretion {
                    resonance_order: sample_resonance(rng),
                },
            }
        }
        _ => {
            // Gas Giant: 5% impact, 35% capture, 60% co-accretion
            match roll {
                r if r < 0.05 => MoonFormation::GiantImpact,
                r if r < 0.40 => MoonFormation::Capture {
                    retrograde: rng.random::<f64>() < 0.3,
                },
                _ => MoonFormation::CoAccretion {
                    resonance_order: sample_resonance(rng),
                },
            }
        }
    }
}

/// Sample resonance order for co-accreted moons
fn sample_resonance(rng: &mut ChaChaRng) -> Option<(u8, u8)> {
    // ~40% of co-accreted moons are in resonance
    if rng.random::<f64>() < 0.4 {
        let resonances = [(2, 1), (3, 2), (4, 3), (4, 2)];
        Some(resonances[rng.random_range(0..resonances.len())])
    } else {
        None
    }
}

/// Sample mass ratio (moon mass / planet mass) based on formation
fn sample_mass_ratio(rng: &mut ChaChaRng, category: usize, formation: &MoonFormation) -> f64 {
    match category {
        0 => {
            // Terrestrial: 0.1% - 5% (Earth-Moon is ~1.2%)
            match formation {
                MoonFormation::GiantImpact => rng.random_range(0.005..0.05),
                MoonFormation::Capture { .. } => rng.random_range(0.001..0.01),
                MoonFormation::CoAccretion { .. } => rng.random_range(0.002..0.02),
            }
        }
        1 => {
            // Super-Earth: 0.05% - 2%
            match formation {
                MoonFormation::GiantImpact => rng.random_range(0.002..0.02),
                MoonFormation::Capture { .. } => rng.random_range(0.0005..0.005),
                MoonFormation::CoAccretion { .. } => rng.random_range(0.001..0.01),
            }
        }
        2 => {
            // Ice Giant: 0.01% - 0.5%
            match formation {
                MoonFormation::GiantImpact => rng.random_range(0.001..0.005),
                MoonFormation::Capture { .. } => rng.random_range(0.0001..0.002),
                MoonFormation::CoAccretion { .. } => rng.random_range(0.0005..0.003),
            }
        }
        _ => {
            // Gas Giant: 0.001% - 0.1% (Ganymede is ~0.008% of Jupiter)
            match formation {
                MoonFormation::GiantImpact => rng.random_range(0.0001..0.001),
                MoonFormation::Capture { .. } => rng.random_range(0.00005..0.0005),
                MoonFormation::CoAccretion { .. } => rng.random_range(0.0001..0.001),
            }
        }
    }
}

/// Sample moon composition based on planet category and tidal heating
fn sample_moon_composition(rng: &mut ChaChaRng, category: usize, heat_flux: f64) -> Composition {
    // High tidal heating drives off volatiles
    if heat_flux > 1.0 {
        // Io-like: rocky/sulfur, no water
        return Composition::new(
            0.10 + rng.random::<f64>() * 0.10,
            0.80 + rng.random::<f64>() * 0.10,
            0.0,
            0.0,
        );
    }

    match category {
        0 | 1 => {
            // Terrestrial/Super-Earth moons: rocky with some iron
            Composition::new(
                0.20 + rng.random::<f64>() * 0.15,
                0.65 + rng.random::<f64>() * 0.15,
                rng.random::<f64>() * 0.10,
                0.0,
            )
        }
        _ => {
            // Ice/Gas giant moons: icy composition
            Composition::new(
                0.05 + rng.random::<f64>() * 0.10,
                0.25 + rng.random::<f64>() * 0.15,
                0.50 + rng.random::<f64>() * 0.20,
                0.0,
            )
        }
    }
}
