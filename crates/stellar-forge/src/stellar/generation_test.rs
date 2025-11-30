use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use units::{Length, Mass, Temperature, Time};

use crate::stellar::generation::{
    black_hole, calculate_radius, calculate_subtype, estimate_lifetime, giant_star,
    main_sequence_star, neutron_star, sample_main_sequence_star, sample_stellar_object,
    solar_analog, spectral_type_from_temp, stellar_object, white_dwarf,
};
use crate::stellar::spectral::{LuminosityClass, SpectralType, VariabilityType};
use crate::stellar::stellar_color::StellarColor;
use crate::stellar::stellar_objects::{
    BlackHole, GiantStar, MainSequenceStar, NeutronStar, StellarObject, WhiteDwarf, WhiteDwarfType,
};

#[test]
fn white_dwarf_test() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let star = white_dwarf(&mut rng);

    assert_eq!(
        star,
        WhiteDwarf {
            mass: Mass::from_solar_masses(0.8228426056215805),
            radius: Length::from_solar_radii(0.009000738823795101),
            luminosity: 0.0013714043427026344,
            temperature: Temperature::from_kelvin(11720.060071150392),
            spectral_type: WhiteDwarfType::DB,
            color: StellarColor::from_temperature(11720.060071150392),
        }
    );
}

#[test]
fn neutron_star_test() {
    let mut rng = ChaChaRng::seed_from_u64(43);
    let star = neutron_star(&mut rng);

    assert_eq!(
        star,
        NeutronStar {
            mass: Mass::from_solar_masses(1.4249211795524308),
            radius: Length::from_km(11.936648469868707),
            magnetic_field: 11.362101818674972,
            magnetar: false,
            pulsar: true,
            color: StellarColor::new(180, 200, 255), // Blue-white (pulsar)
        }
    );
}

#[test]
fn black_hole_test() {
    let mut rng = ChaChaRng::seed_from_u64(43);
    let star = black_hole(&mut rng, 21.0, 1.0);

    assert_eq!(
        star,
        BlackHole {
            mass: Mass::from_solar_masses(4.418330033588915),
            has_accretion: false,
            spin: 0.15309853235773063,
            color: StellarColor::new(20, 20, 30), // Nearly black
        }
    );
}

#[test]
fn main_sequence_test() {
    let star = solar_analog();

    assert_eq!(
        star,
        MainSequenceStar {
            mass: Mass::from_solar_masses(1.0),
            radius: Length::from_solar_radii(0.9924281807372175),
            luminosity: 1.0,
            metallicity: 0.0,
            age: Time::from_myr(1.0),
            temperature: Temperature::from_kelvin(5800.0),
            spectral_type: SpectralType::G,
            subtype: 2,
            luminosity_class: LuminosityClass::V,
            variability: VariabilityType::None,
            color: StellarColor::from_temperature(5800.0),
        }
    );
}

#[test]
fn giant_star_test() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let star = giant_star(&mut rng, 34.0);

    assert_eq!(
        star,
        GiantStar {
            mass: Mass::from_solar_masses(27.21910704224043),
            radius: Length::from_solar_radii(124.40092153463617),
            luminosity: 642222.2222222221,
            temperature: Temperature::from_kelvin(14665.16896210703),
            spectral_type: SpectralType::B,
            subtype: 7,
            luminosity_class: LuminosityClass::IAPLUS,
            variability: VariabilityType::None,
            color: StellarColor::from_temperature(14665.16896210703),
        }
    );
}

#[test]
fn main_sequence_star_from_mass() {
    let star = main_sequence_star(1.5, 0.1, 500.0);

    assert_eq!(star.mass.to_solar_masses(), 1.5);
    assert_eq!(star.metallicity, 0.1);
    assert_eq!(star.age.to_myr(), 500.0);
    assert_eq!(star.luminosity_class, LuminosityClass::V);
}

// ============================================================================
// sample_main_sequence_star tests
// ============================================================================

#[test]
fn sample_main_sequence_star_produces_valid_star() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let star = sample_main_sequence_star(&mut rng);

    // Should produce a valid main sequence star
    assert!(star.mass.to_solar_masses() >= 0.08);
    assert!(star.mass.to_solar_masses() <= 3.0); // Limited to planet-hosting range
    assert!(star.luminosity > 0.0);
    assert!(star.temperature.to_kelvin() > 0.0);
    assert_eq!(star.luminosity_class, LuminosityClass::V);
}

// ============================================================================
// stellar_object tests
// ============================================================================

#[test]
fn stellar_object_young_star_returns_main_sequence() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Young 1 solar mass star (age << lifetime)
    let obj = stellar_object(&mut rng, 1.0, 1.0e9, 0.0);
    assert!(matches!(obj, StellarObject::MainSequence(_)));
}

#[test]
fn stellar_object_old_low_mass_returns_white_dwarf() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Very old low mass star should become white dwarf
    // 1 solar mass lifetime ~10 Gyr, so 15 Gyr should make it a remnant
    let obj = stellar_object(&mut rng, 1.0, 15.0e9, 0.0);
    assert!(matches!(obj, StellarObject::WhiteDwarf(_)));
}

#[test]
fn stellar_object_old_intermediate_mass_returns_neutron_star() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // 10 solar mass star past its lifetime (~10 Myr) should become neutron star
    let obj = stellar_object(&mut rng, 10.0, 1.0e9, 0.0);
    assert!(matches!(obj, StellarObject::NeutronStar(_)));
}

#[test]
fn stellar_object_old_massive_returns_black_hole() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // 25 solar mass star past its lifetime should become black hole
    let obj = stellar_object(&mut rng, 25.0, 1.0e9, 0.0);
    assert!(matches!(obj, StellarObject::BlackHole(_)));
}

#[test]
fn stellar_object_giant_phase() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // 3 solar mass star at ~92% of lifetime should be a giant
    // Lifetime for 3 M☉ is roughly 0.37 Gyr
    let lifetime_gyr = estimate_lifetime(3.0);
    let age_years = lifetime_gyr * 0.92 * 1.0e9;
    let obj = stellar_object(&mut rng, 3.0, age_years, 0.0);
    assert!(matches!(obj, StellarObject::Giant(_)));
}

// ============================================================================
// sample_stellar_object tests
// ============================================================================

#[test]
fn sample_stellar_object_produces_valid_object() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Sample multiple objects
    for _ in 0..10 {
        let obj = sample_stellar_object(&mut rng);
        // Should produce some kind of stellar object
        match obj {
            StellarObject::MainSequence(star) => {
                assert!(star.mass.to_solar_masses() > 0.0);
            }
            StellarObject::Giant(star) => {
                assert!(star.mass.to_solar_masses() > 0.0);
            }
            StellarObject::WhiteDwarf(wd) => {
                assert!(wd.mass.to_solar_masses() > 0.0);
            }
            StellarObject::NeutronStar(ns) => {
                assert!(ns.mass.to_solar_masses() > 0.0);
            }
            StellarObject::BlackHole(bh) => {
                assert!(bh.mass.to_solar_masses() > 0.0);
            }
        }
    }
}

// ============================================================================
// Giant star mass range tests
// ============================================================================

#[test]
fn giant_star_supergiant_range() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // 15 M☉ should produce a supergiant (Ia or Ib)
    let star = giant_star(&mut rng, 15.0);
    assert!(
        star.luminosity_class == LuminosityClass::IA
            || star.luminosity_class == LuminosityClass::IB
    );
}

#[test]
fn giant_star_regular_giant_range() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // 5 M☉ should produce a regular giant (II or III)
    let star = giant_star(&mut rng, 5.0);
    assert!(
        star.luminosity_class == LuminosityClass::II
            || star.luminosity_class == LuminosityClass::III
    );
}

#[test]
fn giant_star_small_giant_range() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // 1.5 M☉ should produce a small giant (III or IV)
    let star = giant_star(&mut rng, 1.5);
    assert!(
        star.luminosity_class == LuminosityClass::III
            || star.luminosity_class == LuminosityClass::IV
    );
}

// ============================================================================
// Main sequence mass range tests
// ============================================================================

#[test]
fn main_sequence_star_very_massive() {
    // > 30 M☉
    let star = main_sequence_star(40.0, 0.0, 1.0);
    assert!(star.temperature.to_kelvin() > 30000.0);
    assert!(star.luminosity > 30000.0);
}

#[test]
fn main_sequence_star_massive() {
    // 8-30 M☉
    let star = main_sequence_star(15.0, 0.0, 1.0);
    assert!(star.temperature.to_kelvin() > 15000.0);
    assert!(star.luminosity > 1000.0);
}

#[test]
fn main_sequence_star_intermediate() {
    // 2-8 M☉
    let star = main_sequence_star(4.0, 0.0, 1.0);
    assert!(star.temperature.to_kelvin() > 7000.0);
    assert!(star.luminosity > 25.0);
}

#[test]
fn main_sequence_star_solar_type() {
    // 0.8-2 M☉
    let star = main_sequence_star(1.2, 0.0, 1.0);
    assert!(star.temperature.to_kelvin() > 5000.0);
    assert!(star.temperature.to_kelvin() < 7000.0);
}

#[test]
fn main_sequence_star_low_mass() {
    // < 0.8 M☉
    let star = main_sequence_star(0.5, 0.0, 1.0);
    assert!(star.temperature.to_kelvin() < 5000.0);
}

#[test]
fn main_sequence_star_very_low_mass() {
    // < 0.45 M☉ - different temperature formula
    let star = main_sequence_star(0.2, 0.0, 1.0);
    assert!(star.temperature.to_kelvin() < 4000.0);
}

// ============================================================================
// Helper function tests
// ============================================================================

#[test]
fn calculate_radius_solar_values() {
    // Solar luminosity (1.0) and solar temp (5778K) should give ~1 solar radius
    let radius = calculate_radius(1.0, 5778.0);
    assert!((radius - 1.0).abs() < 0.01);
}

#[test]
fn calculate_radius_hot_luminous() {
    // Higher luminosity + same temp = larger radius
    let radius = calculate_radius(100.0, 5778.0);
    assert!(radius > 1.0);
}

#[test]
fn spectral_type_from_temp_o_type() {
    assert_eq!(spectral_type_from_temp(35000.0), SpectralType::O);
}

#[test]
fn spectral_type_from_temp_b_type() {
    assert_eq!(spectral_type_from_temp(15000.0), SpectralType::B);
}

#[test]
fn spectral_type_from_temp_a_type() {
    assert_eq!(spectral_type_from_temp(8500.0), SpectralType::A);
}

#[test]
fn spectral_type_from_temp_f_type() {
    assert_eq!(spectral_type_from_temp(6500.0), SpectralType::F);
}

#[test]
fn spectral_type_from_temp_g_type() {
    assert_eq!(spectral_type_from_temp(5500.0), SpectralType::G);
}

#[test]
fn spectral_type_from_temp_k_type() {
    assert_eq!(spectral_type_from_temp(4500.0), SpectralType::K);
}

#[test]
fn spectral_type_from_temp_m_type() {
    assert_eq!(spectral_type_from_temp(3000.0), SpectralType::M);
}

#[test]
fn spectral_type_from_temp_l_type() {
    assert_eq!(spectral_type_from_temp(1800.0), SpectralType::L);
}

#[test]
fn spectral_type_from_temp_t_type() {
    assert_eq!(spectral_type_from_temp(1000.0), SpectralType::T);
}

#[test]
fn spectral_type_from_temp_y_type() {
    assert_eq!(spectral_type_from_temp(400.0), SpectralType::Y);
}

#[test]
fn calculate_subtype_range() {
    // Subtype should always be 0-9
    for temp in [35000.0, 15000.0, 8000.0, 6000.0, 5000.0, 4000.0, 3000.0] {
        let subtype = calculate_subtype(temp);
        assert!(subtype <= 9, "Subtype {} should be <= 9", subtype);
    }
}

#[test]
fn estimate_lifetime_massive_stars_short() {
    // Massive stars have short lifetimes
    let lifetime = estimate_lifetime(20.0);
    assert!(lifetime < 0.1, "20 M☉ star lifetime should be < 0.1 Gyr");
}

#[test]
fn estimate_lifetime_solar_type() {
    // Solar type stars have ~10 Gyr lifetimes
    let lifetime = estimate_lifetime(1.0);
    assert!(
        lifetime > 5.0 && lifetime < 15.0,
        "1 M☉ star lifetime should be ~10 Gyr"
    );
}

#[test]
fn estimate_lifetime_low_mass_long() {
    // Low mass stars have very long lifetimes
    let lifetime = estimate_lifetime(0.3);
    assert!(lifetime > 50.0, "0.3 M☉ star lifetime should be > 50 Gyr");
}
