use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use units::{Length, Mass, Temperature, Time};

use crate::stellar::{
    spectral_class::{LuminosityClass, SpectralType, VariabilityType},
    stellar_object::{
        BlackHole, GiantStar, MainSequenceStar, NeutronStar, WhiteDwarf, WhiteDwarfType,
    },
};

#[test]
fn white_dwarf_test() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let star = WhiteDwarf::new(&mut rng);

    assert_eq!(
        star,
        WhiteDwarf {
            mass: Mass::from_solar_masses(0.8228426056215805),
            radius: Length::from_solar_radii(0.009000738823795101),
            luminosity: 0.0013714043427026344,
            temperature: Temperature::from_kelvin(11720.060071150392),
            spectral_type: WhiteDwarfType::DB,
        }
    );
}

#[test]
fn neutron_star_test() {
    let mut rng = ChaChaRng::seed_from_u64(43);
    let star = NeutronStar::new(&mut rng);

    assert_eq!(
        star,
        NeutronStar {
            mass: Mass::from_solar_masses(1.4249211795524308),
            radius: Length::from_km(11.936648469868707),
            magnetic_field: 11.362101818674972,
            magnetar: false,
            pulsar: true,
        }
    );
}

#[test]
fn black_hole_test() {
    let mut rng = ChaChaRng::seed_from_u64(43);
    let star = BlackHole::new(&mut rng, 21.0, 1.0);

    assert_eq!(
        star,
        BlackHole {
            mass: Mass::from_solar_masses(4.418330033588915),
            has_accretion: false,
            spin: 0.15309853235773063,
        }
    );
}

#[test]
fn main_sequence_test() {
    let star = MainSequenceStar::solar_analog();

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
        }
    );
}

#[test]
fn giant_star_test() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let star = GiantStar::new(&mut rng, 34.0);

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
        }
    );
}
