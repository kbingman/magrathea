// use crate::{
//     planet::earth_analog,
//     system::{
//         HabitableZone, PlanetarySystem, SystemArchitecture, calculate_mutual_hill_radius, snow_line,
//     },
// };
// use rand::SeedableRng;
// use rand_chacha::ChaChaRng;
// use units::Length;

// #[test]
// fn test_architecture_sampling_by_spectral_type() {
//     let mut rng = ChaChaRng::seed_from_u64(42);

//     // M dwarfs should favor compact multi
//     let mut compact_count = 0;
//     for _ in 0..100 {
//         if matches!(
//             SystemArchitecture::sample(&mut rng, "M", 0.0),
//             SystemArchitecture::CompactMulti
//         ) {
//             compact_count += 1;
//         }
//     }
//     assert!(
//         compact_count > 30,
//         "M dwarfs should frequently produce compact systems"
//     );

//     // High metallicity should increase giant probability
//     let mut _giant_low = 0;
//     let mut _giant_high = 0;
//     for _ in 0..100 {
//         if matches!(
//             SystemArchitecture::sample(&mut rng, "G", -0.3),
//             SystemArchitecture::GiantDominated
//         ) {
//             _giant_low += 1;
//         }
//         if matches!(
//             SystemArchitecture::sample(&mut rng, "G", 0.3),
//             SystemArchitecture::GiantDominated
//         ) {
//             _giant_high += 1;
//         }
//     }
//     // Can't guarantee with RNG, but high metallicity should tend toward more giants
// }

// #[test]
// fn test_habitable_zone() {
//     // Sun-like star
//     let hz = HabitableZone::from_luminosity(1.0);
//     assert!(hz.inner_edge < 1.0 && hz.inner_edge > 0.8);
//     assert!(hz.outer_edge > 1.5 && hz.outer_edge < 2.0);
//     assert!(hz.contains(1.0)); // Earth is in HZ

//     // M dwarf (L = 0.01 L☉)
//     let hz_m = HabitableZone::from_luminosity(0.01);
//     assert!(hz_m.inner_edge < 0.15);
//     assert!(hz_m.outer_edge < 0.3);
// }

// #[test]
// fn test_snow_line() {
//     // Sun-like: snow line ~2.7 AU
//     let sl = snow_line(1.0);
//     assert!(sl > 2.0 && sl < 4.0);

//     // M dwarf: much closer
//     let sl_m = snow_line(0.01);
//     assert!(sl_m < 0.5);
// }

// #[test]
// fn test_mutual_hill_radius() {
//     // Earth-Mars at 1.0 and 1.5 AU
//     let r_h = calculate_mutual_hill_radius(1.0, 0.107, 1.0, 1.5, 1.0);
//     // Should be small fraction of separation
//     assert!(r_h < 0.1);
// }

// #[test]
// fn test_system_stability() {
//     let earth = earth_analog();
//     let mut earth2 = earth.clone();
//     earth2.semi_major_axis = Length::from_au(1.01); // Way too close

//     let system = PlanetarySystem::new(
//         1.0,
//         1.0,
//         5778.0,
//         0.0,
//         "G".to_string(),
//         vec![earth, earth2],
//         SystemArchitecture::Mixed,
//     );

//     assert!(
//         !system.is_stable(),
//         "Planets at 1.0 and 1.01 AU should be unstable"
//     );
// }
