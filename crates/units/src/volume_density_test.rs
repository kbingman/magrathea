mod tests {
    use approx::assert_relative_eq;

    use crate::volume_density::VolumeDensity;

    #[test]
    fn test_volume_density_conversions() {
        // Test g/cm³ to kg/m³ conversion
        let density = VolumeDensity::from_grams_per_cm3(5.5);
        assert_relative_eq!(density.to_kg_per_m3(), 5500.0);

        // Round trip test
        let density_kg = VolumeDensity::from_kg_per_m3(5500.0);
        assert_relative_eq!(density_kg.to_grams_per_cm3(), 5.5);

        // Water density
        let water = VolumeDensity::from_kg_per_m3(1000.0);
        assert_relative_eq!(water.to_grams_per_cm3(), 1.0);
    }

    #[test]
    fn test_volume_density_arithmetic() {
        let density1 = VolumeDensity::from_grams_per_cm3(5.0);
        let density2 = VolumeDensity::from_grams_per_cm3(3.0);

        // Test addition and subtraction
        assert_relative_eq!((density1 + density2).to_grams_per_cm3(), 8.0);
        assert_relative_eq!((density1 - density2).to_grams_per_cm3(), 2.0);

        // Test multiplication and division
        let doubled = density1 * 2.0;
        assert_relative_eq!(doubled.to_grams_per_cm3(), 10.0);

        let halved = density1 / 2.0;
        assert_relative_eq!(halved.to_grams_per_cm3(), 2.5);
    }

    #[test]
    fn test_material_constants() {
        // Test that material densities are reasonable
        let iron = VolumeDensity::iron();
        assert_relative_eq!(iron.to_grams_per_cm3(), 7.9);

        let rock = VolumeDensity::silicate_rock();
        assert_relative_eq!(rock.to_grams_per_cm3(), 3.3);

        let ice = VolumeDensity::water_ice();
        assert_relative_eq!(ice.to_grams_per_cm3(), 0.92, epsilon = 0.01);

        let methane = VolumeDensity::methane_ice();
        assert_relative_eq!(methane.to_grams_per_cm3(), 0.47, epsilon = 0.01);
    }

    #[test]
    fn test_weighted_average() {
        // Earth-like: 32% iron, 68% rock (uncompressed)
        let components = vec![
            (VolumeDensity::iron(), 0.32),
            (VolumeDensity::silicate_rock(), 0.68),
        ];

        let bulk = VolumeDensity::weighted_average(&components).unwrap();
        let expected = 7.9 * 0.32 + 3.3 * 0.68;
        assert_relative_eq!(bulk.to_grams_per_cm3(), expected, epsilon = 0.01);

        // Ice-rock mixture (50/50)
        let ice_rock = vec![
            (VolumeDensity::water_ice(), 0.5),
            (VolumeDensity::silicate_rock(), 0.5),
        ];

        let ice_planet = VolumeDensity::weighted_average(&ice_rock).unwrap();
        let expected_ice = 0.92 * 0.5 + 3.3 * 0.5;
        assert_relative_eq!(ice_planet.to_grams_per_cm3(), expected_ice, epsilon = 0.01);
    }

    #[test]
    fn test_weighted_average_invalid() {
        // Fractions don't sum to 1
        let bad_components = vec![
            (VolumeDensity::iron(), 0.5),
            (VolumeDensity::silicate_rock(), 0.3),
            // Only sums to 0.8
        ];

        assert!(VolumeDensity::weighted_average(&bad_components).is_none());
    }

    #[test]
    fn test_realistic_planetary_densities() {
        // Earth: ~5.5 g/cm³
        let earth = VolumeDensity::from_grams_per_cm3(5.5);
        assert!(earth.to_grams_per_cm3() > 5.0);
        assert!(earth.to_grams_per_cm3() < 6.0);

        // Mars: ~3.9 g/cm³ (less iron)
        let mars = VolumeDensity::from_grams_per_cm3(3.9);
        assert!(mars.to_grams_per_cm3() < earth.to_grams_per_cm3());

        // Jupiter: ~1.3 g/cm³ (gas giant)
        let jupiter = VolumeDensity::from_grams_per_cm3(1.3);
        assert!(jupiter.to_grams_per_cm3() < 2.0);

        // Mercury: ~5.4 g/cm³ (high iron content)
        let mercury = VolumeDensity::from_grams_per_cm3(5.4);
        assert!(mercury.to_grams_per_cm3() > 5.0);
    }

    #[test]
    fn test_unit_conversion_consistency() {
        // 1 g/cm³ should equal 1000 kg/m³
        let density = VolumeDensity::from_grams_per_cm3(1.0);
        assert_relative_eq!(density.to_kg_per_m3(), 1000.0);

        // Verify the conversion factor
        let density_kg = VolumeDensity::from_kg_per_m3(1000.0);
        assert_relative_eq!(density_kg.to_grams_per_cm3(), 1.0);
    }
}
