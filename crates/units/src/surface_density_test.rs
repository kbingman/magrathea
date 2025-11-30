mod tests {
    use approx::assert_relative_eq;

    use crate::surface_density::SurfaceDensity;

    #[test]
    fn test_surface_density_conversions() {
        // Test g/cm² to kg/m² conversion
        let sigma = SurfaceDensity::from_grams_per_cm2(100.0);
        assert_relative_eq!(sigma.to_kg_per_m2(), 1000.0);

        // Round trip test
        let sigma_kg = SurfaceDensity::from_kg_per_m2(1000.0);
        assert_relative_eq!(sigma_kg.to_grams_per_cm2(), 100.0);

        // Test solar mass per AU² conversions
        let sigma_solar = SurfaceDensity::from_solar_masses_per_au2(1e-6);
        let g_cm2 = sigma_solar.to_grams_per_cm2();
        let round_trip = SurfaceDensity::from_grams_per_cm2(g_cm2);
        assert_relative_eq!(round_trip.to_solar_masses_per_au2(), 1e-6, epsilon = 1e-10);
    }

    #[test]
    fn test_surface_density_arithmetic() {
        let sigma1 = SurfaceDensity::from_grams_per_cm2(1000.0);
        let sigma2 = SurfaceDensity::from_grams_per_cm2(500.0);

        // Test addition and subtraction
        assert_relative_eq!((sigma1 + sigma2).to_grams_per_cm2(), 1500.0);
        assert_relative_eq!((sigma1 - sigma2).to_grams_per_cm2(), 500.0);

        // Test multiplication and division
        let doubled = sigma1 * 2.0;
        assert_relative_eq!(doubled.to_grams_per_cm2(), 2000.0);

        let halved = sigma1 / 2.0;
        assert_relative_eq!(halved.to_grams_per_cm2(), 500.0);
    }

    #[test]
    fn test_power_law_scaling() {
        // Minimum Mass Solar Nebula at 1 AU
        let sigma_1au = SurfaceDensity::from_grams_per_cm2(1700.0);

        // At 5 AU with p=1.5 power law: Σ(5) = Σ₀ * 5^(-1.5)
        let sigma_5au = sigma_1au.power_law_scaling(5.0, 1.5);
        let expected = 1700.0 / 5.0_f64.powf(1.5);
        assert_relative_eq!(sigma_5au.to_grams_per_cm2(), expected, epsilon = 1e-3);

        // At 0.1 AU (should be much higher)
        let sigma_01au = sigma_1au.power_law_scaling(0.1, 1.5);
        let expected_01 = 1700.0 / 0.1_f64.powf(1.5);
        assert_relative_eq!(sigma_01au.to_grams_per_cm2(), expected_01, epsilon = 1.0);
    }

    #[test]
    fn test_realistic_disk_values() {
        // Minimum Mass Solar Nebula
        let mmsn = SurfaceDensity::from_grams_per_cm2(1700.0);
        assert!(mmsn.to_grams_per_cm2() > 0.0);

        // At 5 AU (Jupiter region), should be ~150-200 g/cm²
        let jupiter_region = mmsn.power_law_scaling(5.0, 1.5);
        assert!(jupiter_region.to_grams_per_cm2() > 100.0);
        assert!(jupiter_region.to_grams_per_cm2() < 300.0);

        // At 30 AU (Neptune region), should be much lower
        let neptune_region = mmsn.power_law_scaling(30.0, 1.5);
        assert!(neptune_region.to_grams_per_cm2() < 50.0);
        assert!(neptune_region.to_grams_per_cm2() > 0.0);

        // Depleted disk (100x lower)
        let depleted = mmsn / 100.0;
        assert_relative_eq!(depleted.to_grams_per_cm2(), 17.0);
    }

    #[test]
    fn test_unit_conversion_consistency() {
        // 1 g/cm² should equal 10 kg/m²
        let sigma = SurfaceDensity::from_grams_per_cm2(1.0);
        assert_relative_eq!(sigma.to_kg_per_m2(), 10.0);

        // Verify the conversion factor
        let sigma_kg = SurfaceDensity::from_kg_per_m2(10.0);
        assert_relative_eq!(sigma_kg.to_grams_per_cm2(), 1.0);
    }
}
