mod tests {
    use approx::assert_relative_eq;

    use crate::mass_rate::MassRate;
    use crate::time::Time;

    #[test]
    fn test_mass_rate_conversions() {
        // Test Earth masses per Myr to solar masses per year
        let rate_em = MassRate::from_earth_masses_per_myr(60.0);
        let solar_rate = rate_em.to_solar_masses_per_year();

        // Round trip test
        let round_trip = MassRate::from_solar_masses_per_year(solar_rate);
        assert_relative_eq!(round_trip.to_earth_masses_per_myr(), 60.0, epsilon = 1e-6);

        // Test grams per year conversions
        let rate_g = MassRate::from_grams_per_year(1e20);
        let solar_g = rate_g.to_solar_masses_per_year();
        let round_trip_g = MassRate::from_solar_masses_per_year(solar_g);
        assert_relative_eq!(round_trip_g.to_grams_per_year(), 1e20, epsilon = 1e10);
    }

    #[test]
    fn test_mass_rate_integration() {
        // Test that integrating mass rate over time gives correct mass
        let flux = MassRate::from_earth_masses_per_myr(60.0);
        let duration = Time::from_years(1_000_000.0); // 1 Myr

        let total_mass = flux.integrate(duration);

        // Should be approximately 60 Earth masses
        assert_relative_eq!(total_mass.to_earth_masses(), 60.0, epsilon = 1e-3);

        // Test with different time scale
        let flux_solar = MassRate::from_solar_masses_per_year(1e-10);
        let duration_short = Time::from_years(1_000.0);
        let mass_solar = flux_solar.integrate(duration_short);

        assert_relative_eq!(
            mass_solar.to_solar_masses(),
            1e-10 * 1000.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn test_mass_rate_arithmetic() {
        let rate1 = MassRate::from_earth_masses_per_myr(60.0);
        let rate2 = MassRate::from_earth_masses_per_myr(40.0);

        // Test addition and subtraction
        assert_relative_eq!(
            (rate1 + rate2).to_earth_masses_per_myr(),
            100.0,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            (rate1 - rate2).to_earth_masses_per_myr(),
            20.0,
            epsilon = 1e-6
        );

        // Test multiplication and division
        let doubled = rate1 * 2.0;
        assert_relative_eq!(doubled.to_earth_masses_per_myr(), 120.0, epsilon = 1e-6);

        let halved = rate1 / 2.0;
        assert_relative_eq!(halved.to_earth_masses_per_myr(), 30.0, epsilon = 1e-6);
    }

    #[test]
    fn test_realistic_planetary_formation_rates() {
        // Typical pebble accretion flux
        let pebble_flux = MassRate::from_earth_masses_per_myr(60.0);

        // Accrete for 100,000 years
        let time = Time::from_years(100_000.0);
        let accreted = pebble_flux.integrate(time);

        // Should accrete 6 Earth masses
        assert_relative_eq!(accreted.to_earth_masses(), 6.0, epsilon = 1e-3);

        // Stellar wind (very small)
        let stellar_wind = MassRate::from_solar_masses_per_year(1e-14);

        // Over 1 billion years
        let long_time = Time::from_years(1e9);
        let mass_lost = stellar_wind.integrate(long_time);

        // Should lose about 1e-5 solar masses
        assert_relative_eq!(mass_lost.to_solar_masses(), 1e-5, epsilon = 1e-7);
    }
}
