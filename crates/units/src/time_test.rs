mod tests {
    use approx::assert_relative_eq;

    use crate::time::{Time, DAYS_PER_YEAR, HOURS_PER_YEAR, SECONDS_PER_YEAR};

    #[test]
    fn test_time_conversions() {
        // Test years to seconds
        let time_years = Time::from_years(1.0);
        assert_relative_eq!(time_years.to_seconds(), SECONDS_PER_YEAR);

        // Test seconds to years
        let time_seconds = Time::from_seconds(SECONDS_PER_YEAR);
        assert_relative_eq!(time_seconds.to_years(), 1.0);

        // Test days
        let days = 30.0;
        let time_days = Time::from_days(days);
        assert_relative_eq!(time_days.to_years(), days / DAYS_PER_YEAR);
        assert_relative_eq!(time_days.to_days(), days);

        // Test hours
        let hours = 240.0;
        let time_hours = Time::from_hours(hours);
        assert_relative_eq!(time_hours.to_years(), hours / HOURS_PER_YEAR);
        assert_relative_eq!(time_hours.to_hours(), hours);

        // Test addition
        let a = Time::from_days(10.0);
        let b = Time::from_days(5.0);
        let sum = a + b;
        assert_relative_eq!(sum.to_days(), 15.0);
        assert_relative_eq!(sum.to_years(), 0.04106776, epsilon = 1e-6);
    }
}
