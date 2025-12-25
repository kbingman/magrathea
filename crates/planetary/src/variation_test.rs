//! Tests for stochastic variation

use crate::variation::PlanetaryVariation;

#[test]
fn test_deterministic_seeding() {
    // Same planet should produce same variation sequence
    let mut var1 = PlanetaryVariation::new_from_planet(1.0, 1.0);
    let mut var2 = PlanetaryVariation::new_from_planet(1.0, 1.0);

    let val1 = var1.vary_composition(0.5, 0.1);
    let val2 = var2.vary_composition(0.5, 0.1);

    assert_eq!(val1, val2, "Same planet should produce same variation");
}

#[test]
fn test_different_planets_different_variation() {
    let mut var1 = PlanetaryVariation::new_from_planet(1.0, 1.0);
    let mut var2 = PlanetaryVariation::new_from_planet(2.0, 1.0);

    let val1 = var1.vary_composition(0.5, 0.1);
    let val2 = var2.vary_composition(0.5, 0.1);

    assert_ne!(
        val1, val2,
        "Different planets should have different variation"
    );
}

#[test]
fn test_no_seed_collision_for_different_inputs() {
    // These pairs would have collided with the old linear seed formula:
    // old: seed = mass * 31000 + distance * 1000
    // (1.0, 31.0) -> 31000 + 31000 = 62000
    // (2.0, 0.0)  -> 62000 + 0     = 62000
    let mut var1 = PlanetaryVariation::new_from_planet(1.0, 31.0);
    let mut var2 = PlanetaryVariation::new_from_planet(2.0, 0.0);

    let val1 = var1.random_unit();
    let val2 = var2.random_unit();

    assert_ne!(val1, val2, "These should no longer collide");
}

#[test]
fn test_vary_composition_bounds() {
    let mut var = PlanetaryVariation::new_from_planet(1.0, 1.0);

    // Test many iterations to ensure bounds are respected
    for _ in 0..100 {
        let result = var.vary_composition(0.5, 0.1);
        assert!(result >= 0.0, "Composition should be >= 0.0");
        assert!(result <= 1.0, "Composition should be <= 1.0");
    }
}

#[test]
fn test_vary_composition_clamping() {
    let mut var = PlanetaryVariation::from_seed(42);

    // Edge case: high base value with variation shouldn't exceed 1.0
    for _ in 0..100 {
        let result = var.vary_composition(0.95, 0.2);
        assert!(result <= 1.0, "Should clamp to 1.0");
    }

    // Edge case: low base value with variation shouldn't go below 0.0
    for _ in 0..100 {
        let result = var.vary_composition(0.05, 0.2);
        assert!(result >= 0.0, "Should clamp to 0.0");
    }
}

#[test]
fn test_vary_property() {
    let mut var = PlanetaryVariation::from_seed(123);

    // 5% variation on base value of 100.0
    let result = var.vary_property(100.0, 5.0);
    assert!(result >= 95.0, "Should be within -5%");
    assert!(result <= 105.0, "Should be within +5%");
}

#[test]
fn test_vary_albedo() {
    let mut var = PlanetaryVariation::from_seed(456);

    // Test albedo variation respects bounds
    for _ in 0..100 {
        let result = var.vary_albedo(0.3, 0.1, 0.8);
        assert!(result >= 0.1, "Should respect min bound");
        assert!(result <= 0.8, "Should respect max bound");
    }
}

#[test]
fn test_vary_color() {
    let mut var = PlanetaryVariation::from_seed(789);

    let base_color = (0.5, 0.6, 0.7);
    let (r, g, b) = var.vary_color(base_color, 0.1);

    // Each component should be within [0, 1]
    assert!(r >= 0.0 && r <= 1.0, "Red should be in [0, 1]");
    assert!(g >= 0.0 && g <= 1.0, "Green should be in [0, 1]");
    assert!(b >= 0.0 && b <= 1.0, "Blue should be in [0, 1]");
}

#[test]
fn test_vary_color_clamping() {
    let mut var = PlanetaryVariation::from_seed(101);

    // Edge case: near-white color shouldn't exceed 1.0
    for _ in 0..100 {
        let (r, g, b) = var.vary_color((0.98, 0.98, 0.98), 0.1);
        assert!(r <= 1.0 && g <= 1.0 && b <= 1.0, "Should clamp to 1.0");
    }

    // Edge case: near-black color shouldn't go below 0.0
    for _ in 0..100 {
        let (r, g, b) = var.vary_color((0.02, 0.02, 0.02), 0.1);
        assert!(r >= 0.0 && g >= 0.0 && b >= 0.0, "Should clamp to 0.0");
    }
}

#[test]
fn test_orbital_variation() {
    let mut var = PlanetaryVariation::from_seed(202);

    for _ in 0..100 {
        let result = var.orbital_variation();
        assert!(result >= -0.05, "Should be >= -0.05");
        assert!(result <= 0.05, "Should be <= 0.05");
    }
}

#[test]
fn test_random_int() {
    let mut var = PlanetaryVariation::from_seed(303);

    for _ in 0..100 {
        let result = var.random_int(1, 10);
        assert!(result >= 1 && result <= 10, "Should be in range [1, 10]");
    }
}

#[test]
fn test_random_bool() {
    let mut var = PlanetaryVariation::from_seed(404);

    // With 0.0 probability, should always be false
    for _ in 0..10 {
        assert!(!var.random_bool(0.0), "0.0 probability should be false");
    }

    // With 1.0 probability, should always be true
    for _ in 0..10 {
        assert!(var.random_bool(1.0), "1.0 probability should be true");
    }

    // With 0.5 probability, we should see some of each (statistically)
    let mut var2 = PlanetaryVariation::from_seed(505);
    let trues: i32 = (0..100)
        .map(|_| if var2.random_bool(0.5) { 1 } else { 0 })
        .sum();
    assert!(
        trues > 20 && trues < 80,
        "0.5 probability should give ~50% true"
    );
}

#[test]
fn test_random_unit() {
    let mut var = PlanetaryVariation::from_seed(606);

    for _ in 0..100 {
        let result = var.random_unit();
        assert!(result >= 0.0, "Should be >= 0.0");
        assert!(result < 1.0, "Should be < 1.0");
    }
}

#[test]
fn test_random_range() {
    let mut var = PlanetaryVariation::from_seed(707);

    for _ in 0..100 {
        let result = var.random_range(5.0, 10.0);
        assert!(
            result >= 5.0 && result <= 10.0,
            "Should be in range [5, 10]"
        );
    }
}

#[test]
fn test_from_seed() {
    // Same seed should produce same sequence
    let mut var1 = PlanetaryVariation::from_seed(12345);
    let mut var2 = PlanetaryVariation::from_seed(12345);

    assert_eq!(
        var1.random_unit(),
        var2.random_unit(),
        "Same seed should produce same sequence"
    );
}

#[test]
fn test_sequence_reproducibility() {
    // Generate a sequence, then regenerate and compare
    let mut var1 = PlanetaryVariation::new_from_planet(5.5, 2.3);
    let seq1: Vec<f64> = (0..10).map(|_| var1.random_unit()).collect();

    let mut var2 = PlanetaryVariation::new_from_planet(5.5, 2.3);
    let seq2: Vec<f64> = (0..10).map(|_| var2.random_unit()).collect();

    assert_eq!(seq1, seq2, "Sequences should be identical for same planet");
}
