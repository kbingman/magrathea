use approx::assert_relative_eq;

use crate::composition::Composition;

#[test]
fn test_fractions_sum_to_one() {
    let compositions = [
        Composition::earth_like(),
        Composition::iron_rich(),
        Composition::water_world(),
        Composition::mini_neptune(),
        Composition::ice_giant(),
        Composition::gas_giant(),
    ];

    for comp in compositions {
        let sum = comp.iron + comp.silicate + comp.water + comp.h_he_gas;

        assert_relative_eq!(sum.abs(), 1.0, epsilon = 1e-3);
    }
}

#[test]
fn test_stripped_envelope() {
    let mini_nep = Composition::mini_neptune();
    assert!(mini_nep.h_he_gas > 0.0);

    let stripped = mini_nep.stripped_envelope();
    assert_eq!(stripped.h_he_gas, 0.0);

    let sum = stripped.iron + stripped.silicate + stripped.water;
    assert_relative_eq!(sum.abs(), 1.0, epsilon = 1e-10);
}
