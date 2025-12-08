import { describe, expect, test } from "vitest";
import {
  scaleRadius,
  scaleSemiMajorAxis,
  scalePlanetaryMass,
  scaleStellarRadius,
  SEMI_MAJOR_AXIS_SCALE,
  PLANETARY_MASS_SCALE,
  STELLAR_RADIUS_SCALE,
  scaleUnifiedMass,
  UNIFIED_MASS_SCALE,
} from "../src/scale";

describe("scaleRadius", () => {
  test("clamps to a min radius", () => {
    expect(scaleRadius(0.0000003)).toBe(1);
  });

  test("clamps to a custom min radius", () => {
    expect(scaleRadius(0.0000003, { minOutput: 0.5 })).toBe(0.5);
  });

  test("clamps to a max radius", () => {
    expect(scaleRadius(256)).toBeLessThan(12);
  });

  test("clamps to a custom max radius", () => {
    const radius = scaleRadius(256, { maxOutput: 16 });
    expect(radius).toBeGreaterThan(8);
    expect(radius).toBeLessThan(24);
  });
});

describe("scaleSemiMajorAxis", () => {
  test("scales Mercury-like orbit (0.39 AU)", () => {
    const scaled = scaleSemiMajorAxis(0.39);
    expect(scaled).toBeGreaterThan(10);
    expect(scaled).toBeLessThan(100);
  });

  test("scales Earth-like orbit (1.0 AU)", () => {
    const scaled = scaleSemiMajorAxis(1.0);
    expect(scaled).toBeGreaterThan(50);
    expect(scaled).toBeLessThan(200);
  });

  test("scales distant orbit (50 AU)", () => {
    const scaled = scaleSemiMajorAxis(50);
    expect(scaled).toBeGreaterThan(300);
    expect(scaled).toBeLessThan(400);
  });

  test("allows custom options", () => {
    const scaled = scaleSemiMajorAxis(1.0, {
      minOutput: 100,
      maxOutput: 500,
    });
    expect(scaled).toBeGreaterThan(100);
    expect(scaled).toBeLessThan(500);
  });
});

describe("scalePlanetaryMass", () => {
  test("scales small moon mass (0.001 Earth masses)", () => {
    const scaled = scalePlanetaryMass(0.001);
    expect(scaled).toBeGreaterThan(2);
    expect(scaled).toBeLessThan(8);
  });

  test("scales Earth-like mass (1.0 Earth masses)", () => {
    const scaled = scalePlanetaryMass(1.0);
    expect(scaled).toBeGreaterThan(8);
    expect(scaled).toBeLessThan(14);
  });

  test("scales Jupiter-like mass (10.0 Earth masses)", () => {
    const scaled = scalePlanetaryMass(10.0);
    expect(scaled).toBeGreaterThan(12);
    expect(scaled).toBeLessThan(15);
  });

  test("allows custom options", () => {
    const scaled = scalePlanetaryMass(1.0, { minOutput: 5, maxOutput: 25 });
    expect(scaled).toBeGreaterThan(5);
    expect(scaled).toBeLessThan(25);
  });
});

describe("scaleStellarRadius", () => {
  test("uses existing stellar radius defaults", () => {
    const scaled = scaleStellarRadius(1.0);
    expect(scaled).toBeGreaterThan(1);
    expect(scaled).toBeLessThan(8);
  });
});

describe("Scale configurations", () => {
  test("SEMI_MAJOR_AXIS_SCALE has correct defaults", () => {
    expect(SEMI_MAJOR_AXIS_SCALE.minInput).toBe(0.1);
    expect(SEMI_MAJOR_AXIS_SCALE.maxInput).toBe(100.0);
    expect(SEMI_MAJOR_AXIS_SCALE.minOutput).toBe(10);
    expect(SEMI_MAJOR_AXIS_SCALE.maxOutput).toBe(400);
  });

  test("PLANETARY_MASS_SCALE has correct defaults", () => {
    expect(PLANETARY_MASS_SCALE.minInput).toBe(0.000001);
    expect(PLANETARY_MASS_SCALE.maxInput).toBe(20.0);
    expect(PLANETARY_MASS_SCALE.minOutput).toBe(2);
    expect(PLANETARY_MASS_SCALE.maxOutput).toBe(15);
  });

  test("STELLAR_RADIUS_SCALE has correct defaults", () => {
    expect(STELLAR_RADIUS_SCALE.minInput).toBe(0.001);
    expect(STELLAR_RADIUS_SCALE.maxInput).toBe(650);
    expect(STELLAR_RADIUS_SCALE.minOutput).toBe(1);
    expect(STELLAR_RADIUS_SCALE.maxOutput).toBe(8);
  });
});

describe("scaleUnifiedMass", () => {
  test("scales Pluto-sized mass", () => {
    const scaled = scaleUnifiedMass(0.0000000066); // Pluto-sized object
    expect(scaled).toBeGreaterThanOrEqual(1);
    expect(scaled).toBeLessThan(5);
  });

  test("scales Earth-like mass", () => {
    const scaled = scaleUnifiedMass(0.000003); // ~1 Earth mass
    expect(scaled).toBeGreaterThan(1);
    expect(scaled).toBeLessThan(16);
  });

  test("scales Jupiter-like mass", () => {
    const scaled = scaleUnifiedMass(0.001); // ~Jupiter mass
    expect(scaled).toBeGreaterThan(20);
    expect(scaled).toBeLessThan(25);
  });

  test("scales stellar mass", () => {
    const scaled = scaleUnifiedMass(1.0); // Sun mass
    expect(scaled).toBeGreaterThan(30);
    expect(scaled).toBeLessThan(40);
  });

  test("scales large stellar mass", () => {
    const scaled = scaleUnifiedMass(5.0); // Large star
    expect(scaled).toBeGreaterThan(38);
    expect(scaled).toBeLessThan(40);
  });

  test("allows custom options", () => {
    const scaled = scaleUnifiedMass(1.0, { minOutput: 25, maxOutput: 50 });
    expect(scaled).toBeGreaterThan(25);
    expect(scaled).toBeLessThan(50);
  });
});

describe("UNIFIED_MASS_SCALE configuration", () => {
  test("has correct defaults for unified scaling", () => {
    expect(UNIFIED_MASS_SCALE.minInput).toBe(0.0000000066);
    expect(UNIFIED_MASS_SCALE.maxInput).toBe(10.0);
    expect(UNIFIED_MASS_SCALE.minOutput).toBe(1);
    expect(UNIFIED_MASS_SCALE.maxOutput).toBe(40);
  });
});
