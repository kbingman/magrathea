type ScaleOptions = {
  minInput?: number;
  maxInput?: number;
  minOutput?: number;
  maxOutput?: number;
};

export function scaleRadius(radius: number, options?: ScaleOptions): number {
  const minInput = options?.minInput || 0.001;
  const maxInput = options?.maxInput || 10;
  const minOutput = options?.minOutput || 1;
  const maxOutput = options?.maxOutput || 8;

  // Use log scaling since stellar radii have a large range
  const logInput = Math.log(Math.max(radius, minInput));
  const logMinInput = Math.log(minInput);
  const logMaxInput = Math.log(maxInput);

  // Calculate the scaled value using log scale
  const scaled = (logInput - logMinInput) / (logMaxInput - logMinInput);

  // Map to output range
  return minOutput + scaled * (maxOutput - minOutput);
}

// Predefined scale configurations for common use cases

// Configuration for planetary semi-major axis (in AU)
export const SEMI_MAJOR_AXIS_SCALE: ScaleOptions = {
  minInput: 0.1, // Closest planets like Mercury
  maxInput: 100.0, // Distant planets/objects
  minOutput: 10, // Minimum display distance
  maxOutput: 400, // Maximum display distance
};

// Configuration for planetary mass (in Earth masses)
export const PLANETARY_MASS_SCALE: ScaleOptions = {
  minInput: 0.000001, // Small asteroids/moons
  maxInput: 20.0, // Super-Jupiter sized planets
  minOutput: 2, // Minimum display radius
  maxOutput: 15, // Maximum display radius
};

// Configuration for stellar radii (in solar radii) - existing defaults
export const STELLAR_RADIUS_SCALE: ScaleOptions = {
  minInput: 0.001,
  maxInput: 650,
  minOutput: 1,
  maxOutput: 8,
};

// Convenience functions using predefined configurations
export function scaleSemiMajorAxis(
  semiMajorAxis: number,
  customOptions?: Partial<ScaleOptions>
): number {
  return scaleRadius(semiMajorAxis, {
    ...SEMI_MAJOR_AXIS_SCALE,
    ...customOptions,
  });
}

export function scalePlanetaryMass(
  mass: number,
  customOptions?: Partial<ScaleOptions>
): number {
  return scaleRadius(mass, { ...PLANETARY_MASS_SCALE, ...customOptions });
}

export function scaleStellarRadius(
  radius: number,
  customOptions?: Partial<ScaleOptions>
): number {
  return scaleRadius(radius, { ...STELLAR_RADIUS_SCALE, ...customOptions });
}

// Configuration for unified stellar and planetary mass scaling (in solar masses)
export const UNIFIED_MASS_SCALE: ScaleOptions = {
  minInput: 0.0000000066, // Pluto-sized objects (~0.0022 Earth masses)
  maxInput: 10.0, // Large stars in solar masses
  minOutput: 1, // Minimum display radius
  maxOutput: 40, // Maximum display radius
};

// Unified scaling function for both stellar and planetary masses
export function scaleUnifiedMass(
  mass: number,
  customOptions?: Partial<ScaleOptions>
): number {
  return scaleRadius(mass, { ...UNIFIED_MASS_SCALE, ...customOptions });
}
