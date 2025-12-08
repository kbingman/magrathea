export const SOLAR_TO_EARTH_MASS = 332_946.05;
export const SOLAR_TO_JUPITER_MASS = 1_047.35;
export const AU_TO_KM = 1.496e8;

export const toEarthMasses = (solar: number) => solar * SOLAR_TO_EARTH_MASS;
export const toJupiterMasses = (solar: number) => solar * SOLAR_TO_JUPITER_MASS;
