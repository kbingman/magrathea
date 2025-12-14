import type {
  Planet,
  PlanetarySystem,
  PlanetClass,
} from "@magrathea/planetary-wasm";

export function collectPlanets(systems: PlanetarySystem[]): Planet[] {
  return systems
    .reduce((acc, { planets }) => [...acc, ...planets], [] as Planet[])
    .filter(
      ({ planetType }) =>
        planetType.type !== "kuiperBeltObject" &&
        planetType.type !== "dwarfPlanet"
    );
}

export function filterByClass(
  planets: Planet[],
  planetClass: PlanetClass
): Planet[] {
  return planets.filter((planet) => planet.class === planetClass);
}
