import type {
  Planet,
  PlanetarySystem,
  PlanetClass,
} from "@magrathea/planetary-wasm";

export function collectPlanets(systems: PlanetarySystem[]): Planet[] {
  return systems
    .flatMap(({ planets }) => planets)
    .filter(
      ({ planetType }) =>
        planetType.type !== "kuiperBeltObject" &&
        planetType.type !== "dwarfPlanet",
    );
}

export function filterByClass(
  planets: Planet[],
  planetClass: PlanetClass,
): Planet[] {
  return planets.filter((planet) => planet.class === planetClass);
}
