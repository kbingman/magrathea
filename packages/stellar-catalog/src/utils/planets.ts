import type {
  Planet,
  PlanetarySystem,
  PlanetClass,
} from "@magrathea/planetary-wasm";

export function collectPlanets(systems: PlanetarySystem[]): Planet[] {
  return systems.reduce(
    (acc, { planets }) => [...acc, ...planets],
    [] as Planet[]
  );
}

export function filterByClass(
  planets: Planet[],
  planetClass: PlanetClass
): Planet[] {
  return planets.filter((planet) => planet.class === planetClass);
}

export function createPlanetClassDictionary(systems: PlanetarySystem[]) {
  const planets = collectPlanets(systems);

  return {
    compact: filterByClass(planets, "Compact"),
    transitional: filterByClass(planets, "Transitional"),
    volatile: filterByClass(planets, "Volatile"),
    giant: filterByClass(planets, "Giant"),
  };
}
