import { useMemo } from "react";

import type { PlanetarySystem } from "@magrathea/planetary-wasm";
import { collectPlanets, filterByClass } from "../utils/planets";

export function usePlanets(systems: PlanetarySystem[]) {
  return useMemo(() => {
    const planets = collectPlanets(systems);

    return {
      compact: filterByClass(planets, "Compact"),
      transitional: filterByClass(planets, "Transitional"),
      volatile: filterByClass(planets, "Volatile"),
      giant: filterByClass(planets, "Giant"),
    };
  }, [systems]);
}
