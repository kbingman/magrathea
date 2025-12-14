import { useMemo } from "react";

import type { PlanetarySystem } from "@magrathea/planetary-wasm";
import { filterBySpectralType } from "../utils/stars";

export function useStars(systems: PlanetarySystem[]) {
  return useMemo(
    () => ({
      M: filterBySpectralType(systems, "M"),
      K: filterBySpectralType(systems, "K"),
      G: filterBySpectralType(systems, "G"),
      F: filterBySpectralType(systems, "F"),
      A: filterBySpectralType(systems, "A"),
    }),
    [systems]
  );
}
