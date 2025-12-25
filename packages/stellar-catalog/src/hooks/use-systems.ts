import init, {
  generateSystemsBatch,
  type PlanetarySystem,
} from "@magrathea/planetary-wasm";
import { useEffect, useState } from "react";
// import { filterBySpectralType } from "../utils/stars";

export function useGenerateSystems(seed: string, count: number) {
  const [systems, setSystems] = useState<PlanetarySystem[]>([]);
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    init().then(() => setInitialized(true));
  }, []);

  useEffect(() => {
    if (initialized) {
      Promise.resolve(generateSystemsBatch(seed, count)).then((systems) => {
        setSystems(systems);
      });
    }
  }, [seed, count, initialized]);

  return {
    systems,
    initialized,
  };
}
