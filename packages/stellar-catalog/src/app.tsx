import { useGenerateSystems } from "./hooks/use-systems.ts";
import { SystemCard } from "./components/system-card.tsx";
import { DistanceLabel } from "./components/distance-label.tsx";
import { useMemo } from "react";
import { collectPlanets, filterByClass } from "./utils/planets.ts";
import { filterBySpectralType } from "./utils/stars.ts";

const MAX_WIDTH = 720;

function App() {
  const { systems } = useGenerateSystems("gamma-quadrant", 1000);

  const stars = useMemo(
    () => ({
      M: filterBySpectralType(systems, "M"),
      K: filterBySpectralType(systems, "K"),
      G: filterBySpectralType(systems, "G"),
      F: filterBySpectralType(systems, "F"),
      A: filterBySpectralType(systems, "A"),
    }),
    [systems]
  );

  const planets = useMemo(() => {
    const planets = collectPlanets(systems);

    return {
      compact: filterByClass(planets, "Compact"),
      transitional: filterByClass(planets, "Transitional"),
      volatile: filterByClass(planets, "Volatile"),
      giant: filterByClass(planets, "Giant"),
    };
  }, [systems]);

  return (
    <div className="p-4 pb-16 min-h-screen text-amber-600 bg-neutral-950 font-mono">
      <div className="fixed top-0 left-0 w-full bg-neutral-950 z-10 px-4 py-2 shadow-lg">
        <h1 className="text-lg font-semibold uppercase">
          Magrathea Systems Catalog
        </h1>
        <div className="relative h-3">
          <div className="text-xs text-neutral-500 absolute top-0">Name</div>
          <DistanceLabel semiMajorAxis={1} maxOutput={MAX_WIDTH}>
            1 AU
          </DistanceLabel>
          <DistanceLabel semiMajorAxis={10} maxOutput={MAX_WIDTH}>
            10 AU
          </DistanceLabel>
          <DistanceLabel semiMajorAxis={100} maxOutput={MAX_WIDTH}>
            100 AU
          </DistanceLabel>
          <DistanceLabel semiMajorAxis={1000} maxOutput={MAX_WIDTH}>
            1000 AU
          </DistanceLabel>
        </div>
      </div>

      <div className="grid grid-cols-[1fr_420px] mt-12">
        <div className="grid gap-2">
          {systems
            // .filter(({ stars }) => stars[0]?.spectralType === "G")
            .map(({ metadata, planets, stars }) => (
              <SystemCard
                key={`system-${metadata.id}`}
                id={metadata.id}
                star={stars[0]}
                planets={planets}
                name={metadata.catalogName || ""}
                maxWidth={MAX_WIDTH}
              />
            ))}
        </div>
        <div className="border-amber-950 border p-4 fixed right-4 w-[420px] h-[calc(100vh-80px)]">
          <h2 className="text-base uppercase">Stars</h2>
          <div className="text-xs uppercase mb-2">
            <div>M-Type: {stars.M.length}</div>
            <div>K-Type: {stars.K.length}</div>
            <div>G-Type: {stars.G.length}</div>
            <div>F-Type: {stars.F.length}</div>
            <div>A-Type: {stars.A.length}</div>
          </div>
          <h2 className="text-base uppercase">Planets</h2>
          <div className="text-xs uppercase mb-2">
            <div>Compact: {planets.compact.length}</div>
            <div>Transitional: {planets.transitional.length}</div>
            <div>Volatile: {planets.volatile.length}</div>
            <div>Giant: {planets.giant.length}</div>
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;
