import { useGenerateSystems } from "./hooks/use-systems.ts";
import { DistanceLabel } from "./components/distance-label.tsx";
import { SectorInfo } from "./components/sector-info.tsx";
import { SystemCatalog } from "./components/systems-catalog.tsx";

const MAX_WIDTH = 720;

function App() {
  const { systems } = useGenerateSystems("gamma-quadrant", 1000);

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

      <div className="grid grid-cols-[1fr_420px] mt-12 bg-neutral-950">
        <SystemCatalog systems={systems} maxWidth={MAX_WIDTH} />
        <div className="border-amber-950 border p-4 fixed right-4 w-[420px] h-[calc(100vh-80px)]">
          <SectorInfo systems={systems} />
        </div>
      </div>
    </div>
  );
}

export default App;
