import { renderDisk } from "./render-disk";
import "./style.css";

document.querySelector<HTMLDivElement>("#app")!.innerHTML = `
  <div class="h-screen w-screen fixed l-0 top-0">
    <h1 class="absolute left-4 top-4 z-10 text-4xl font-light text-neutral-100">Protodisk</h1>
    <canvas id="simulation" class="block fixed bg-neutral-900 top-0 left-0 h-screen w-screen" />
  </div>
`;

renderDisk(document.querySelector<HTMLCanvasElement>("#simulation")).catch(
  console.error
);
