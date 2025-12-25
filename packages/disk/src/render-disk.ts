import init, {
  grid_disk_create,
  grid_disk_evolve,
  grid_disk_profile,
  solar_analog,
} from "@magrathea/magrathea-wasm";

interface Star {
  color: { r: number; g: number; b: number };
  radius_solar: number;
  mass_solar: number;
}

interface DiskProfile {
  radiiAu: number[];
  sigmaGCm2: number[];
  innerAu: number;
  outerAu: number;
  totalMassSolar: number;
}

interface SimulationState {
  diskId: number;
  star: Star;
  elapsedYears: number;
  animationId: number | null;
  initialMaxSigma: number; // Fixed reference for normalization
}

let state: SimulationState | null = null;

export async function renderDisk(canvas: HTMLCanvasElement | null) {
  if (!canvas) return;

  await init();

  const ctx = canvas.getContext("2d");
  if (!ctx) return;

  // Create star and disk
  const star = solar_analog() as Star;
  const diskId = grid_disk_create(star, 200);
  const initialProfile = grid_disk_profile(diskId) as DiskProfile;

  state = {
    diskId,
    star,
    elapsedYears: 0,
    animationId: null,
    initialMaxSigma: Math.max(...initialProfile.sigmaGCm2),
  };

  // Set canvas size to match display
  const resize = () => {
    canvas.width = window.innerWidth * window.devicePixelRatio;
    canvas.height = window.innerHeight * window.devicePixelRatio;
    ctx.setTransform(
      window.devicePixelRatio,
      0,
      0,
      window.devicePixelRatio,
      0,
      0
    );
    draw(
      ctx,
      canvas.width / window.devicePixelRatio,
      canvas.height / window.devicePixelRatio
    );
  };

  window.addEventListener("resize", resize);
  resize();

  // Start animation
  animate(ctx, canvas);
}

function animate(ctx: CanvasRenderingContext2D, canvas: HTMLCanvasElement) {
  if (!state) return;

  // Evolve many steps per frame to see visible changes
  // Viscous timescale is ~10^5-10^6 years, so we need lots of steps
  const stepsPerFrame = 10000;
  const elapsedThisFrame = grid_disk_evolve(state.diskId, stepsPerFrame);
  state.elapsedYears += elapsedThisFrame;

  const width = canvas.width / window.devicePixelRatio;
  const height = canvas.height / window.devicePixelRatio;
  draw(ctx, width, height);

  state.animationId = requestAnimationFrame(() => animate(ctx, canvas));
}

function draw(ctx: CanvasRenderingContext2D, width: number, height: number) {
  if (!state) return;

  const profile = grid_disk_profile(state.diskId) as DiskProfile;
  const centerX = width / 2;
  const centerY = height / 2;

  // Clear canvas
  ctx.fillStyle = "#171717"; // neutral-900
  ctx.fillRect(0, 0, width, height);

  // Draw disk
  drawDisk(ctx, centerX, centerY, profile);

  // Draw star
  drawStar(ctx, centerX, centerY, state.star);

  // Draw info overlay
  drawInfo(ctx, profile, state.elapsedYears);
}

function drawDisk(
  ctx: CanvasRenderingContext2D,
  cx: number,
  cy: number,
  profile: DiskProfile
) {
  if (!state) return;

  const { radiiAu, sigmaGCm2 } = profile;

  // Use current min/max for good contrast
  const maxSigma = Math.max(...sigmaGCm2);
  const minSigma = Math.min(...sigmaGCm2);

  // Use log scale for radius to handle the large dynamic range (0.1 to 100 AU)
  const logInner = Math.log10(profile.innerAu);
  const logOuter = Math.log10(profile.outerAu);
  const logRange = logOuter - logInner;

  const maxRadius =
    (Math.min(ctx.canvas.width, ctx.canvas.height) * 0.4) /
    window.devicePixelRatio;

  const auToPixel = (au: number) => {
    const logR = Math.log10(au);
    const t = (logR - logInner) / logRange;
    return t * maxRadius;
  };

  const outerPx = auToPixel(profile.outerAu);

  // Create radial gradient using grid zones directly
  const gradient = ctx.createRadialGradient(cx, cy, 0, cx, cy, outerPx);

  // Transparent center
  gradient.addColorStop(0, "rgba(200, 80, 60, 0)");

  // Sample grid zones for smooth gradient (every 10th point)
  for (let i = 0; i < radiiAu.length; i += 10) {
    const rPx = auToPixel(radiiAu[i]);
    const gradientPos = rPx / outerPx;
    const sigma = sigmaGCm2[i];

    const normalized = (sigma - minSigma) / (maxSigma - minSigma);
    const alpha = 0.05 + normalized * 0.45;

    gradient.addColorStop(
      gradientPos,
      `rgba(200, 80, 60, ${alpha.toFixed(3)})`
    );
  }
  // Ensure we include the last point
  const lastIdx = radiiAu.length - 1;
  const lastNormalized =
    (sigmaGCm2[lastIdx] - minSigma) / (maxSigma - minSigma);
  const lastAlpha = 0.05 + lastNormalized * 0.45;
  gradient.addColorStop(1, `rgba(200, 80, 60, ${lastAlpha.toFixed(3)})`);

  // Draw the disk
  ctx.beginPath();
  ctx.arc(cx, cy, outerPx, 0, Math.PI * 2);
  ctx.fillStyle = gradient;
  ctx.fill();

  // Draw inner boundary
  ctx.strokeStyle = "rgba(255, 150, 100, 0.3)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.arc(cx, cy, auToPixel(profile.innerAu), 0, Math.PI * 2);
  ctx.stroke();
}

function drawStar(
  ctx: CanvasRenderingContext2D,
  cx: number,
  cy: number,
  star: Star
) {
  const { r, g, b } = star.color;

  // Fixed visible size (realistic scale would be invisible)
  const starRadiusPx = 2;

  // Draw glow
  const gradient = ctx.createRadialGradient(
    cx,
    cy,
    0,
    cx,
    cy,
    starRadiusPx * 3
  );
  gradient.addColorStop(0, `rgba(${r}, ${g}, ${b}, 0.8)`);
  gradient.addColorStop(0.3, `rgba(${r}, ${g}, ${b}, 0.3)`);
  gradient.addColorStop(1, `rgba(${r}, ${g}, ${b}, 0)`);

  ctx.fillStyle = gradient;
  ctx.beginPath();
  ctx.arc(cx, cy, starRadiusPx * 3, 0, Math.PI * 2);
  ctx.fill();

  // Draw star core
  ctx.fillStyle = `rgb(${r}, ${g}, ${b})`;
  ctx.beginPath();
  ctx.arc(cx, cy, starRadiusPx, 0, Math.PI * 2);
  ctx.fill();
}

function drawInfo(
  ctx: CanvasRenderingContext2D,
  profile: DiskProfile,
  elapsedYears: number
) {
  ctx.fillStyle = "#f5f5f5";
  ctx.font = "14px ui-monospace, monospace";

  const formatYears = (years: number) => {
    if (years < 1000) return `${years.toFixed(0)} yr`;
    if (years < 1e6) return `${(years / 1000).toFixed(1)} kyr`;
    return `${(years / 1e6).toFixed(2)} Myr`;
  };

  const lines = [
    `Time: ${formatYears(elapsedYears)}`,
    `Disk mass: ${(profile.totalMassSolar * 100).toFixed(2)}% M_star`,
    `Extent: ${profile.innerAu.toFixed(2)} - ${profile.outerAu.toFixed(1)} AU`,
  ];

  lines.forEach((line, i) => {
    ctx.fillText(line, 20, 80 + i * 20);
  });
}
