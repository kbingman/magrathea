/**
 * Takes a Canvas element and resizes it according to the devicePixelRatio,
 * then uses CSS to reduce the size back to the original. This creates a
 * Canvas element with Retina resolution.
 *
 * @param canvas
 * @param width
 * @param height
 * @param pixelRatio
 *
 * @returns void
 */
export function resizeCanvas(
  canvas: HTMLCanvasElement,
  width: number,
  height: number,
  pixelRatio: number = 1
) {
  canvas.width = width * pixelRatio;
  canvas.height = height * pixelRatio;
  canvas.style.height = `${height}px`;
  canvas.style.width = `${width}px`;
  canvas.style.display = "block";
}

/**
 * Returns a properly sized Retina ready Canvas element
 *
 * @param width
 * @param height
 * @param devicePixelRatio
 *
 * @returns HTMLCanvasElement
 */
export function createCanvas(
  width: number,
  height: number,
  devicePixelRatio: number
) {
  const canvas = document.createElement("CANVAS") as HTMLCanvasElement;
  resizeCanvas(canvas, width, height, devicePixelRatio);

  return canvas;
}

type BaseOptions = {
  fill?: string | CanvasGradient | CanvasPattern;
  lineWidth?: number;
  strokeStyle?: string | CanvasGradient | CanvasPattern;
};

type CircleOptions = {
  x: number;
  y: number;
  r: number;
} & BaseOptions;

type LineOptions = {
  x1: number;
  y1: number;
  x2: number;
  y2: number;
  strokeStyle?: string | CanvasGradient | CanvasPattern;
  lineWidth?: number;
  opacity?: number;
};

/**
 * @param ctx
 * @param options
 *
 * @returns void
 */
export function drawCircle(
  ctx: CanvasRenderingContext2D,
  { x, y, r, fill, strokeStyle, lineWidth }: CircleOptions
) {
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2 * Math.PI);
  if (fill) {
    ctx.fillStyle = fill;
    ctx.fill();
  }
  if (lineWidth) {
    ctx.lineWidth = lineWidth;
  }
  if (strokeStyle) {
    ctx.strokeStyle = strokeStyle;
    ctx.stroke();
  }
}

/**
 * @param ctx
 * @param options
 *
 * @returns void
 */
export function drawLine(
  ctx: CanvasRenderingContext2D,
  { x1, y1, x2, y2, strokeStyle, lineWidth, opacity }: LineOptions
) {
  ctx.save();

  if (opacity !== undefined) {
    ctx.globalAlpha = opacity;
  }

  if (lineWidth) {
    ctx.lineWidth = lineWidth;
  }

  if (strokeStyle) {
    ctx.strokeStyle = strokeStyle;
  }

  ctx.beginPath();
  ctx.moveTo(x1, y1);
  ctx.lineTo(x2, y2);
  ctx.stroke();

  ctx.restore();
}

type RectOptions = {
  x: number;
  y: number;
  w: number;
  h: number;
} & BaseOptions;

/**
 * @param ctx
 * @param options
 *
 * @returns void
 */
export function drawRect(
  ctx: CanvasRenderingContext2D,
  { x, y, w, h, fill, strokeStyle, lineWidth }: RectOptions
) {
  if (fill) {
    ctx.fillStyle = fill;
    ctx.fillRect(x, y, w, h);
  }
  if (lineWidth) {
    ctx.lineWidth = lineWidth;
  }
  if (strokeStyle) {
    ctx.strokeStyle = strokeStyle;
    ctx.strokeRect(x, y, w, h);
  }
}

/**
 * Fill the entire canvas with the given color
 *
 * @param ctx
 * @param fill
 *
 * @returns void
 */
export function fillCanvas(
  ctx: CanvasRenderingContext2D,
  fill: string | CanvasGradient | CanvasPattern
) {
  drawRect(ctx, {
    x: 0,
    y: 0,
    w: ctx.canvas.width,
    h: ctx.canvas.height,
    fill,
  });
}
