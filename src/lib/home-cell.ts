type Point = [number, number];
type Curve = [Point, Point, Point];
type Contour = { start: Point; curves: Curve[] };

const clamp = (value: number) => Math.min(1, Math.max(0, value));
const ease = (value: number) => {
  const t = clamp(value);
  return t * t * (3 - 2 * t);
};
const mix = (a: Point, b: Point, t: number): Point => [
  a[0] + (b[0] - a[0]) * t,
  a[1] + (b[1] - a[1]) * t,
];

function arcs(cx: number, cy: number, rx: number, ry: number, angles: number[]): Contour {
  const point = (angle: number): Point => [cx + rx * Math.cos(angle), cy + ry * Math.sin(angle)];
  const radians = angles.map((angle) => angle * Math.PI / 180);
  return {
    start: point(radians[0]),
    curves: radians.slice(1).map((end, index): Curve => {
      const start = radians[index];
      const k = 4 / 3 * Math.tan((end - start) / 4);
      const a = point(start);
      const b = point(end);
      return [
        [a[0] - k * rx * Math.sin(start), a[1] + k * ry * Math.cos(start)],
        [b[0] + k * rx * Math.sin(end), b[1] - k * ry * Math.cos(end)],
        b,
      ];
    }),
  };
}

const circle = arcs(160, 160, 158, 158, [0, -45, -90, -180, -270, -315, -360]);
// Begin with the same circular outer ends as the daughter cells and an open middle.
const openCell: Contour = {
  start: [360, 2],
  curves: [
    [[308, 2], [301.354, 2], [271.723, 2]],
    [[242.092, 2], [201.904, 2], [160, 2]],
    ...circle.curves.slice(2, 4),
    [[201.904, 318], [242.092, 318], [271.723, 318]],
    [[301.354, 318], [308, 318], [360, 318]],
  ],
};

function joined(neck: number): Contour {
  // The outer arcs remain circular; equal tangent handles keep the shoulders smooth.
  const upper = circle.curves[0][2];
  const lower = circle.curves[4][2];
  const upperHandle = circle.curves[1][0];
  const lowerHandle = circle.curves[4][1];
  return {
    start: [360, 160 - neck],
    curves: [
      [[308, 160 - neck], [2 * upper[0] - upperHandle[0], 2 * upper[1] - upperHandle[1]], upper],
      ...circle.curves.slice(1, 5),
      [[2 * lower[0] - lowerHandle[0], 2 * lower[1] - lowerHandle[1]], [308, 160 + neck], [360, 160 + neck]],
    ],
  };
}

function interpolate(a: Contour, b: Contour, t: number): Contour {
  return {
    start: mix(a.start, b.start, t),
    curves: a.curves.map((curve, index) => curve.map((point, j) => mix(point, b.curves[index][j], t)) as Curve),
  };
}

function reflectAndReverse(contour: Contour): Contour {
  const reflect = ([x, y]: Point): Point => [720 - x, y];
  const starts = [contour.start, ...contour.curves.map((curve) => curve[2])];
  return {
    start: reflect(starts[starts.length - 1]),
    curves: contour.curves.map((curve, index): Curve => [reflect(curve[1]), reflect(curve[0]), reflect(starts[index])]).reverse(),
  };
}

const pointString = (point: Point) => point.map((value) => Number(value.toFixed(3))).join(' ');
const curveString = (contour: Contour) => contour.curves.map((curve) => `C ${curve.map(pointString).join(' ')}`).join(' ');

export function cellPath(progress: number): string {
  const t = clamp(progress);
  let left: Contour;
  if (t < 0.48) left = interpolate(openCell, joined(23), ease(t / 0.48));
  else if (t < 0.72) left = joined(23 * (1 - ease((t - 0.48) / 0.24)));
  else left = interpolate(joined(0), circle, ease((t - 0.72) / 0.28));

  const right = reflectAndReverse(left);
  // At zero neck width the shared contour becomes two closed daughter cells.
  return `M ${pointString(left.start)} ${curveString(left)} ${t >= 0.72 ? `Z M ${pointString(right.start)}` : ''} ${curveString(right)} Z`;
}
