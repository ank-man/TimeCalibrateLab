/**
 * Highest Posterior Density (HPD) interval computation.
 * Also includes kernel density estimation for posterior visualization.
 */

/**
 * Compute the 95% HPD interval from sorted samples.
 * The HPD is the shortest interval containing 95% of the posterior mass.
 */
export function computeHPD(samples: number[], credMass: number = 0.95): [number, number] {
  const sorted = [...samples].sort((a, b) => a - b);
  const n = sorted.length;
  const intervalSize = Math.ceil(credMass * n);
  let minWidth = Infinity;
  let hpdMin = sorted[0];
  let hpdMax = sorted[intervalSize - 1];

  for (let i = 0; i <= n - intervalSize; i++) {
    const width = sorted[i + intervalSize - 1] - sorted[i];
    if (width < minWidth) {
      minWidth = width;
      hpdMin = sorted[i];
      hpdMax = sorted[i + intervalSize - 1];
    }
  }
  return [hpdMin, hpdMax];
}

/**
 * Compute posterior mean from samples.
 */
export function posteriorMean(samples: number[]): number {
  if (samples.length === 0) return 0;
  return samples.reduce((a, b) => a + b, 0) / samples.length;
}

/**
 * Compute posterior median from samples.
 */
export function posteriorMedian(samples: number[]): number {
  if (samples.length === 0) return 0;
  const sorted = [...samples].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0
    ? (sorted[mid - 1] + sorted[mid]) / 2
    : sorted[mid];
}

/**
 * Gaussian kernel density estimation.
 * Uses Silverman's rule of thumb for bandwidth selection.
 */
export function kernelDensityEstimate(
  samples: number[],
  numPoints: number = 200,
  bandwidth?: number
): { x: number; y: number }[] {
  if (samples.length === 0) return [];

  const n = samples.length;
  const sorted = [...samples].sort((a, b) => a - b);
  const min = sorted[0];
  const max = sorted[n - 1];

  if (min === max) {
    return [{ x: min, y: 1 }];
  }

  // Silverman's rule of thumb
  const sd = standardDeviation(samples);
  const iqr = sorted[Math.floor(0.75 * n)] - sorted[Math.floor(0.25 * n)];
  const h = bandwidth ?? 0.9 * Math.min(sd, iqr / 1.34) * Math.pow(n, -0.2);
  const effectiveH = Math.max(h, (max - min) / 100);

  const padding = 3 * effectiveH;
  const xMin = min - padding;
  const xMax = max + padding;
  const step = (xMax - xMin) / (numPoints - 1);

  const points: { x: number; y: number }[] = [];
  for (let i = 0; i < numPoints; i++) {
    const x = xMin + i * step;
    let density = 0;
    for (let j = 0; j < n; j++) {
      const z = (x - samples[j]) / effectiveH;
      density += Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
    }
    density /= n * effectiveH;
    points.push({ x, y: density });
  }
  return points;
}

/**
 * Standard deviation of samples.
 */
function standardDeviation(samples: number[]): number {
  const n = samples.length;
  if (n <= 1) return 0;
  const mean = samples.reduce((a, b) => a + b, 0) / n;
  const variance = samples.reduce((sum, x) => sum + (x - mean) ** 2, 0) / (n - 1);
  return Math.sqrt(variance);
}

/**
 * Effective sample size estimation using autocorrelation.
 */
export function effectiveSampleSize(samples: number[]): number {
  const n = samples.length;
  if (n < 10) return n;
  const mean = samples.reduce((a, b) => a + b, 0) / n;
  const variance = samples.reduce((sum, x) => sum + (x - mean) ** 2, 0) / n;
  if (variance === 0) return n;

  let sumAcf = 0;
  for (let lag = 1; lag < Math.min(n / 2, 100); lag++) {
    let acf = 0;
    for (let i = 0; i < n - lag; i++) {
      acf += (samples[i] - mean) * (samples[i + lag] - mean);
    }
    acf /= (n - lag) * variance;
    if (acf < 0.05) break;
    sumAcf += acf;
  }
  return Math.max(1, n / (1 + 2 * sumAcf));
}
