/**
 * Prior distribution functions for Bayesian phylogenetic time calibration.
 * Supports Uniform, Lognormal, Exponential, and Gamma distributions.
 */

/** Uniform distribution density: p(x) = 1/(max-min) for min <= x <= max */
export function uniformPDF(x: number, min: number, max: number): number {
  if (x < min || x > max) return 0;
  return 1 / (max - min);
}

/** Uniform distribution log-density */
export function uniformLogPDF(x: number, min: number, max: number): number {
  if (x < min || x > max) return -Infinity;
  return -Math.log(max - min);
}

/** Sample from Uniform(min, max) */
export function uniformSample(min: number, max: number): number {
  return min + Math.random() * (max - min);
}

/** Lognormal PDF: p(x) = (1/(x*sigma*sqrt(2pi))) * exp(-(ln(x)-mu)^2 / (2*sigma^2)) */
export function lognormalPDF(x: number, mu: number, sigma: number): number {
  if (x <= 0) return 0;
  const logx = Math.log(x);
  const exponent = -((logx - mu) ** 2) / (2 * sigma ** 2);
  return (1 / (x * sigma * Math.sqrt(2 * Math.PI))) * Math.exp(exponent);
}

/** Lognormal log-density */
export function lognormalLogPDF(x: number, mu: number, sigma: number): number {
  if (x <= 0) return -Infinity;
  const logx = Math.log(x);
  return -Math.log(x) - Math.log(sigma) - 0.5 * Math.log(2 * Math.PI)
    - ((logx - mu) ** 2) / (2 * sigma ** 2);
}

/** Sample from LogNormal(mu, sigma) */
export function lognormalSample(mu: number, sigma: number): number {
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return Math.exp(mu + sigma * z);
}

/** Lognormal with offset: T = offset + LogNormal(mu, sigma) */
export function lognormalOffsetPDF(x: number, offset: number, mu: number, sigma: number): number {
  return lognormalPDF(x - offset, mu, sigma);
}

export function lognormalOffsetLogPDF(x: number, offset: number, mu: number, sigma: number): number {
  return lognormalLogPDF(x - offset, mu, sigma);
}

export function lognormalOffsetSample(offset: number, mu: number, sigma: number): number {
  return offset + lognormalSample(mu, sigma);
}

/** Exponential PDF: p(x) = lambda * exp(-lambda * x) for x >= 0 */
export function exponentialPDF(x: number, lambda: number): number {
  if (x < 0) return 0;
  return lambda * Math.exp(-lambda * x);
}

/** Exponential log-density */
export function exponentialLogPDF(x: number, lambda: number): number {
  if (x < 0) return -Infinity;
  return Math.log(lambda) - lambda * x;
}

/** Sample from Exponential(lambda) */
export function exponentialSample(lambda: number): number {
  return -Math.log(Math.random()) / lambda;
}

/** Exponential with offset: T = offset + Exp(lambda) */
export function exponentialOffsetPDF(x: number, offset: number, lambda: number): number {
  return exponentialPDF(x - offset, lambda);
}

export function exponentialOffsetLogPDF(x: number, offset: number, lambda: number): number {
  return exponentialLogPDF(x - offset, lambda);
}

export function exponentialOffsetSample(offset: number, lambda: number): number {
  return offset + exponentialSample(lambda);
}

/** Gamma PDF: p(x) = (beta^alpha / Gamma(alpha)) * x^(alpha-1) * exp(-beta*x) */
export function gammaPDF(x: number, alpha: number, beta: number): number {
  if (x <= 0) return 0;
  return Math.exp(gammaLogPDF(x, alpha, beta));
}

/** Gamma log-density */
export function gammaLogPDF(x: number, alpha: number, beta: number): number {
  if (x <= 0) return -Infinity;
  return alpha * Math.log(beta) - logGamma(alpha) + (alpha - 1) * Math.log(x) - beta * x;
}

/** Sample from Gamma(alpha, beta) using Marsaglia and Tsang's method */
export function gammaSample(alpha: number, beta: number): number {
  if (alpha < 1) {
    return gammaSample(alpha + 1, beta) * Math.pow(Math.random(), 1 / alpha);
  }
  const d = alpha - 1 / 3;
  const c = 1 / Math.sqrt(9 * d);
  while (true) {
    let x: number, v: number;
    do {
      x = normalSample(0, 1);
      v = 1 + c * x;
    } while (v <= 0);
    v = v * v * v;
    const u = Math.random();
    if (u < 1 - 0.0331 * (x * x) * (x * x)) return (d * v) / beta;
    if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return (d * v) / beta;
  }
}

/** Standard normal sample (Box-Muller) */
export function normalSample(mu: number, sigma: number): number {
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return mu + sigma * z;
}

/** Normal PDF */
export function normalPDF(x: number, mu: number, sigma: number): number {
  const z = (x - mu) / sigma;
  return (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
}

/** Normal log-PDF */
export function normalLogPDF(x: number, mu: number, sigma: number): number {
  const z = (x - mu) / sigma;
  return -0.5 * Math.log(2 * Math.PI) - Math.log(sigma) - 0.5 * z * z;
}

/** Log-gamma function (Lanczos approximation) */
export function logGamma(x: number): number {
  const g = 7;
  const coefficients = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
  ];
  if (x < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
  }
  x -= 1;
  let a = coefficients[0];
  const t = x + g + 0.5;
  for (let i = 1; i < g + 2; i++) {
    a += coefficients[i] / (x + i);
  }
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

export type CalibrationDistribution =
  | { type: 'uniform'; min: number; max: number }
  | { type: 'lognormal'; offset: number; mu: number; sigma: number }
  | { type: 'exponential'; offset: number; lambda: number };

/** Evaluate the log-prior density for a calibration distribution at value x */
export function calibrationLogPDF(x: number, dist: CalibrationDistribution): number {
  switch (dist.type) {
    case 'uniform':
      return uniformLogPDF(x, dist.min, dist.max);
    case 'lognormal':
      return lognormalOffsetLogPDF(x, dist.offset, dist.mu, dist.sigma);
    case 'exponential':
      return exponentialOffsetLogPDF(x, dist.offset, dist.lambda);
  }
}

/** Evaluate the prior density for a calibration distribution at value x */
export function calibrationPDF(x: number, dist: CalibrationDistribution): number {
  switch (dist.type) {
    case 'uniform':
      return uniformPDF(x, dist.min, dist.max);
    case 'lognormal':
      return lognormalOffsetPDF(x, dist.offset, dist.mu, dist.sigma);
    case 'exponential':
      return exponentialOffsetPDF(x, dist.offset, dist.lambda);
  }
}

/** Sample from a calibration distribution */
export function calibrationSample(dist: CalibrationDistribution): number {
  switch (dist.type) {
    case 'uniform':
      return uniformSample(dist.min, dist.max);
    case 'lognormal':
      return lognormalOffsetSample(dist.offset, dist.mu, dist.sigma);
    case 'exponential':
      return exponentialOffsetSample(dist.offset, dist.lambda);
  }
}

/** Generate density curve data points for a calibration distribution */
export function calibrationDensityCurve(
  dist: CalibrationDistribution,
  numPoints: number = 200
): { x: number; y: number }[] {
  let xMin: number, xMax: number;
  switch (dist.type) {
    case 'uniform':
      xMin = dist.min - (dist.max - dist.min) * 0.1;
      xMax = dist.max + (dist.max - dist.min) * 0.1;
      break;
    case 'lognormal': {
      const mean = dist.offset + Math.exp(dist.mu + dist.sigma ** 2 / 2);
      const sd = Math.sqrt((Math.exp(dist.sigma ** 2) - 1) * Math.exp(2 * dist.mu + dist.sigma ** 2));
      xMin = Math.max(dist.offset, mean - 3 * sd);
      xMax = mean + 3 * sd;
      break;
    }
    case 'exponential': {
      const expMean = dist.offset + 1 / dist.lambda;
      xMin = dist.offset;
      xMax = expMean + 4 / dist.lambda;
      break;
    }
  }
  const points: { x: number; y: number }[] = [];
  const step = (xMax - xMin) / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    const x = xMin + i * step;
    points.push({ x, y: calibrationPDF(x, dist) });
  }
  return points;
}
