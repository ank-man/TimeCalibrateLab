/**
 * Likelihood computation for phylogenetic time calibration.
 * Implements JC69 model and distance-based approximation.
 */

/**
 * JC69 probability of observing the same base after distance d.
 * P(same) = 1/4 + 3/4 * exp(-4d/3)
 */
export function jc69ProbSame(d: number): number {
  if (d <= 0) return 1;
  return 0.25 + 0.75 * Math.exp((-4 * d) / 3);
}

/**
 * JC69 probability of observing a different base after distance d.
 * P(diff) = 1/4 - 1/4 * exp(-4d/3)
 */
export function jc69ProbDiff(d: number): number {
  if (d <= 0) return 0;
  return 0.25 - 0.25 * Math.exp((-4 * d) / 3);
}

/**
 * JC69 log-likelihood for a pairwise comparison.
 * Given observed proportion of differences p_obs and alignment length L,
 * compute log-likelihood of distance d.
 *
 * L(d) = nSame * log(P(same|d)) + nDiff * log(P(diff|d))
 */
export function jc69LogLikelihood(d: number, pObs: number, L: number): number {
  const nDiff = Math.round(pObs * L);
  const nSame = L - nDiff;
  const pSame = jc69ProbSame(d);
  const pDiffModel = jc69ProbDiff(d);

  if (pSame <= 0 || pDiffModel <= 0) return -Infinity;
  return nSame * Math.log(pSame) + nDiff * Math.log(pDiffModel);
}

/**
 * Distance-based approximation log-likelihood.
 * observed distance ~ Normal(expected, variance/L)
 * This is a fast Gaussian approximation for the likelihood.
 */
export function distanceLogLikelihood(
  dExpected: number,
  dObserved: number,
  L: number
): number {
  const variance = dExpected / L;
  if (variance <= 0) return dExpected === dObserved ? 0 : -Infinity;
  const diff = dObserved - dExpected;
  return -0.5 * Math.log(2 * Math.PI * variance) - (diff * diff) / (2 * variance);
}

/**
 * Simulate observed pairwise distance given true distance and alignment length.
 * Uses binomial sampling under JC69.
 */
export function simulateObservedDistance(trueDistance: number, L: number): number {
  const pDiff = jc69ProbDiff(trueDistance);
  let nDiff = 0;
  for (let i = 0; i < L; i++) {
    if (Math.random() < pDiff) nDiff++;
  }
  return nDiff / L;
}

/**
 * Fast simulation using normal approximation of binomial.
 */
export function simulateObservedDistanceFast(trueDistance: number, L: number): number {
  const pDiff = jc69ProbDiff(trueDistance);
  const mean = pDiff;
  const sd = Math.sqrt((pDiff * (1 - pDiff)) / L);
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return Math.max(0, Math.min(0.74, mean + sd * z));
}

/**
 * Compute expected substitution distance: d = rate * time
 */
export function expectedDistance(rate: number, time: number): number {
  return rate * time;
}

/**
 * Compute log-likelihood for a branch given:
 * - nodeAge: age of child node (or 0 for tips)
 * - parentAge: age of parent node
 * - rate: substitution rate for this branch
 * - observedDist: observed pairwise distance for this branch
 * - alignmentLength: number of sites
 * - mode: 'jc69' or 'distance'
 */
export function branchLogLikelihood(
  nodeAge: number,
  parentAge: number,
  rate: number,
  observedDist: number,
  alignmentLength: number,
  mode: 'jc69' | 'distance' = 'jc69'
): number {
  const branchTime = parentAge - nodeAge;
  if (branchTime < 0) return -Infinity;
  const d = expectedDistance(rate, branchTime);
  if (mode === 'jc69') {
    return jc69LogLikelihood(d, observedDist, alignmentLength);
  } else {
    return distanceLogLikelihood(d, observedDist, alignmentLength);
  }
}

/**
 * Generate likelihood curve data points for a node age.
 */
export function likelihoodCurve(
  parentAge: number,
  childAge: number,
  rate: number,
  observedDist: number,
  alignmentLength: number,
  mode: 'jc69' | 'distance',
  numPoints: number = 200
): { x: number; y: number }[] {
  const xMin = childAge;
  const xMax = parentAge;
  if (xMax <= xMin) return [];
  const points: { x: number; y: number }[] = [];
  const step = (xMax - xMin) / (numPoints - 1);

  let maxLL = -Infinity;
  const rawPoints: { x: number; ll: number }[] = [];

  for (let i = 0; i < numPoints; i++) {
    const age = xMin + i * step;
    const branchTime = parentAge - age;
    const d = expectedDistance(rate, branchTime);
    const ll = mode === 'jc69'
      ? jc69LogLikelihood(d, observedDist, alignmentLength)
      : distanceLogLikelihood(d, observedDist, alignmentLength);
    rawPoints.push({ x: age, ll });
    if (ll > maxLL) maxLL = ll;
  }

  for (const p of rawPoints) {
    points.push({ x: p.x, y: Math.exp(p.ll - maxLL) });
  }
  return points;
}
