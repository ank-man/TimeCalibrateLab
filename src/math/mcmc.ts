/**
 * Lightweight Metropolis-Hastings MCMC engine for Bayesian time calibration.
 * Runs entirely in browser; designed for educational demonstration.
 */

import { calibrationLogPDF, type CalibrationDistribution, gammaLogPDF, lognormalLogPDF } from './priors';
import { branchLogLikelihood, simulateObservedDistanceFast } from './likelihood';

export interface TreeNode {
  id: string;
  name: string;
  parent: string | null;
  children: string[];
  isTip: boolean;
  trueAge: number;
  age: number;
  calibration: CalibrationDistribution | null;
  observedDistances: Record<string, number>;
}

export interface MCMCConfig {
  iterations: number;
  burnin: number;
  proposalSD: number;
  rateProposalSD: number;
  clockType: 'strict' | 'relaxed';
  rateAlpha: number;
  rateBeta: number;
  relaxedSigma: number;
  globalRate: number;
  alignmentLength: number;
  likelihoodMode: 'jc69' | 'distance';
}

export interface MCMCSamples {
  nodeAges: Record<string, number[]>;
  rates: Record<string, number[]>;
  globalRate: number[];
  logPosterior: number[];
}

export interface MCMCResult {
  samples: MCMCSamples;
  acceptanceRate: number;
}

const DEFAULT_CONFIG: MCMCConfig = {
  iterations: 3000,
  burnin: 750,
  proposalSD: 2.0,
  rateProposalSD: 0.005,
  clockType: 'strict',
  rateAlpha: 2.0,
  rateBeta: 200.0,
  relaxedSigma: 0.3,
  globalRate: 0.01,
  alignmentLength: 1000,
  likelihoodMode: 'jc69',
};

/**
 * Compute log-prior for all node ages given calibrations.
 */
function computeLogPrior(
  nodes: Record<string, TreeNode>,
  nodeAges: Record<string, number>
): number {
  let logPrior = 0;
  for (const id of Object.keys(nodes)) {
    const node = nodes[id];
    if (node.isTip) continue;

    // Parent must be older than children
    for (const childId of node.children) {
      const childAge = nodes[childId].isTip ? 0 : (nodeAges[childId] ?? 0);
      if (nodeAges[id] <= childAge) return -Infinity;
    }

    // If node has parent, parent must be older
    if (node.parent) {
      const parentAge = nodeAges[node.parent];
      if (parentAge !== undefined && parentAge <= nodeAges[id]) return -Infinity;
    }

    // Calibration prior
    if (node.calibration) {
      logPrior += calibrationLogPDF(nodeAges[id], node.calibration);
    }
  }
  return logPrior;
}

/**
 * Compute log-prior for rate parameters.
 */
function computeRateLogPrior(
  globalRate: number,
  branchRates: Record<string, number>,
  config: MCMCConfig
): number {
  let logPrior = gammaLogPDF(globalRate, config.rateAlpha, config.rateBeta);
  if (config.clockType === 'relaxed') {
    const logMean = Math.log(globalRate) - 0.5 * config.relaxedSigma ** 2;
    for (const rate of Object.values(branchRates)) {
      logPrior += lognormalLogPDF(rate, logMean, config.relaxedSigma);
    }
  }
  return logPrior;
}

/**
 * Compute total log-likelihood across all branches.
 */
function computeLogLikelihood(
  nodes: Record<string, TreeNode>,
  nodeAges: Record<string, number>,
  branchRates: Record<string, number>,
  globalRate: number,
  config: MCMCConfig
): number {
  let logLik = 0;
  for (const id of Object.keys(nodes)) {
    const node = nodes[id];
    if (!node.parent) continue;

    const childAge = node.isTip ? 0 : (nodeAges[id] ?? 0);
    const parentAge = nodeAges[node.parent] ?? 0;
    const rate = config.clockType === 'strict' ? globalRate : (branchRates[id] ?? globalRate);

    const obsDist = node.observedDistances[node.parent] ?? 0;
    const ll = branchLogLikelihood(
      childAge, parentAge, rate, obsDist,
      config.alignmentLength, config.likelihoodMode
    );
    if (!isFinite(ll)) return -Infinity;
    logLik += ll;
  }
  return logLik;
}

/**
 * Run Metropolis-Hastings MCMC for time calibration.
 */
export function runMCMC(
  nodes: Record<string, TreeNode>,
  config: Partial<MCMCConfig> = {},
  onProgress?: (pct: number) => void
): MCMCResult {
  const cfg = { ...DEFAULT_CONFIG, ...config };

  // Initialize current state
  const currentAges: Record<string, number> = {};
  const branchRates: Record<string, number> = {};
  const internalNodeIds: string[] = [];

  for (const id of Object.keys(nodes)) {
    const node = nodes[id];
    if (!node.isTip) {
      currentAges[id] = node.age;
      internalNodeIds.push(id);
    }
    if (node.parent) {
      branchRates[id] = cfg.globalRate;
    }
  }

  let currentRate = cfg.globalRate;

  // Samples storage
  const samples: MCMCSamples = {
    nodeAges: {},
    rates: {},
    globalRate: [],
    logPosterior: [],
  };
  for (const id of internalNodeIds) {
    samples.nodeAges[id] = [];
  }
  for (const id of Object.keys(branchRates)) {
    samples.rates[id] = [];
  }

  let currentLogPrior = computeLogPrior(nodes, currentAges);
  let currentRateLogPrior = computeRateLogPrior(currentRate, branchRates, cfg);
  let currentLogLik = computeLogLikelihood(nodes, currentAges, branchRates, currentRate, cfg);
  let currentLogPost = currentLogPrior + currentRateLogPrior + currentLogLik;

  let accepted = 0;
  const totalProposals = cfg.iterations;

  for (let iter = 0; iter < cfg.iterations; iter++) {
    // Propose new node ages
    const proposedAges = { ...currentAges };
    const nodeToUpdate = internalNodeIds[Math.floor(Math.random() * internalNodeIds.length)];
    proposedAges[nodeToUpdate] += (Math.random() - 0.5) * 2 * cfg.proposalSD;

    // Ensure positive
    if (proposedAges[nodeToUpdate] <= 0) {
      proposedAges[nodeToUpdate] = currentAges[nodeToUpdate];
    }

    // Propose rate changes
    let proposedRate = currentRate;
    const proposedBranchRates = { ...branchRates };

    if (Math.random() < 0.3) {
      proposedRate = currentRate + (Math.random() - 0.5) * 2 * cfg.rateProposalSD;
      if (proposedRate <= 0) proposedRate = currentRate;
    }

    if (cfg.clockType === 'relaxed' && Math.random() < 0.3) {
      const branchIds = Object.keys(proposedBranchRates);
      const branchToUpdate = branchIds[Math.floor(Math.random() * branchIds.length)];
      proposedBranchRates[branchToUpdate] += (Math.random() - 0.5) * 2 * cfg.rateProposalSD;
      if (proposedBranchRates[branchToUpdate] <= 0) {
        proposedBranchRates[branchToUpdate] = branchRates[branchToUpdate];
      }
    }

    // Compute proposed posterior
    const proposedLogPrior = computeLogPrior(nodes, proposedAges);
    const proposedRateLogPrior = computeRateLogPrior(proposedRate, proposedBranchRates, cfg);
    const proposedLogLik = computeLogLikelihood(nodes, proposedAges, proposedBranchRates, proposedRate, cfg);
    const proposedLogPost = proposedLogPrior + proposedRateLogPrior + proposedLogLik;

    // Acceptance ratio
    const logAlpha = proposedLogPost - currentLogPost;
    if (Math.log(Math.random()) < logAlpha) {
      // Accept
      Object.assign(currentAges, proposedAges);
      Object.assign(branchRates, proposedBranchRates);
      currentRate = proposedRate;
      currentLogPrior = proposedLogPrior;
      currentRateLogPrior = proposedRateLogPrior;
      currentLogLik = proposedLogLik;
      currentLogPost = proposedLogPost;
      accepted++;
    }

    // Store samples (after burnin)
    if (iter >= cfg.burnin) {
      for (const id of internalNodeIds) {
        samples.nodeAges[id].push(currentAges[id]);
      }
      for (const id of Object.keys(branchRates)) {
        samples.rates[id].push(branchRates[id]);
      }
      samples.globalRate.push(currentRate);
      samples.logPosterior.push(currentLogPost);
    }

    // Progress callback
    if (onProgress && iter % 100 === 0) {
      onProgress(iter / cfg.iterations);
    }
  }

  return {
    samples,
    acceptanceRate: accepted / totalProposals,
  };
}

/**
 * Create default example tree with 8 taxa.
 */
export function createExampleTree(): Record<string, TreeNode> {
  const nodes: Record<string, TreeNode> = {};

  // Tips (age = 0, extant)
  const tipNames = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon', 'OWMonkey', 'NWMonkey', 'Lemur'];
  for (const name of tipNames) {
    nodes[name] = {
      id: name,
      name,
      parent: null,
      children: [],
      isTip: true,
      trueAge: 0,
      age: 0,
      calibration: null,
      observedDistances: {},
    };
  }

  // Internal nodes with approximate ages (Ma)
  const internals: Array<{ id: string; name: string; age: number; children: string[] }> = [
    { id: 'N1', name: 'Human-Chimp', age: 7, children: ['Human', 'Chimp'] },
    { id: 'N2', name: 'African Apes', age: 10, children: ['N1', 'Gorilla'] },
    { id: 'N3', name: 'Great Apes', age: 14, children: ['N2', 'Orangutan'] },
    { id: 'N4', name: 'Hominoidea', age: 20, children: ['N3', 'Gibbon'] },
    { id: 'N5', name: 'Catarrhini', age: 30, children: ['N4', 'OWMonkey'] },
    { id: 'N6', name: 'Anthropoidea', age: 45, children: ['N5', 'NWMonkey'] },
    { id: 'N7', name: 'Primates', age: 65, children: ['N6', 'Lemur'] },
  ];

  for (const n of internals) {
    nodes[n.id] = {
      id: n.id,
      name: n.name,
      parent: null,
      children: n.children,
      isTip: false,
      trueAge: n.age,
      age: n.age,
      calibration: null,
      observedDistances: {},
    };
    for (const childId of n.children) {
      nodes[childId].parent = n.id;
    }
  }

  // Set default calibrations on some nodes
  nodes['N1'].calibration = { type: 'lognormal', offset: 5.0, mu: 0.5, sigma: 0.5 };
  nodes['N5'].calibration = { type: 'uniform', min: 25, max: 35 };
  nodes['N7'].calibration = { type: 'exponential', offset: 55, lambda: 0.1 };

  // Simulate observed distances
  const rate = 0.01; // substitutions per site per Myr
  for (const id of Object.keys(nodes)) {
    const node = nodes[id];
    if (node.parent) {
      const parentNode = nodes[node.parent];
      const childAge = node.isTip ? 0 : node.trueAge;
      const trueDistance = rate * (parentNode.trueAge - childAge);
      const obsDist = simulateObservedDistanceFast(trueDistance, 1000);
      node.observedDistances[node.parent] = obsDist;
    }
  }

  return nodes;
}
