/**
 * Global state management with Zustand for TimeCalibrateLab.
 */

import { create } from 'zustand';
import type { CalibrationDistribution } from '../math/priors';
import type { TreeNode, MCMCConfig, MCMCSamples } from '../math/mcmc';
import { createExampleTree } from '../math/mcmc';
import { computeHPD, posteriorMean, posteriorMedian, kernelDensityEstimate } from '../math/hpd';

export interface NodeSummary {
  mean: number;
  median: number;
  hpd95: [number, number];
  density: { x: number; y: number }[];
}

export type AppMode = 'beginner' | 'advanced';
export type Theme = 'light' | 'dark';
export type BranchLengthUnit = 'substitutions' | 'time';
export type ExampleScenario = 'wide_calibration' | 'strong_data' | 'conflicting_calibrations' | 'high_relaxed_variance';

export interface AppState {
  // Theme & mode
  theme: Theme;
  mode: AppMode;
  setTheme: (theme: Theme) => void;
  toggleTheme: () => void;
  setMode: (mode: AppMode) => void;

  // Tree
  nodes: Record<string, TreeNode>;
  selectedNodeId: string | null;
  branchLengthUnit: BranchLengthUnit;
  setNodes: (nodes: Record<string, TreeNode>) => void;
  selectNode: (id: string | null) => void;
  setBranchLengthUnit: (unit: BranchLengthUnit) => void;
  updateNodeAge: (id: string, age: number) => void;
  setCalibration: (nodeId: string, cal: CalibrationDistribution | null) => void;

  // Clock model
  clockType: 'strict' | 'relaxed';
  globalRate: number;
  rateAlpha: number;
  rateBeta: number;
  relaxedSigma: number;
  setClockType: (type: 'strict' | 'relaxed') => void;
  setGlobalRate: (rate: number) => void;
  setRateAlpha: (alpha: number) => void;
  setRateBeta: (beta: number) => void;
  setRelaxedSigma: (sigma: number) => void;

  // Data
  alignmentLength: number;
  substitutionRate: number;
  likelihoodMode: 'jc69' | 'distance';
  setAlignmentLength: (len: number) => void;
  setSubstitutionRate: (rate: number) => void;
  setLikelihoodMode: (mode: 'jc69' | 'distance') => void;

  // MCMC
  mcmcConfig: MCMCConfig;
  samples: MCMCSamples | null;
  posteriorSummaries: Record<string, NodeSummary>;
  isRunning: boolean;
  progress: number;
  acceptanceRate: number;
  setMCMCConfig: (config: Partial<MCMCConfig>) => void;
  setSamples: (samples: MCMCSamples, acceptanceRate: number) => void;
  setRunning: (running: boolean) => void;
  setProgress: (progress: number) => void;
  invalidateSamples: () => void;
  runMCMC: () => void;

  // Calibration conflicts
  conflicts: string[];
  checkConflicts: () => void;

  // Scenarios
  loadScenario: (scenario: ExampleScenario) => void;

  // Initialize
  initialize: () => void;
}

function computeSummaries(samples: MCMCSamples): Record<string, NodeSummary> {
  const summaries: Record<string, NodeSummary> = {};
  for (const [id, ageSamples] of Object.entries(samples.nodeAges)) {
    if (ageSamples.length === 0) continue;
    summaries[id] = {
      mean: posteriorMean(ageSamples),
      median: posteriorMedian(ageSamples),
      hpd95: computeHPD(ageSamples, 0.95),
      density: kernelDensityEstimate(ageSamples),
    };
  }
  return summaries;
}

function detectConflicts(nodes: Record<string, TreeNode>): string[] {
  const conflicts: string[] = [];
  for (const node of Object.values(nodes)) {
    if (node.isTip || !node.calibration) continue;
    // Check if calibration conflicts with parent/child constraints
    for (const childId of node.children) {
      const child = nodes[childId];
      if (!child.isTip && child.calibration) {
        // Both have calibrations — check overlap
        const parentRange = getCalibrationRange(node.calibration);
        const childRange = getCalibrationRange(child.calibration);
        if (childRange[1] > parentRange[0] && childRange[0] < parentRange[1]) {
          // Potential conflict if child could be older than parent
          if (childRange[1] > parentRange[1] * 0.9) {
            conflicts.push(`${child.name} calibration may conflict with ${node.name}`);
          }
        }
      }
    }
  }
  return conflicts;
}

function getCalibrationRange(cal: CalibrationDistribution): [number, number] {
  switch (cal.type) {
    case 'uniform':
      return [cal.min, cal.max];
    case 'lognormal': {
      const mean = cal.offset + Math.exp(cal.mu + cal.sigma ** 2 / 2);
      const sd = Math.sqrt((Math.exp(cal.sigma ** 2) - 1) * Math.exp(2 * cal.mu + cal.sigma ** 2));
      return [cal.offset, mean + 3 * sd];
    }
    case 'exponential': {
      return [cal.offset, cal.offset + 5 / cal.lambda];
    }
  }
}

export const useStore = create<AppState>((set, get) => ({
  // Theme & mode
  theme: 'dark',
  mode: 'beginner',
  setTheme: (theme) => {
    set({ theme });
    if (theme === 'dark') {
      document.documentElement.classList.add('dark');
    } else {
      document.documentElement.classList.remove('dark');
    }
  },
  toggleTheme: () => {
    const current = get().theme;
    get().setTheme(current === 'dark' ? 'light' : 'dark');
  },
  setMode: (mode) => set({ mode }),

  // Tree
  nodes: {},
  selectedNodeId: null,
  branchLengthUnit: 'time',
  setNodes: (nodes) => set({ nodes }),
  selectNode: (id) => set({ selectedNodeId: id }),
  setBranchLengthUnit: (unit) => set({ branchLengthUnit: unit }),
  updateNodeAge: (id, age) => {
    const nodes = { ...get().nodes };
    if (nodes[id]) {
      nodes[id] = { ...nodes[id], age };
      set({ nodes });
      get().invalidateSamples();
    }
  },
  setCalibration: (nodeId, cal) => {
    const nodes = { ...get().nodes };
    if (nodes[nodeId]) {
      nodes[nodeId] = { ...nodes[nodeId], calibration: cal };
      set({ nodes });
      get().checkConflicts();
      get().invalidateSamples();
    }
  },

  // Clock model
  clockType: 'strict',
  globalRate: 0.01,
  rateAlpha: 2.0,
  rateBeta: 200.0,
  relaxedSigma: 0.3,
  setClockType: (clockType) => { set({ clockType }); get().invalidateSamples(); },
  setGlobalRate: (globalRate) => { set({ globalRate }); get().invalidateSamples(); },
  setRateAlpha: (rateAlpha) => { set({ rateAlpha }); get().invalidateSamples(); },
  setRateBeta: (rateBeta) => { set({ rateBeta }); get().invalidateSamples(); },
  setRelaxedSigma: (relaxedSigma) => { set({ relaxedSigma }); get().invalidateSamples(); },

  // Data
  alignmentLength: 1000,
  substitutionRate: 0.01,
  likelihoodMode: 'jc69',
  setAlignmentLength: (alignmentLength) => { set({ alignmentLength }); get().invalidateSamples(); },
  setSubstitutionRate: (substitutionRate) => { set({ substitutionRate }); get().invalidateSamples(); },
  setLikelihoodMode: (likelihoodMode) => { set({ likelihoodMode }); get().invalidateSamples(); },

  // MCMC
  mcmcConfig: {
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
  },
  samples: null,
  posteriorSummaries: {},
  isRunning: false,
  progress: 0,
  acceptanceRate: 0,
  setMCMCConfig: (config) => set((s) => ({ mcmcConfig: { ...s.mcmcConfig, ...config } })),
  setSamples: (samples, acceptanceRate) => {
    const summaries = computeSummaries(samples);
    set({ samples, posteriorSummaries: summaries, acceptanceRate, isRunning: false, progress: 1 });
  },
  setRunning: (isRunning) => set({ isRunning }),
  setProgress: (progress) => set({ progress }),
  invalidateSamples: () => set({ samples: null, posteriorSummaries: {}, progress: 0 }),

  runMCMC: () => {
    const state = get();
    if (state.isRunning) return;

    set({ isRunning: true, progress: 0 });

    const config: MCMCConfig = {
      iterations: state.mcmcConfig.iterations,
      burnin: state.mcmcConfig.burnin,
      proposalSD: state.mcmcConfig.proposalSD,
      rateProposalSD: state.mcmcConfig.rateProposalSD,
      clockType: state.clockType,
      rateAlpha: state.rateAlpha,
      rateBeta: state.rateBeta,
      relaxedSigma: state.relaxedSigma,
      globalRate: state.globalRate,
      alignmentLength: state.alignmentLength,
      likelihoodMode: state.likelihoodMode,
    };

    // Use Web Worker if available
    try {
      const worker = new Worker(
        new URL('../workers/mcmcWorker.ts', import.meta.url),
        { type: 'module' }
      );

      worker.onmessage = (e) => {
        const { type, data } = e.data;
        if (type === 'progress') {
          set({ progress: data });
        } else if (type === 'result') {
          get().setSamples(data.samples, data.acceptanceRate);
          worker.terminate();
        } else if (type === 'error') {
          console.error('MCMC Worker error:', data);
          set({ isRunning: false });
          worker.terminate();
        }
      };

      worker.onerror = (err) => {
        console.error('Worker error:', err);
        // Fallback to main thread
        runOnMainThread(state.nodes, config, set, get);
      };

      worker.postMessage({ nodes: state.nodes, config });
    } catch {
      // Fallback to main thread if workers unavailable
      runOnMainThread(state.nodes, config, set, get);
    }
  },

  // Conflicts
  conflicts: [],
  checkConflicts: () => {
    const conflicts = detectConflicts(get().nodes);
    set({ conflicts });
  },

  // Scenarios
  loadScenario: (scenario) => {
    const nodes = createExampleTree();
    switch (scenario) {
      case 'wide_calibration':
        nodes['N7'].calibration = { type: 'uniform', min: 50, max: 150 };
        set({ nodes, alignmentLength: 200, clockType: 'strict' as const });
        break;
      case 'strong_data':
        set({ nodes, alignmentLength: 10000, clockType: 'strict' as const });
        break;
      case 'conflicting_calibrations':
        nodes['N4'].calibration = { type: 'uniform', min: 18, max: 22 };
        nodes['N3'].calibration = { type: 'uniform', min: 20, max: 25 };
        set({ nodes });
        break;
      case 'high_relaxed_variance':
        set({ nodes, clockType: 'relaxed' as const, relaxedSigma: 1.0 });
        break;
    }
    get().invalidateSamples();
    get().checkConflicts();
  },

  // Initialize
  initialize: () => {
    const nodes = createExampleTree();
    set({ nodes });
    document.documentElement.classList.add('dark');
    get().checkConflicts();
  },
}));

function runOnMainThread(
  nodes: Record<string, TreeNode>,
  config: MCMCConfig,
  set: (partial: Partial<AppState> | ((state: AppState) => Partial<AppState>)) => void,
  get: () => AppState
) {
  import('../math/mcmc').then(({ runMCMC }) => {
    const result = runMCMC(nodes, config, (pct) => {
      set({ progress: pct });
    });
    get().setSamples(result.samples, result.acceptanceRate);
  });
}
