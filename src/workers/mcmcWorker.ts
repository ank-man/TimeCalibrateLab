/**
 * Web Worker for running MCMC sampling off the main thread.
 * Prevents UI freezing during computation.
 */

import { runMCMC } from '../math/mcmc';
import type { TreeNode, MCMCConfig } from '../math/mcmc';

self.onmessage = (e: MessageEvent<{ nodes: Record<string, TreeNode>; config: MCMCConfig }>) => {
  try {
    const { nodes, config } = e.data;
    const result = runMCMC(nodes, config, (pct: number) => {
      self.postMessage({ type: 'progress', data: pct });
    });
    self.postMessage({ type: 'result', data: result });
  } catch (err) {
    self.postMessage({ type: 'error', data: String(err) });
  }
};
