/**
 * Clock Model Panel - Controls for strict/relaxed molecular clock.
 * Shows rate parameters and branch rate heatmap for relaxed clock.
 */

import { useStore } from '../store/useStore';
import { Info } from 'lucide-react';
import { useState } from 'react';

const CLOCK_TOOLTIPS = {
  strict: 'A strict molecular clock assumes a single, constant substitution rate across all branches. All lineages evolve at the same speed.',
  relaxed: 'A relaxed molecular clock allows each branch to have its own rate, drawn from a lognormal distribution. The sigma parameter controls how much rates vary among branches.',
  rate: 'The substitution rate determines how fast sequences evolve. Typical mammalian rates are ~0.01 substitutions/site/Myr.',
  sigma: 'Sigma controls the amount of rate variation among branches. Higher sigma = more rate heterogeneity. sigma=0 is equivalent to a strict clock.',
};

export default function ClockPanel() {
  const clockType = useStore((s) => s.clockType);
  const globalRate = useStore((s) => s.globalRate);
  const rateAlpha = useStore((s) => s.rateAlpha);
  const rateBeta = useStore((s) => s.rateBeta);
  const relaxedSigma = useStore((s) => s.relaxedSigma);
  const setClockType = useStore((s) => s.setClockType);
  const setGlobalRate = useStore((s) => s.setGlobalRate);
  const setRateAlpha = useStore((s) => s.setRateAlpha);
  const setRateBeta = useStore((s) => s.setRateBeta);
  const setRelaxedSigma = useStore((s) => s.setRelaxedSigma);
  const mode = useStore((s) => s.mode);
  const samples = useStore((s) => s.samples);
  const nodes = useStore((s) => s.nodes);

  const [showTooltip, setShowTooltip] = useState<string | null>(null);

  // Compute branch rate summary for heatmap
  const branchRateSummary: Record<string, number> = {};
  if (samples && clockType === 'relaxed') {
    for (const [id, rateSamples] of Object.entries(samples.rates)) {
      if (rateSamples.length > 0) {
        branchRateSummary[id] = rateSamples.reduce((a, b) => a + b, 0) / rateSamples.length;
      }
    }
  }

  const rateValues = Object.values(branchRateSummary);
  const minRate = rateValues.length > 0 ? Math.min(...rateValues) : 0;
  const maxRate = rateValues.length > 0 ? Math.max(...rateValues) : 0.02;

  function rateToColor(rate: number): string {
    if (maxRate === minRate) return '#3b82f6';
    const t = (rate - minRate) / (maxRate - minRate);
    const r = Math.round(59 + t * (239 - 59));
    const g = Math.round(130 + t * (68 - 130));
    const b = Math.round(246 + t * (68 - 246));
    return `rgb(${r},${g},${b})`;
  }

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-bold uppercase tracking-wider text-slate-400">
          Clock Model
        </h3>
      </div>

      {/* Clock type toggle */}
      <div className="grid grid-cols-2 gap-1">
        <button
          onClick={() => setClockType('strict')}
          className={`px-3 py-2 rounded-lg text-xs font-medium transition-all ${
            clockType === 'strict'
              ? 'bg-blue-600 text-white shadow-lg shadow-blue-600/20'
              : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
          }`}
          title="Strict molecular clock"
        >
          Strict Clock
          {mode === 'beginner' && (
            <button
              onClick={(e) => {
                e.stopPropagation();
                setShowTooltip(showTooltip === 'strict' ? null : 'strict');
              }}
              className="ml-1"
              title="Info about strict clock"
            >
              <Info size={10} className="inline text-slate-300" />
            </button>
          )}
        </button>
        <button
          onClick={() => setClockType('relaxed')}
          className={`px-3 py-2 rounded-lg text-xs font-medium transition-all ${
            clockType === 'relaxed'
              ? 'bg-purple-600 text-white shadow-lg shadow-purple-600/20'
              : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
          }`}
          title="Relaxed molecular clock"
        >
          Relaxed Clock
          {mode === 'beginner' && (
            <button
              onClick={(e) => {
                e.stopPropagation();
                setShowTooltip(showTooltip === 'relaxed' ? null : 'relaxed');
              }}
              className="ml-1"
              title="Info about relaxed clock"
            >
              <Info size={10} className="inline text-slate-300" />
            </button>
          )}
        </button>
      </div>

      {showTooltip && CLOCK_TOOLTIPS[showTooltip as keyof typeof CLOCK_TOOLTIPS] && (
        <div className="p-2 rounded-lg bg-indigo-900/30 border border-indigo-700/50 text-xs text-indigo-200 animate-fade-in">
          {CLOCK_TOOLTIPS[showTooltip as keyof typeof CLOCK_TOOLTIPS]}
        </div>
      )}

      {/* Global rate */}
      <div>
        <div className="flex justify-between items-center mb-1">
          <label className="text-xs text-slate-400">
            Global Rate (subst/site/Myr)
            {mode === 'beginner' && (
              <button
                onClick={() => setShowTooltip(showTooltip === 'rate' ? null : 'rate')}
                className="ml-1"
                title="Info about substitution rate"
              >
                <Info size={10} className="inline text-slate-500" />
              </button>
            )}
          </label>
          <span className="text-xs font-mono text-white">{globalRate.toFixed(4)}</span>
        </div>
        <input
          type="range"
          min={0.001}
          max={0.05}
          step={0.001}
          value={globalRate}
          onChange={(e) => setGlobalRate(parseFloat(e.target.value))}
          className="w-full bg-slate-600 accent-blue-500"
          title="Global substitution rate"
        />
      </div>

      {/* Rate prior (Gamma) */}
      {mode === 'advanced' && (
        <div className="space-y-2 p-2 rounded bg-slate-800/50 border border-slate-700">
          <div className="text-xs text-slate-400 font-semibold">Rate Prior: Gamma(alpha, beta)</div>
          <div>
            <div className="flex justify-between">
              <label className="text-xs text-slate-400">Alpha</label>
              <span className="text-xs font-mono text-white">{rateAlpha.toFixed(1)}</span>
            </div>
            <input
              type="range"
              min={0.5}
              max={10}
              step={0.5}
              value={rateAlpha}
              onChange={(e) => setRateAlpha(parseFloat(e.target.value))}
              className="w-full bg-slate-600 accent-blue-500"
              title="Gamma alpha parameter"
            />
          </div>
          <div>
            <div className="flex justify-between">
              <label className="text-xs text-slate-400">Beta</label>
              <span className="text-xs font-mono text-white">{rateBeta.toFixed(0)}</span>
            </div>
            <input
              type="range"
              min={10}
              max={1000}
              step={10}
              value={rateBeta}
              onChange={(e) => setRateBeta(parseFloat(e.target.value))}
              className="w-full bg-slate-600 accent-blue-500"
              title="Gamma beta parameter"
            />
          </div>
        </div>
      )}

      {/* Relaxed clock sigma */}
      {clockType === 'relaxed' && (
        <div>
          <div className="flex justify-between items-center mb-1">
            <label className="text-xs text-slate-400">
              Rate Variance (sigma)
              {mode === 'beginner' && (
                <button
                  onClick={() => setShowTooltip(showTooltip === 'sigma' ? null : 'sigma')}
                  className="ml-1"
                  title="Info about sigma parameter"
                >
                  <Info size={10} className="inline text-slate-500" />
                </button>
              )}
            </label>
            <span className="text-xs font-mono text-white">{relaxedSigma.toFixed(2)}</span>
          </div>
          <input
            type="range"
            min={0.05}
            max={2.0}
            step={0.05}
            value={relaxedSigma}
            onChange={(e) => setRelaxedSigma(parseFloat(e.target.value))}
            className="w-full bg-slate-600 accent-purple-500"
            title="Relaxed clock sigma"
          />
        </div>
      )}

      {/* Branch rate heatmap */}
      {clockType === 'relaxed' && Object.keys(branchRateSummary).length > 0 && (
        <div className="space-y-1">
          <div className="text-xs text-slate-400 font-semibold">Branch Rate Heatmap</div>
          <div className="grid gap-0.5">
            {Object.entries(branchRateSummary).map(([id, rate]) => {
              const node = nodes[id];
              if (!node) return null;
              return (
                <div key={id} className="flex items-center gap-2 text-xs">
                  <div
                    className="w-4 h-3 rounded-sm"
                    style={{ backgroundColor: rateToColor(rate) }}
                  />
                  <span className="text-slate-400 truncate flex-1">
                    {node.name}
                  </span>
                  <span className="font-mono text-white">
                    {rate.toFixed(4)}
                  </span>
                </div>
              );
            })}
          </div>
          <div className="flex items-center gap-1 mt-1">
            <div className="text-[10px] text-slate-500">{minRate.toFixed(4)}</div>
            <div className="flex-1 h-2 rounded" style={{
              background: `linear-gradient(to right, rgb(59,130,246), rgb(239,68,68))`
            }} />
            <div className="text-[10px] text-slate-500">{maxRate.toFixed(4)}</div>
          </div>
        </div>
      )}
    </div>
  );
}
