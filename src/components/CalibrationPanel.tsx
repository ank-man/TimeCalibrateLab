/**
 * Calibration Control Panel - Right sidebar.
 * Allows users to select nodes and configure calibration distributions.
 */

import { useState } from 'react';
import { useStore } from '../store/useStore';
import type { CalibrationDistribution } from '../math/priors';
import DensityPlot from './DensityPlot';
import { Info, AlertTriangle, Trash2 } from 'lucide-react';

interface TooltipInfo {
  title: string;
  text: string;
}

const TOOLTIPS: Record<string, TooltipInfo> = {
  prior: {
    title: 'What is a Prior?',
    text: 'A prior distribution represents our belief about a parameter (like node age) before seeing sequence data. It is often informed by fossil evidence.',
  },
  uniform: {
    title: 'Uniform Distribution',
    text: 'All values between min and max are equally likely. Use when you only know the range of possible ages.',
  },
  lognormal: {
    title: 'Lognormal Distribution',
    text: 'A right-skewed distribution useful when the minimum age is well-constrained (e.g., by a fossil) but the maximum is less certain.',
  },
  exponential: {
    title: 'Exponential Distribution',
    text: 'Places most probability mass near the offset. Useful when you expect the true age to be close to the minimum fossil age.',
  },
  hpd: {
    title: '95% HPD Interval',
    text: 'The Highest Posterior Density interval is the shortest interval containing 95% of the posterior probability. It represents our credible range for the node age.',
  },
};

export default function CalibrationPanel() {
  const nodes = useStore((s) => s.nodes);
  const selectedNodeId = useStore((s) => s.selectedNodeId);
  const selectNode = useStore((s) => s.selectNode);
  const setCalibration = useStore((s) => s.setCalibration);
  const posteriorSummaries = useStore((s) => s.posteriorSummaries);
  const conflicts = useStore((s) => s.conflicts);
  const mode = useStore((s) => s.mode);

  const [activeTooltip, setActiveTooltip] = useState<string | null>(null);

  const internalNodes = Object.values(nodes).filter((n) => !n.isTip);
  const selectedNode = selectedNodeId ? nodes[selectedNodeId] : null;
  const summary = selectedNodeId ? posteriorSummaries[selectedNodeId] : null;

  function handleCalTypeChange(type: 'uniform' | 'lognormal' | 'exponential') {
    if (!selectedNodeId) return;
    let cal: CalibrationDistribution;
    const currentAge = selectedNode?.age ?? 50;
    switch (type) {
      case 'uniform':
        cal = { type: 'uniform', min: currentAge * 0.5, max: currentAge * 1.5 };
        break;
      case 'lognormal':
        cal = { type: 'lognormal', offset: currentAge * 0.7, mu: 1.0, sigma: 0.5 };
        break;
      case 'exponential':
        cal = { type: 'exponential', offset: currentAge * 0.8, lambda: 0.1 };
        break;
    }
    setCalibration(selectedNodeId, cal);
  }

  function updateCalParam(key: string, value: number) {
    if (!selectedNodeId || !selectedNode?.calibration) return;
    const cal = { ...selectedNode.calibration, [key]: value } as CalibrationDistribution;
    setCalibration(selectedNodeId, cal);
  }

  function removeCalibration() {
    if (!selectedNodeId) return;
    setCalibration(selectedNodeId, null);
  }

  return (
    <div className="h-full flex flex-col overflow-y-auto">
      <div className="flex items-center justify-between mb-3">
        <h2 className="text-sm font-bold uppercase tracking-wider text-slate-400">
          Calibration
        </h2>
        {mode === 'beginner' && (
          <button
            onClick={() => setActiveTooltip(activeTooltip === 'prior' ? null : 'prior')}
            className="p-1 rounded hover:bg-slate-700 transition-colors"
            title="What is a prior?"
            aria-label="Prior info"
          >
            <Info size={14} className="text-slate-400" />
          </button>
        )}
      </div>

      {activeTooltip && TOOLTIPS[activeTooltip] && (
        <div className="mb-3 p-3 rounded-lg bg-blue-900/30 border border-blue-700/50 text-xs text-blue-200 animate-fade-in">
          <div className="font-semibold mb-1">{TOOLTIPS[activeTooltip].title}</div>
          {TOOLTIPS[activeTooltip].text}
        </div>
      )}

      {/* Conflict warnings */}
      {conflicts.length > 0 && (
        <div className="mb-3 p-2 rounded-lg bg-red-900/30 border border-red-700/50 text-xs text-red-300">
          <div className="flex items-center gap-1 font-semibold mb-1">
            <AlertTriangle size={12} /> Calibration Conflicts
          </div>
          {conflicts.map((c, i) => (
            <div key={i} className="ml-4">{c}</div>
          ))}
        </div>
      )}

      {/* Node selector */}
      <div className="mb-3">
        <label className="text-xs text-slate-400 block mb-1">Select Node</label>
        <select
          value={selectedNodeId ?? ''}
          onChange={(e) => selectNode(e.target.value || null)}
          className="w-full px-2 py-1.5 rounded bg-slate-700 border border-slate-600 text-sm text-white focus:outline-none focus:border-blue-500"
          title="Select a node to calibrate"
          aria-label="Node selector"
        >
          <option value="">-- Select a node --</option>
          {internalNodes.map((n) => (
            <option key={n.id} value={n.id}>
              {n.name} ({n.age.toFixed(1)} Ma)
              {n.calibration ? ' 🦴' : ''}
            </option>
          ))}
        </select>
      </div>

      {selectedNode && !selectedNode.isTip && (
        <div className="space-y-3 animate-fade-in">
          {/* Node info */}
          <div className="p-2 rounded bg-slate-800/50 border border-slate-700">
            <div className="text-xs text-slate-400">Current Age</div>
            <div className="text-lg font-mono font-bold text-white">
              {selectedNode.age.toFixed(1)} Ma
            </div>
            {summary && (
              <div className="mt-1 space-y-0.5">
                <div className="text-xs text-blue-400">
                  Posterior Mean: {summary.mean.toFixed(2)} Ma
                </div>
                <div className="text-xs text-blue-400">
                  95% HPD: [{summary.hpd95[0].toFixed(2)}, {summary.hpd95[1].toFixed(2)}]
                  {mode === 'beginner' && (
                    <button
                      onClick={() => setActiveTooltip(activeTooltip === 'hpd' ? null : 'hpd')}
                      className="ml-1 inline"
                    >
                      <Info size={10} className="inline text-slate-500" />
                    </button>
                  )}
                </div>
              </div>
            )}
          </div>

          {/* Calibration type */}
          <div>
            <label className="text-xs text-slate-400 block mb-1">Calibration Type</label>
            <div className="grid grid-cols-3 gap-1">
              {(['uniform', 'lognormal', 'exponential'] as const).map((type) => (
                <button
                  key={type}
                  onClick={() => handleCalTypeChange(type)}
                  className={`px-2 py-1.5 rounded text-xs font-medium transition-colors ${
                    selectedNode.calibration?.type === type
                      ? 'bg-amber-600 text-white'
                      : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
                  }`}
                >
                  {type.charAt(0).toUpperCase() + type.slice(1)}
                  {mode === 'beginner' && (
                    <button
                      onClick={(e) => {
                        e.stopPropagation();
                        setActiveTooltip(activeTooltip === type ? null : type);
                      }}
                      className="ml-0.5"
                    >
                      <Info size={8} className="inline text-slate-400" />
                    </button>
                  )}
                </button>
              ))}
            </div>
          </div>

          {/* Calibration parameters */}
          {selectedNode.calibration && (
            <div className="space-y-2">
              {selectedNode.calibration.type === 'uniform' && (
                <>
                  <ParamSlider
                    label="Min Age (Ma)"
                    value={selectedNode.calibration.min}
                    min={0}
                    max={selectedNode.age * 2}
                    step={0.5}
                    onChange={(v) => updateCalParam('min', v)}
                  />
                  <ParamSlider
                    label="Max Age (Ma)"
                    value={selectedNode.calibration.max}
                    min={0}
                    max={selectedNode.age * 3}
                    step={0.5}
                    onChange={(v) => updateCalParam('max', v)}
                  />
                </>
              )}
              {selectedNode.calibration.type === 'lognormal' && (
                <>
                  <ParamSlider
                    label="Offset (Ma)"
                    value={selectedNode.calibration.offset}
                    min={0}
                    max={selectedNode.age * 2}
                    step={0.5}
                    onChange={(v) => updateCalParam('offset', v)}
                  />
                  <ParamSlider
                    label="Mean (mu)"
                    value={selectedNode.calibration.mu}
                    min={0}
                    max={5}
                    step={0.1}
                    onChange={(v) => updateCalParam('mu', v)}
                  />
                  <ParamSlider
                    label="Std Dev (sigma)"
                    value={selectedNode.calibration.sigma}
                    min={0.1}
                    max={3}
                    step={0.1}
                    onChange={(v) => updateCalParam('sigma', v)}
                  />
                </>
              )}
              {selectedNode.calibration.type === 'exponential' && (
                <>
                  <ParamSlider
                    label="Offset (Ma)"
                    value={selectedNode.calibration.offset}
                    min={0}
                    max={selectedNode.age * 2}
                    step={0.5}
                    onChange={(v) => updateCalParam('offset', v)}
                  />
                  <ParamSlider
                    label="Lambda"
                    value={selectedNode.calibration.lambda}
                    min={0.01}
                    max={2}
                    step={0.01}
                    onChange={(v) => updateCalParam('lambda', v)}
                  />
                </>
              )}

              <button
                onClick={removeCalibration}
                className="flex items-center gap-1 px-2 py-1 text-xs text-red-400 hover:text-red-300 hover:bg-red-900/20 rounded transition-colors"
              >
                <Trash2 size={12} /> Remove Calibration
              </button>
            </div>
          )}

          {/* Density plot */}
          {selectedNodeId && (selectedNode.calibration || summary) && (
            <div className="mt-2">
              <div className="text-xs text-slate-400 mb-1">Prior vs Posterior</div>
              <div className="bg-slate-800/50 rounded-lg border border-slate-700 p-1">
                <DensityPlot nodeId={selectedNodeId} width={300} height={180} />
              </div>
            </div>
          )}
        </div>
      )}

      {!selectedNode && (
        <div className="text-xs text-slate-500 text-center mt-8">
          Click a node on the tree or select one above to configure calibration.
        </div>
      )}
    </div>
  );
}

function ParamSlider({
  label,
  value,
  min,
  max,
  step,
  onChange,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  step: number;
  onChange: (v: number) => void;
}) {
  return (
    <div>
      <div className="flex justify-between items-center mb-0.5">
        <label className="text-xs text-slate-400">{label}</label>
        <span className="text-xs font-mono text-white">{value.toFixed(2)}</span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="w-full bg-slate-600 accent-amber-500"
      />
    </div>
  );
}
