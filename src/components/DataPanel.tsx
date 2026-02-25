/**
 * Data Panel - Controls alignment length, substitution rate, likelihood mode.
 * Shows likelihood curve for selected node.
 */

import { useStore } from '../store/useStore';
import { useRef, useEffect } from 'react';
import * as d3 from 'd3';
import { likelihoodCurve } from '../math/likelihood';
import { Info } from 'lucide-react';
import { useState } from 'react';

const DATA_TOOLTIPS = {
  alignment: 'Alignment length is the number of nucleotide sites compared. More sites = more data = tighter likelihood. With 10,000 sites, the data strongly constrains divergence times.',
  jc69: 'The JC69 (Jukes-Cantor 1969) model assumes equal base frequencies and equal substitution rates. It is the simplest nucleotide substitution model.',
  distance: 'The distance approximation uses a Gaussian likelihood centered on the observed pairwise distance. It is faster but less accurate than the full JC69 model.',
  likelihood: 'The likelihood measures how probable the observed sequence data is given the model parameters (divergence times and rates). It is the data-driven component of the posterior.',
};

export default function DataPanel() {
  const alignmentLength = useStore((s) => s.alignmentLength);
  const substitutionRate = useStore((s) => s.substitutionRate);
  const likelihoodMode = useStore((s) => s.likelihoodMode);
  const setAlignmentLength = useStore((s) => s.setAlignmentLength);
  const setSubstitutionRate = useStore((s) => s.setSubstitutionRate);
  const setLikelihoodMode = useStore((s) => s.setLikelihoodMode);
  const selectedNodeId = useStore((s) => s.selectedNodeId);
  const nodes = useStore((s) => s.nodes);
  const globalRate = useStore((s) => s.globalRate);
  const mode = useStore((s) => s.mode);
  const theme = useStore((s) => s.theme);

  const [showTooltip, setShowTooltip] = useState<string | null>(null);

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-bold uppercase tracking-wider text-slate-400">
          Sequence Data
        </h3>
        {mode === 'beginner' && (
          <button
            onClick={() => setShowTooltip(showTooltip === 'likelihood' ? null : 'likelihood')}
            className="p-1 rounded hover:bg-slate-700 transition-colors"
            title="Info about likelihood"
          >
            <Info size={14} className="text-slate-400" />
          </button>
        )}
      </div>

      {showTooltip && DATA_TOOLTIPS[showTooltip as keyof typeof DATA_TOOLTIPS] && (
        <div className="p-2 rounded-lg bg-green-900/30 border border-green-700/50 text-xs text-green-200 animate-fade-in">
          {DATA_TOOLTIPS[showTooltip as keyof typeof DATA_TOOLTIPS]}
        </div>
      )}

      {/* Alignment length */}
      <div>
        <div className="flex justify-between items-center mb-1">
          <label className="text-xs text-slate-400">
            Alignment Length
            {mode === 'beginner' && (
              <button
                onClick={() => setShowTooltip(showTooltip === 'alignment' ? null : 'alignment')}
                className="ml-1"
                title="Info about alignment length"
              >
                <Info size={10} className="inline text-slate-500" />
              </button>
            )}
          </label>
          <span className="text-xs font-mono text-white">{alignmentLength.toLocaleString()} sites</span>
        </div>
        <input
          type="range"
          min={100}
          max={10000}
          step={100}
          value={alignmentLength}
          onChange={(e) => setAlignmentLength(parseInt(e.target.value))}
          className="w-full bg-slate-600 accent-green-500"
          title="Alignment length slider"
        />
        <div className="flex justify-between text-[10px] text-slate-500 mt-0.5">
          <span>100</span>
          <span>10,000</span>
        </div>
      </div>

      {/* Substitution rate */}
      <div>
        <div className="flex justify-between items-center mb-1">
          <label className="text-xs text-slate-400">Substitution Rate</label>
          <span className="text-xs font-mono text-white">{substitutionRate.toFixed(4)}</span>
        </div>
        <input
          type="range"
          min={0.001}
          max={0.05}
          step={0.001}
          value={substitutionRate}
          onChange={(e) => setSubstitutionRate(parseFloat(e.target.value))}
          className="w-full bg-slate-600 accent-green-500"
          title="Substitution rate slider"
        />
      </div>

      {/* Likelihood mode */}
      <div>
        <label className="text-xs text-slate-400 block mb-1">Likelihood Model</label>
        <div className="grid grid-cols-2 gap-1">
          <button
            onClick={() => setLikelihoodMode('jc69')}
            className={`px-3 py-1.5 rounded text-xs font-medium transition-colors ${
              likelihoodMode === 'jc69'
                ? 'bg-green-600 text-white'
                : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
            }`}
            title="JC69 substitution model"
          >
            JC69
            {mode === 'beginner' && (
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  setShowTooltip(showTooltip === 'jc69' ? null : 'jc69');
                }}
                className="ml-1"
                title="Info about JC69 model"
              >
                <Info size={8} className="inline" />
              </button>
            )}
          </button>
          <button
            onClick={() => setLikelihoodMode('distance')}
            className={`px-3 py-1.5 rounded text-xs font-medium transition-colors ${
              likelihoodMode === 'distance'
                ? 'bg-green-600 text-white'
                : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
            }`}
            title="Distance approximation model"
          >
            Distance
            {mode === 'beginner' && (
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  setShowTooltip(showTooltip === 'distance' ? null : 'distance');
                }}
                className="ml-1"
                title="Info about distance model"
              >
                <Info size={8} className="inline" />
              </button>
            )}
          </button>
        </div>
      </div>

      {/* Data strength indicator */}
      <div className="p-2 rounded bg-slate-800/50 border border-slate-700">
        <div className="text-xs text-slate-400 mb-1">Data Strength</div>
        <div className="w-full bg-slate-700 rounded-full h-2">
          <div
            className="bg-gradient-to-r from-green-600 to-green-400 h-2 rounded-full transition-all duration-300"
            style={{ width: `${Math.min(100, (alignmentLength / 10000) * 100)}%` }}
          />
        </div>
        <div className="text-[10px] text-slate-500 mt-1">
          {alignmentLength < 500 ? 'Weak — prior dominates' :
           alignmentLength < 2000 ? 'Moderate — balanced inference' :
           alignmentLength < 5000 ? 'Strong — data constrains' :
           'Very strong — data overwhelms prior'}
        </div>
      </div>

      {/* Likelihood curve for selected node */}
      {selectedNodeId && nodes[selectedNodeId] && !nodes[selectedNodeId].isTip && nodes[selectedNodeId].parent && (
        <div>
          <div className="text-xs text-slate-400 mb-1">Likelihood Curve</div>
          <div className="bg-slate-800/50 rounded-lg border border-slate-700 p-1">
            <LikelihoodCurveViz
              nodeId={selectedNodeId}
              parentId={nodes[selectedNodeId].parent!}
              rate={globalRate}
              alignmentLength={alignmentLength}
              likelihoodMode={likelihoodMode}
              theme={theme}
            />
          </div>
        </div>
      )}
    </div>
  );
}

function LikelihoodCurveViz({
  nodeId,
  parentId,
  rate,
  alignmentLength,
  likelihoodMode,
  theme,
}: {
  nodeId: string;
  parentId: string;
  rate: number;
  alignmentLength: number;
  likelihoodMode: 'jc69' | 'distance';
  theme: string;
}) {
  const svgRef = useRef<SVGSVGElement>(null);
  const nodes = useStore((s) => s.nodes);

  useEffect(() => {
    if (!svgRef.current) return;
    const node = nodes[nodeId];
    const parent = nodes[parentId];
    if (!node || !parent) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const width = 280;
    const height = 120;
    const margin = { top: 10, right: 10, bottom: 25, left: 35 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    const isDark = theme === 'dark';
    const textColor = isDark ? '#cbd5e1' : '#334155';
    const gridColor = isDark ? '#334155' : '#e2e8f0';

    const childAge = node.isTip ? 0 : node.age;
    const parentAge = parent.age;
    const obsDist = node.observedDistances[parentId] ?? rate * (parentAge - childAge);

    const data = likelihoodCurve(parentAge, childAge, rate, obsDist, alignmentLength, likelihoodMode);
    if (data.length === 0) return;

    const g = svg
      .attr('width', width)
      .attr('height', height)
      .append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    const xScale = d3.scaleLinear()
      .domain([d3.min(data, (d) => d.x) ?? 0, d3.max(data, (d) => d.x) ?? 1])
      .range([0, w]);

    const yScale = d3.scaleLinear()
      .domain([0, d3.max(data, (d) => d.y) ?? 1])
      .range([h, 0]);

    g.append('g')
      .attr('transform', `translate(0,${h})`)
      .call(d3.axisBottom(xScale).ticks(4).tickSize(-h))
      .call((g) => g.select('.domain').remove())
      .call((g) => g.selectAll('.tick line').attr('stroke', gridColor).attr('stroke-dasharray', '2,2'))
      .call((g) => g.selectAll('.tick text').attr('fill', textColor).attr('font-size', 9));

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(3).tickSize(-w))
      .call((g) => g.select('.domain').remove())
      .call((g) => g.selectAll('.tick line').attr('stroke', gridColor).attr('stroke-dasharray', '2,2'))
      .call((g) => g.selectAll('.tick text').attr('fill', textColor).attr('font-size', 9));

    const area = d3.area<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y0(h)
      .y1((d) => yScale(d.y))
      .curve(d3.curveBasis);

    const line = d3.line<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y((d) => yScale(d.y))
      .curve(d3.curveBasis);

    g.append('path')
      .datum(data)
      .attr('d', area)
      .attr('fill', '#22c55e')
      .attr('opacity', 0.15);

    g.append('path')
      .datum(data)
      .attr('d', line)
      .attr('fill', 'none')
      .attr('stroke', '#22c55e')
      .attr('stroke-width', 2);

    g.append('text')
      .attr('x', w / 2)
      .attr('y', h + 20)
      .attr('text-anchor', 'middle')
      .attr('fill', textColor)
      .attr('font-size', 9)
      .text('Node Age (Ma)');

  }, [nodeId, parentId, rate, alignmentLength, likelihoodMode, theme, nodes]);

  return <svg ref={svgRef} className="w-full" />;
}
