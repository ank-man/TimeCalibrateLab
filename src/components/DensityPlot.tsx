/**
 * D3-based density plot for prior vs posterior visualization.
 * Shows overlaid distributions with HPD shading.
 */

import { useRef, useEffect } from 'react';
import * as d3 from 'd3';
import { useStore } from '../store/useStore';
import { calibrationDensityCurve } from '../math/priors';

interface DensityPlotProps {
  nodeId: string;
  width?: number;
  height?: number;
  showPrior?: boolean;
  showPosterior?: boolean;
  showLikelihood?: boolean;
}

export default function DensityPlot({
  nodeId,
  width = 360,
  height = 200,
  showPrior = true,
  showPosterior = true,
}: DensityPlotProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const nodes = useStore((s) => s.nodes);
  const posteriorSummaries = useStore((s) => s.posteriorSummaries);
  const theme = useStore((s) => s.theme);

  const node = nodes[nodeId];
  const summary = posteriorSummaries[nodeId];

  useEffect(() => {
    if (!svgRef.current || !node) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const margin = { top: 15, right: 15, bottom: 30, left: 40 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    const g = svg
      .attr('width', width)
      .attr('height', height)
      .append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    const isDark = theme === 'dark';
    const textColor = isDark ? '#cbd5e1' : '#334155';
    const gridColor = isDark ? '#334155' : '#e2e8f0';

    // Collect all data
    let priorData: { x: number; y: number }[] = [];
    let posteriorData: { x: number; y: number }[] = [];

    if (showPrior && node.calibration) {
      priorData = calibrationDensityCurve(node.calibration);
    }

    if (showPosterior && summary) {
      posteriorData = summary.density;
    }

    if (priorData.length === 0 && posteriorData.length === 0) {
      g.append('text')
        .attr('x', w / 2)
        .attr('y', h / 2)
        .attr('text-anchor', 'middle')
        .attr('fill', textColor)
        .attr('font-size', 12)
        .text('No distribution data');
      return;
    }

    // Compute x domain
    const allX = [...priorData.map((d) => d.x), ...posteriorData.map((d) => d.x)];
    const xMin = d3.min(allX) ?? 0;
    const xMax = d3.max(allX) ?? 100;

    const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, w]);

    const allY = [...priorData.map((d) => d.y), ...posteriorData.map((d) => d.y)];
    const yMax = d3.max(allY) ?? 1;

    const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([h, 0]);

    // Grid lines
    g.append('g')
      .attr('transform', `translate(0,${h})`)
      .call(d3.axisBottom(xScale).ticks(5).tickSize(-h))
      .call((g) => g.select('.domain').remove())
      .call((g) => g.selectAll('.tick line').attr('stroke', gridColor).attr('stroke-dasharray', '2,2'))
      .call((g) => g.selectAll('.tick text').attr('fill', textColor).attr('font-size', 10));

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(3).tickSize(-w))
      .call((g) => g.select('.domain').remove())
      .call((g) => g.selectAll('.tick line').attr('stroke', gridColor).attr('stroke-dasharray', '2,2'))
      .call((g) => g.selectAll('.tick text').attr('fill', textColor).attr('font-size', 10));

    // X-axis label
    g.append('text')
      .attr('x', w / 2)
      .attr('y', h + 25)
      .attr('text-anchor', 'middle')
      .attr('fill', textColor)
      .attr('font-size', 10)
      .text('Age (Ma)');

    const line = d3
      .line<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y((d) => yScale(d.y))
      .curve(d3.curveBasis);

    const area = d3
      .area<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y0(h)
      .y1((d) => yScale(d.y))
      .curve(d3.curveBasis);

    // HPD shading
    if (showPosterior && summary) {
      const [hpdLow, hpdHigh] = summary.hpd95;
      g.append('rect')
        .attr('x', xScale(hpdLow))
        .attr('y', 0)
        .attr('width', Math.max(0, xScale(hpdHigh) - xScale(hpdLow)))
        .attr('height', h)
        .attr('fill', '#3b82f6')
        .attr('opacity', 0.1);
    }

    // Prior area + line
    if (priorData.length > 0) {
      g.append('path')
        .datum(priorData)
        .attr('d', area)
        .attr('fill', '#f59e0b')
        .attr('opacity', 0.15);

      g.append('path')
        .datum(priorData)
        .attr('d', line)
        .attr('fill', 'none')
        .attr('stroke', '#f59e0b')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '4,3');
    }

    // Posterior area + line
    if (posteriorData.length > 0) {
      g.append('path')
        .datum(posteriorData)
        .attr('d', area)
        .attr('fill', '#3b82f6')
        .attr('opacity', 0.2);

      g.append('path')
        .datum(posteriorData)
        .attr('d', line)
        .attr('fill', 'none')
        .attr('stroke', '#3b82f6')
        .attr('stroke-width', 2.5);
    }

    // Posterior mean line
    if (summary) {
      g.append('line')
        .attr('x1', xScale(summary.mean))
        .attr('x2', xScale(summary.mean))
        .attr('y1', 0)
        .attr('y2', h)
        .attr('stroke', '#3b82f6')
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '3,3');
    }

    // Legend
    const legendY = 5;
    if (priorData.length > 0) {
      g.append('line')
        .attr('x1', w - 80)
        .attr('x2', w - 65)
        .attr('y1', legendY)
        .attr('y2', legendY)
        .attr('stroke', '#f59e0b')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '4,3');
      g.append('text')
        .attr('x', w - 60)
        .attr('y', legendY + 4)
        .attr('fill', textColor)
        .attr('font-size', 9)
        .text('Prior');
    }
    if (posteriorData.length > 0) {
      g.append('line')
        .attr('x1', w - 80)
        .attr('x2', w - 65)
        .attr('y1', legendY + 15)
        .attr('y2', legendY + 15)
        .attr('stroke', '#3b82f6')
        .attr('stroke-width', 2.5);
      g.append('text')
        .attr('x', w - 60)
        .attr('y', legendY + 19)
        .attr('fill', textColor)
        .attr('font-size', 9)
        .text('Posterior');
    }
  }, [node, summary, width, height, showPrior, showPosterior, theme, nodeId]);

  if (!node) return null;

  return (
    <div className="animate-fade-in">
      <svg ref={svgRef} className="w-full" />
    </div>
  );
}
