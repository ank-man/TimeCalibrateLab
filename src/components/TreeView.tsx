/**
 * D3-based phylogenetic tree visualization with interactive nodes.
 * Supports time/substitution branch length display, fossil icons,
 * and animated posterior updates.
 */

import { useRef, useEffect, useCallback } from 'react';
import * as d3 from 'd3';
import { useStore } from '../store/useStore';
import type { TreeNode } from '../math/mcmc';

interface LayoutNode {
  id: string;
  x: number;
  y: number;
  node: TreeNode;
}

interface LayoutLink {
  source: LayoutNode;
  target: LayoutNode;
  rate?: number;
}

export default function TreeView() {
  const svgRef = useRef<SVGSVGElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  const nodes = useStore((s) => s.nodes);
  const selectedNodeId = useStore((s) => s.selectedNodeId);
  const selectNode = useStore((s) => s.selectNode);
  const posteriorSummaries = useStore((s) => s.posteriorSummaries);
  const branchLengthUnit = useStore((s) => s.branchLengthUnit);
  const clockType = useStore((s) => s.clockType);
  const globalRate = useStore((s) => s.globalRate);
  const theme = useStore((s) => s.theme);
  const samples = useStore((s) => s.samples);

  const computeLayout = useCallback((): { layoutNodes: LayoutNode[]; layoutLinks: LayoutLink[] } => {
    const nodeList = Object.values(nodes);
    if (nodeList.length === 0) return { layoutNodes: [], layoutLinks: [] };

    const root = nodeList.find((n) => !n.parent);
    if (!root) return { layoutNodes: [], layoutLinks: [] };

    // Compute tip order (left to right)
    const tipOrder: string[] = [];
    function orderTips(nodeId: string) {
      const node = nodes[nodeId];
      if (node.isTip) {
        tipOrder.push(nodeId);
        return;
      }
      for (const childId of node.children) {
        orderTips(childId);
      }
    }
    orderTips(root.id);

    // Get displayed ages (posterior mean if available, else current age)
    function getDisplayAge(nodeId: string): number {
      const node = nodes[nodeId];
      if (node.isTip) return 0;
      const summary = posteriorSummaries[nodeId];
      return summary ? summary.mean : node.age;
    }

    // Find max age for scaling
    const maxAge = Math.max(...nodeList.filter((n) => !n.isTip).map((n) => getDisplayAge(n.id)), 1);

    const margin = { top: 30, right: 120, bottom: 30, left: 30 };
    const svgEl = svgRef.current;
    const totalWidth = svgEl ? svgEl.clientWidth : 600;
    const totalHeight = Math.max(400, tipOrder.length * 45);
    const w = totalWidth - margin.left - margin.right;
    const h = totalHeight - margin.top - margin.bottom;

    const xScale = d3.scaleLinear().domain([maxAge * 1.05, 0]).range([0, w]);
    const yScale = d3.scaleLinear().domain([0, tipOrder.length - 1]).range([margin.top, h]);

    // Assign y positions
    const yPositions: Record<string, number> = {};
    tipOrder.forEach((id, i) => {
      yPositions[id] = yScale(i);
    });

    function computeY(nodeId: string): number {
      if (yPositions[nodeId] !== undefined) return yPositions[nodeId];
      const node = nodes[nodeId];
      const childYs = node.children.map((cid) => computeY(cid));
      const avg = childYs.reduce((a, b) => a + b, 0) / childYs.length;
      yPositions[nodeId] = avg;
      return avg;
    }
    computeY(root.id);

    const layoutNodes: LayoutNode[] = [];
    const layoutLinks: LayoutLink[] = [];

    for (const node of nodeList) {
      const age = getDisplayAge(node.id);
      const lNode: LayoutNode = {
        id: node.id,
        x: xScale(age) + margin.left,
        y: yPositions[node.id],
        node,
      };
      layoutNodes.push(lNode);
    }

    const nodeMap = new Map(layoutNodes.map((n) => [n.id, n]));

    for (const lNode of layoutNodes) {
      if (lNode.node.parent) {
        const parent = nodeMap.get(lNode.node.parent);
        if (parent) {
          const branchRate = samples?.rates[lNode.id]
            ? samples.rates[lNode.id].reduce((a, b) => a + b, 0) / samples.rates[lNode.id].length
            : globalRate;
          layoutLinks.push({ source: parent, target: lNode, rate: branchRate });
        }
      }
    }

    return { layoutNodes, layoutLinks };
  }, [nodes, posteriorSummaries, globalRate, samples]);

  useEffect(() => {
    if (!svgRef.current) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const { layoutNodes, layoutLinks } = computeLayout();
    if (layoutNodes.length === 0) return;

    const isDark = theme === 'dark';
    const textColor = isDark ? '#cbd5e1' : '#334155';
    const linkColor = isDark ? '#64748b' : '#94a3b8';
    const selectedColor = '#3b82f6';
    const calibratedColor = '#f59e0b';
    const nodeColor = isDark ? '#1e293b' : '#f8fafc';
    const nodeBorder = isDark ? '#475569' : '#cbd5e1';

    const g = svg.append('g');

    // Draw links (right-angle phylogram style)
    for (const link of layoutLinks) {
      // Horizontal from parent to vertical connector
      g.append('path')
        .attr('d', `M ${link.source.x},${link.source.y} H ${link.target.x} V ${link.target.y}`)
        .attr('fill', 'none')
        .attr('stroke', linkColor)
        .attr('stroke-width', clockType === 'relaxed' && link.rate ? Math.max(1, Math.min(6, link.rate * 300)) : 2)
        .attr('opacity', 0.8);

      // Branch length label
      if (branchLengthUnit === 'time') {
        const parentAge = layoutNodes.find((n) => n.id === link.source.id)?.node;
        const childAge = layoutNodes.find((n) => n.id === link.target.id)?.node;
        if (parentAge && childAge) {
          const pAge = posteriorSummaries[parentAge.id]?.mean ?? parentAge.age;
          const cAge = childAge.isTip ? 0 : (posteriorSummaries[childAge.id]?.mean ?? childAge.age);
          const branchLen = (pAge - cAge).toFixed(1);
          const mx = (link.source.x + link.target.x) / 2;
          g.append('text')
            .attr('x', mx)
            .attr('y', link.target.y - 5)
            .attr('text-anchor', 'middle')
            .attr('fill', textColor)
            .attr('font-size', 8)
            .attr('opacity', 0.6)
            .text(`${branchLen} Ma`);
        }
      }
    }

    // Draw nodes
    for (const lNode of layoutNodes) {
      const isSelected = lNode.id === selectedNodeId;
      const hasCalibration = lNode.node.calibration !== null;
      const hasPosterior = posteriorSummaries[lNode.id] !== undefined;

      const nodeGroup = g.append('g')
        .attr('transform', `translate(${lNode.x},${lNode.y})`)
        .attr('cursor', 'pointer')
        .on('click', () => {
          selectNode(lNode.id);
        })
        .on('mouseenter', (event) => {
          if (tooltipRef.current) {
            const summary = posteriorSummaries[lNode.id];
            let html = `<div class="font-semibold">${lNode.node.name}</div>`;
            if (!lNode.node.isTip) {
              html += `<div class="text-xs mt-1">Age: ${lNode.node.age.toFixed(1)} Ma</div>`;
            }
            if (lNode.node.calibration) {
              html += `<div class="text-xs text-amber-400">Calibrated (${lNode.node.calibration.type})</div>`;
            }
            if (summary) {
              html += `<div class="text-xs text-blue-400 mt-1">
                Mean: ${summary.mean.toFixed(2)} Ma<br/>
                95% HPD: [${summary.hpd95[0].toFixed(2)}, ${summary.hpd95[1].toFixed(2)}]
              </div>`;
            }
            tooltipRef.current.innerHTML = html;
            tooltipRef.current.style.display = 'block';
            tooltipRef.current.style.left = `${event.pageX + 12}px`;
            tooltipRef.current.style.top = `${event.pageY - 10}px`;
          }
        })
        .on('mouseleave', () => {
          if (tooltipRef.current) {
            tooltipRef.current.style.display = 'none';
          }
        });

      if (lNode.node.isTip) {
        // Tip: small circle + name
        nodeGroup.append('circle')
          .attr('r', 4)
          .attr('fill', isSelected ? selectedColor : (isDark ? '#94a3b8' : '#64748b'))
          .attr('stroke', isSelected ? selectedColor : 'none')
          .attr('stroke-width', 2);

        nodeGroup.append('text')
          .attr('x', 10)
          .attr('y', 4)
          .attr('fill', textColor)
          .attr('font-size', 11)
          .attr('font-style', 'italic')
          .text(lNode.node.name);
      } else {
        // Internal node
        const radius = isSelected ? 10 : 7;

        // HPD bar
        if (hasPosterior) {
          const summary = posteriorSummaries[lNode.id];
          const nodeList = Object.values(nodes);
          const maxAge = Math.max(...nodeList.filter((n) => !n.isTip).map((n) => {
            const s = posteriorSummaries[n.id];
            return s ? s.mean : n.age;
          }), 1);
          const svgEl = svgRef.current!;
          const totalWidth = svgEl.clientWidth;
          const margin = { left: 30, right: 120 };
          const w = totalWidth - margin.left - margin.right;
          const xScaleLocal = d3.scaleLinear().domain([maxAge * 1.05, 0]).range([0, w]);

          const hpdLeft = xScaleLocal(summary.hpd95[1]) + margin.left - lNode.x;
          const hpdRight = xScaleLocal(summary.hpd95[0]) + margin.left - lNode.x;

          nodeGroup.append('rect')
            .attr('x', hpdLeft)
            .attr('y', -3)
            .attr('width', Math.max(0, hpdRight - hpdLeft))
            .attr('height', 6)
            .attr('rx', 3)
            .attr('fill', '#3b82f6')
            .attr('opacity', 0.3);
        }

        // Node circle
        nodeGroup.append('circle')
          .attr('r', radius)
          .attr('fill', isSelected ? selectedColor : nodeColor)
          .attr('stroke', isSelected ? selectedColor : (hasCalibration ? calibratedColor : nodeBorder))
          .attr('stroke-width', isSelected ? 3 : (hasCalibration ? 2.5 : 1.5));

        // Fossil icon for calibrated nodes
        if (hasCalibration) {
          nodeGroup.append('text')
            .attr('x', 0)
            .attr('y', 3.5)
            .attr('text-anchor', 'middle')
            .attr('font-size', radius > 8 ? 10 : 8)
            .text('🦴');
        }

        // Node name label above
        nodeGroup.append('text')
          .attr('x', 0)
          .attr('y', -radius - 5)
          .attr('text-anchor', 'middle')
          .attr('fill', textColor)
          .attr('font-size', 9)
          .attr('font-weight', isSelected ? 'bold' : 'normal')
          .text(lNode.node.name);
      }
    }

    // Time axis at bottom
    const nodeList = Object.values(nodes);
    const maxAge = Math.max(...nodeList.filter((n) => !n.isTip).map((n) => {
      const s = posteriorSummaries[n.id];
      return s ? s.mean : n.age;
    }), 1);
    const svgEl = svgRef.current!;
    const totalWidth = svgEl.clientWidth;
    const margin = { left: 30, right: 120, bottom: 30 };
    const w = totalWidth - margin.left - margin.right;
    const totalHeight = Math.max(400, Object.values(nodes).filter((n) => n.isTip).length * 45);
    const xScale = d3.scaleLinear().domain([maxAge * 1.05, 0]).range([0, w]);

    const axisG = svg.append('g')
      .attr('transform', `translate(${margin.left},${totalHeight - margin.bottom})`);

    axisG.call(d3.axisBottom(xScale).ticks(8))
      .call((g) => g.selectAll('text').attr('fill', textColor).attr('font-size', 10))
      .call((g) => g.selectAll('line').attr('stroke', textColor))
      .call((g) => g.select('.domain').attr('stroke', textColor));

    axisG.append('text')
      .attr('x', w / 2)
      .attr('y', 28)
      .attr('text-anchor', 'middle')
      .attr('fill', textColor)
      .attr('font-size', 11)
      .text('Time (Ma)');

  }, [nodes, selectedNodeId, posteriorSummaries, branchLengthUnit, clockType, globalRate, theme, computeLayout, samples, selectNode]);

  const tipCount = Object.values(nodes).filter((n) => n.isTip).length;
  const svgHeight = Math.max(400, tipCount * 45 + 30);

  return (
    <div className="relative w-full h-full overflow-auto">
      <svg
        ref={svgRef}
        width="100%"
        height={svgHeight}
        className="select-none"
      />
      <div
        ref={tooltipRef}
        className="d3-tooltip bg-slate-800 text-white border border-slate-600 shadow-xl"
        style={{ display: 'none' }}
      />
    </div>
  );
}
