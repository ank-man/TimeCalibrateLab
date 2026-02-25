/**
 * TimeCalibrateLab — Main Application Layout
 * Interactive Bayesian Phylogenetic Time Calibration Simulator
 */

import { useEffect, useState, useRef } from 'react';
import * as d3 from 'd3';
import { useStore } from './store/useStore';
import TreeView from './components/TreeView';
import CalibrationPanel from './components/CalibrationPanel';
import ClockPanel from './components/ClockPanel';
import DataPanel from './components/DataPanel';
import DensityPlot from './components/DensityPlot';
import LearningLab from './components/LearningLab';
import {
  Sun, Moon, Play, Download, FileJson, FileSpreadsheet,
  Image, BookOpen, Beaker, GraduationCap, RotateCcw,
  ChevronDown, Activity, Settings, FlaskConical,
} from 'lucide-react';
import { samplesToCSV, treeToJSON, downloadFile, downloadSVGasPNG } from './utils/sampling';

function App() {
  const initialize = useStore((s) => s.initialize);
  const theme = useStore((s) => s.theme);
  const toggleTheme = useStore((s) => s.toggleTheme);
  const mode = useStore((s) => s.mode);
  const setMode = useStore((s) => s.setMode);
  const runMCMC = useStore((s) => s.runMCMC);
  const isRunning = useStore((s) => s.isRunning);
  const progress = useStore((s) => s.progress);
  const acceptanceRate = useStore((s) => s.acceptanceRate);
  const samples = useStore((s) => s.samples);
  const nodes = useStore((s) => s.nodes);
  const posteriorSummaries = useStore((s) => s.posteriorSummaries);
  const loadScenario = useStore((s) => s.loadScenario);
  const selectedNodeId = useStore((s) => s.selectedNodeId);
  const invalidateSamples = useStore((s) => s.invalidateSamples);

  const [showAbout, setShowAbout] = useState(false);
  const [showScenarios, setShowScenarios] = useState(false);
  const [showLearningLab, setShowLearningLab] = useState(false);
  const [activeRightTab, setActiveRightTab] = useState<'calibration' | 'clock' | 'data'>('calibration');

  useEffect(() => {
    initialize();
  }, [initialize]);

  function handleExportJSON() {
    downloadFile(treeToJSON(nodes), 'timecalibratelab-tree.json', 'application/json');
  }

  function handleExportCSV() {
    if (!samples) return;
    const csv = samplesToCSV(samples.nodeAges, samples.globalRate, samples.rates);
    downloadFile(csv, 'posterior-samples.csv', 'text/csv');
  }

  function handleExportPNG() {
    const svg = document.querySelector('.tree-svg-container svg') as SVGSVGElement;
    if (svg) downloadSVGasPNG(svg, 'phylogeny.png');
  }

  function handleReset() {
    initialize();
    invalidateSamples();
  }

  const scenarioList = [
    { id: 'wide_calibration' as const, name: 'Wide Calibration', desc: 'Uniform(50,150) — posterior ≈ prior when data is weak' },
    { id: 'strong_data' as const, name: 'Strong Data', desc: 'L=10,000 — data overwhelms prior' },
    { id: 'conflicting_calibrations' as const, name: 'Conflicting Calibrations', desc: 'Two nodes with incompatible priors' },
    { id: 'high_relaxed_variance' as const, name: 'High Relaxed Variance', desc: 'σ=1.0 — wide posterior spread' },
  ];

  return (
    <div className={`min-h-screen flex flex-col ${theme === 'dark' ? 'bg-slate-900 text-slate-100' : 'bg-gray-50 text-gray-900'}`}>
      {/* Header */}
      <header className={`flex items-center justify-between px-4 py-2 border-b ${
        theme === 'dark' ? 'bg-slate-800/80 border-slate-700' : 'bg-white border-gray-200'
      } backdrop-blur-sm sticky top-0 z-50`}>
        <div className="flex items-center gap-3">
          <div className="flex items-center gap-2">
            <Beaker size={20} className="text-blue-500" />
            <h1 className="text-base font-bold tracking-tight">
              TimeCalibraté<span className="text-blue-500">Lab</span>
            </h1>
          </div>
          <span className={`text-[10px] px-1.5 py-0.5 rounded-full ${
            theme === 'dark' ? 'bg-slate-700 text-slate-400' : 'bg-gray-200 text-gray-500'
          }`}>v1.0</span>
        </div>

        <div className="flex items-center gap-2">
          {/* Mode toggle */}
          <div className={`flex rounded-lg overflow-hidden border ${
            theme === 'dark' ? 'border-slate-600' : 'border-gray-300'
          }`}>
            <button
              onClick={() => setMode('beginner')}
              className={`px-2.5 py-1 text-xs font-medium transition-colors flex items-center gap-1 ${
                mode === 'beginner'
                  ? 'bg-blue-600 text-white'
                  : theme === 'dark' ? 'bg-slate-700 text-slate-300 hover:bg-slate-600' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'
              }`}
              title="Beginner mode with guided tooltips"
            >
              <GraduationCap size={12} /> Beginner
            </button>
            <button
              onClick={() => setMode('advanced')}
              className={`px-2.5 py-1 text-xs font-medium transition-colors flex items-center gap-1 ${
                mode === 'advanced'
                  ? 'bg-purple-600 text-white'
                  : theme === 'dark' ? 'bg-slate-700 text-slate-300 hover:bg-slate-600' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'
              }`}
              title="Advanced mode with full parameter editing"
            >
              <Settings size={12} /> Advanced
            </button>
          </div>

          {/* Learning Lab */}
          <button
            onClick={() => setShowLearningLab(true)}
            className={`px-2.5 py-1 text-xs rounded-lg border flex items-center gap-1 transition-colors font-medium ${
              theme === 'dark' ? 'border-purple-600 bg-purple-900/30 text-purple-300 hover:bg-purple-900/50' : 'border-purple-300 bg-purple-50 text-purple-700 hover:bg-purple-100'
            }`}
            title="Open the Learning Lab — guided interactive lessons"
          >
            <FlaskConical size={12} /> Learning Lab
          </button>

          {/* Scenarios dropdown */}
          <div className="relative">
            <button
              onClick={() => setShowScenarios(!showScenarios)}
              className={`px-2.5 py-1 text-xs rounded-lg border flex items-center gap-1 transition-colors ${
                theme === 'dark' ? 'border-slate-600 bg-slate-700 text-slate-300 hover:bg-slate-600' : 'border-gray-300 bg-white text-gray-600 hover:bg-gray-100'
              }`}
              title="Load example scenario"
            >
              <BookOpen size={12} /> Scenarios <ChevronDown size={10} />
            </button>
            {showScenarios && (
              <div className={`absolute right-0 top-full mt-1 w-72 rounded-lg shadow-xl border z-50 ${
                theme === 'dark' ? 'bg-slate-800 border-slate-600' : 'bg-white border-gray-200'
              }`}>
                {scenarioList.map((s) => (
                  <button
                    key={s.id}
                    onClick={() => { loadScenario(s.id); setShowScenarios(false); }}
                    className={`w-full text-left px-3 py-2 text-xs transition-colors ${
                      theme === 'dark' ? 'hover:bg-slate-700' : 'hover:bg-gray-50'
                    }`}
                    title={s.desc}
                  >
                    <div className="font-semibold">{s.name}</div>
                    <div className={theme === 'dark' ? 'text-slate-400' : 'text-gray-500'}>{s.desc}</div>
                  </button>
                ))}
              </div>
            )}
          </div>

          {/* Reset */}
          <button
            onClick={handleReset}
            className={`p-1.5 rounded-lg transition-colors ${
              theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
            }`}
            title="Reset to defaults"
          >
            <RotateCcw size={14} />
          </button>

          {/* Theme */}
          <button
            onClick={toggleTheme}
            className={`p-1.5 rounded-lg transition-colors ${
              theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
            }`}
            title="Toggle theme"
          >
            {theme === 'dark' ? <Sun size={14} /> : <Moon size={14} />}
          </button>

          {/* About */}
          <button
            onClick={() => setShowAbout(!showAbout)}
            className={`px-2 py-1 text-xs rounded-lg transition-colors ${
              theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
            }`}
            title="About TimeCalibrateLab"
          >
            About
          </button>
        </div>
      </header>

      {/* About modal */}
      {showAbout && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/50 backdrop-blur-sm" onClick={() => setShowAbout(false)}>
          <div
            className={`max-w-lg w-full mx-4 rounded-xl p-6 shadow-2xl ${
              theme === 'dark' ? 'bg-slate-800 border border-slate-700' : 'bg-white border border-gray-200'
            }`}
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-lg font-bold mb-3">About TimeCalibrateLab</h2>
            <div className={`text-sm space-y-2 ${theme === 'dark' ? 'text-slate-300' : 'text-gray-600'}`}>
              <p>
                <strong>TimeCalibrateLab</strong> is an interactive educational tool for understanding
                Bayesian phylogenetic time calibration.
              </p>
              <p>It demonstrates the core equation:</p>
              <div className={`p-3 rounded-lg text-center font-mono text-base ${
                theme === 'dark' ? 'bg-slate-900' : 'bg-gray-100'
              }`}>
                Posterior ∝ Likelihood × Prior
              </div>
              <p>You can interactively explore how:</p>
              <ul className="list-disc list-inside space-y-1 ml-2">
                <li>Changing priors shifts the posterior</li>
                <li>Adding sequence data tightens the likelihood</li>
                <li>Relaxed clock changes node time variance</li>
                <li>Wide calibrations dominate when data is weak</li>
                <li>Strong data overwhelms priors</li>
              </ul>
              <p className="text-xs mt-3 opacity-60">
                This is a teaching tool, not a replacement for BEAST or MrBayes.
                Built with React, D3.js, and a browser-based MCMC engine.
              </p>
            </div>
            <button
              onClick={() => setShowAbout(false)}
              className="mt-4 px-4 py-2 bg-blue-600 text-white rounded-lg text-sm hover:bg-blue-700 transition-colors"
              title="Close about dialog"
            >
              Close
            </button>
          </div>
        </div>
      )}

      {/* Main content */}
      <div className="flex-1 flex overflow-hidden">
        {/* Left panel: Tree */}
        <div className={`flex-1 flex flex-col border-r ${
          theme === 'dark' ? 'border-slate-700' : 'border-gray-200'
        }`}>
          {/* Tree toolbar */}
          <div className={`flex items-center justify-between px-3 py-2 border-b ${
            theme === 'dark' ? 'border-slate-700 bg-slate-800/50' : 'border-gray-200 bg-gray-50'
          }`}>
            <div className="flex items-center gap-3">
              <h2 className="text-xs font-bold uppercase tracking-wider text-slate-400">Phylogeny</h2>

              {/* Bayesian equation reminder */}
              {mode === 'beginner' && (
                <div className={`text-[10px] px-2 py-0.5 rounded-full font-mono ${
                  theme === 'dark' ? 'bg-blue-900/30 text-blue-300 border border-blue-800/50' : 'bg-blue-50 text-blue-600 border border-blue-200'
                }`}>
                  Posterior ∝ Likelihood × Prior
                </div>
              )}
            </div>

            <div className="flex items-center gap-2">
              {/* Run MCMC button */}
              <button
                onClick={runMCMC}
                disabled={isRunning}
                className={`px-3 py-1.5 rounded-lg text-xs font-semibold flex items-center gap-1.5 transition-all ${
                  isRunning
                    ? 'bg-slate-600 text-slate-400 cursor-not-allowed'
                    : 'bg-blue-600 text-white hover:bg-blue-700 shadow-lg shadow-blue-600/20 hover:shadow-blue-600/30'
                }`}
                title="Run MCMC sampling"
              >
                {isRunning ? (
                  <>
                    <Activity size={12} className="animate-pulse-soft" />
                    Running... {Math.round(progress * 100)}%
                  </>
                ) : (
                  <>
                    <Play size={12} />
                    Run MCMC
                  </>
                )}
              </button>

              {/* Export buttons */}
              {samples && (
                <div className="flex items-center gap-1">
                  <button
                    onClick={handleExportJSON}
                    className={`p-1.5 rounded transition-colors ${
                      theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
                    }`}
                    title="Export tree as JSON"
                  >
                    <FileJson size={14} />
                  </button>
                  <button
                    onClick={handleExportCSV}
                    className={`p-1.5 rounded transition-colors ${
                      theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
                    }`}
                    title="Export posterior samples as CSV"
                  >
                    <FileSpreadsheet size={14} />
                  </button>
                  <button
                    onClick={handleExportPNG}
                    className={`p-1.5 rounded transition-colors ${
                      theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
                    }`}
                    title="Download tree as PNG"
                  >
                    <Image size={14} />
                  </button>
                  <button
                    onClick={() => {
                      const config = {
                        nodes,
                        clockType: useStore.getState().clockType,
                        globalRate: useStore.getState().globalRate,
                        alignmentLength: useStore.getState().alignmentLength,
                        likelihoodMode: useStore.getState().likelihoodMode,
                        relaxedSigma: useStore.getState().relaxedSigma,
                      };
                      downloadFile(JSON.stringify(config, null, 2), 'session-config.json', 'application/json');
                    }}
                    className={`p-1.5 rounded transition-colors ${
                      theme === 'dark' ? 'hover:bg-slate-700 text-slate-400' : 'hover:bg-gray-200 text-gray-500'
                    }`}
                    title="Save session configuration"
                  >
                    <Download size={14} />
                  </button>
                </div>
              )}
            </div>
          </div>

          {/* Progress bar */}
          {isRunning && (
            <div className="w-full h-1 bg-slate-700">
              <div
                className="h-1 bg-blue-500 transition-all duration-200"
                style={{ width: `${Math.round(progress * 100)}%` }}
              />
            </div>
          )}

          {/* MCMC summary bar */}
          {samples && !isRunning && (
            <div className={`flex items-center gap-4 px-3 py-1.5 text-[10px] border-b ${
              theme === 'dark' ? 'bg-slate-800/30 border-slate-700 text-slate-400' : 'bg-gray-50 border-gray-200 text-gray-500'
            }`}>
              <span>Acceptance: <strong className="text-green-400">{(acceptanceRate * 100).toFixed(1)}%</strong></span>
              <span>Post-burnin samples: <strong>{samples.globalRate.length}</strong></span>
              <span>Mean rate: <strong>{(samples.globalRate.reduce((a, b) => a + b, 0) / samples.globalRate.length).toFixed(5)}</strong></span>
            </div>
          )}

          {/* Tree view */}
          <div className="flex-1 overflow-auto tree-svg-container">
            <TreeView />
          </div>

          {/* Node posterior summaries */}
          {Object.keys(posteriorSummaries).length > 0 && (
            <div className={`border-t px-3 py-2 max-h-40 overflow-y-auto ${
              theme === 'dark' ? 'border-slate-700 bg-slate-800/30' : 'border-gray-200 bg-gray-50'
            }`}>
              <div className="text-[10px] uppercase tracking-wider text-slate-500 mb-1.5">Posterior Summaries</div>
              <div className="grid grid-cols-2 lg:grid-cols-3 gap-1.5">
                {Object.entries(posteriorSummaries).map(([id, summary]) => {
                  const node = nodes[id];
                  if (!node) return null;
                  return (
                    <div
                      key={id}
                      className={`px-2 py-1.5 rounded text-xs cursor-pointer transition-colors ${
                        selectedNodeId === id
                          ? 'bg-blue-600/20 border border-blue-500/30'
                          : theme === 'dark' ? 'bg-slate-800 hover:bg-slate-700' : 'bg-white hover:bg-gray-100 border border-gray-200'
                      }`}
                      onClick={() => useStore.getState().selectNode(id)}
                    >
                      <div className="font-semibold truncate">{node.name}</div>
                      <div className={theme === 'dark' ? 'text-slate-400' : 'text-gray-500'}>
                        {summary.mean.toFixed(2)} Ma [{summary.hpd95[0].toFixed(1)}, {summary.hpd95[1].toFixed(1)}]
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          )}

          {/* Selected node density plot at bottom of tree panel */}
          {selectedNodeId && posteriorSummaries[selectedNodeId] && (
            <div className={`border-t px-3 py-2 ${
              theme === 'dark' ? 'border-slate-700 bg-slate-800/30' : 'border-gray-200 bg-gray-50'
            }`}>
              <div className="text-[10px] uppercase tracking-wider text-slate-500 mb-1">
                {nodes[selectedNodeId]?.name} — Prior vs Posterior
              </div>
              <DensityPlot nodeId={selectedNodeId} width={500} height={160} />
            </div>
          )}
        </div>

        {/* Right sidebar: Tabbed panels */}
        <div className={`w-80 flex flex-col ${
          theme === 'dark' ? 'bg-slate-800/50' : 'bg-white'
        }`}>
          {/* Tabs */}
          <div className={`flex border-b ${
            theme === 'dark' ? 'border-slate-700' : 'border-gray-200'
          }`}>
            {[
              { key: 'calibration' as const, label: 'Calibration', icon: '🦴' },
              { key: 'clock' as const, label: 'Clock', icon: '⏱' },
              { key: 'data' as const, label: 'Data', icon: '🧬' },
            ].map((tab) => (
              <button
                key={tab.key}
                onClick={() => setActiveRightTab(tab.key)}
                className={`flex-1 px-2 py-2 text-xs font-medium transition-colors border-b-2 ${
                  activeRightTab === tab.key
                    ? 'border-blue-500 text-blue-500'
                    : `border-transparent ${theme === 'dark' ? 'text-slate-400 hover:text-slate-300' : 'text-gray-500 hover:text-gray-700'}`
                }`}
                title={`${tab.label} panel`}
              >
                <span className="mr-1">{tab.icon}</span>
                {tab.label}
              </button>
            ))}
          </div>

          {/* Panel content */}
          <div className="flex-1 overflow-y-auto p-3">
            {activeRightTab === 'calibration' && <CalibrationPanel />}
            {activeRightTab === 'clock' && <ClockPanel />}
            {activeRightTab === 'data' && <DataPanel />}
          </div>

          {/* Advanced: Trace plot section */}
          {mode === 'advanced' && samples && (
            <div className={`border-t p-3 ${
              theme === 'dark' ? 'border-slate-700' : 'border-gray-200'
            }`}>
              <div className="text-[10px] uppercase tracking-wider text-slate-500 mb-1">
                Log-Posterior Trace
              </div>
              <TracePlot data={samples.logPosterior} theme={theme} />
            </div>
          )}
        </div>
      </div>

      {/* Learning Lab overlay */}
      {showLearningLab && (
        <LearningLab onClose={() => setShowLearningLab(false)} />
      )}

      {/* Footer */}
      <footer className={`px-4 py-1.5 text-[10px] border-t flex justify-between ${
        theme === 'dark' ? 'bg-slate-800/50 border-slate-700 text-slate-500' : 'bg-gray-50 border-gray-200 text-gray-400'
      }`}>
        <span>TimeCalibrateLab — Interactive Bayesian Phylogenetic Time Calibration Simulator</span>
        <span>Educational tool | Not a replacement for BEAST</span>
      </footer>
    </div>
  );
}

/** Simple trace plot for log-posterior values */
function TracePlot({ data, theme }: { data: number[]; theme: string }) {
  const svgRef = useRef<SVGSVGElement>(null);

  useEffect(() => {
    if (!svgRef.current || data.length === 0) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const width = 260;
    const height = 60;
    const margin = { top: 5, right: 5, bottom: 15, left: 35 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    const isDark = theme === 'dark';
    const textColor = isDark ? '#64748b' : '#94a3b8';

    const g = svg.attr('width', width).attr('height', height)
      .append('g').attr('transform', `translate(${margin.left},${margin.top})`);

    const validData = data.filter((d) => isFinite(d));
    if (validData.length === 0) return;

    const xScale = d3.scaleLinear().domain([0, validData.length - 1]).range([0, w]);
    const yScale = d3.scaleLinear()
      .domain([d3.min(validData) ?? 0, d3.max(validData) ?? 0])
      .range([h, 0]);

    const line = d3.line<number>()
      .x((_, i) => xScale(i))
      .y((d) => yScale(d))
      .curve(d3.curveLinear);

    g.append('path')
      .datum(validData)
      .attr('d', line)
      .attr('fill', 'none')
      .attr('stroke', '#3b82f6')
      .attr('stroke-width', 0.8)
      .attr('opacity', 0.7);

    g.append('g')
      .attr('transform', `translate(0,${h})`)
      .call(d3.axisBottom(xScale).ticks(3).tickSize(0))
      .call((g) => g.select('.domain').attr('stroke', textColor))
      .call((g) => g.selectAll('text').attr('fill', textColor).attr('font-size', 8));

    g.append('g')
      .call(d3.axisLeft(yScale).ticks(2).tickSize(0).tickFormat(d3.format('.0f')))
      .call((g) => g.select('.domain').attr('stroke', textColor))
      .call((g) => g.selectAll('text').attr('fill', textColor).attr('font-size', 8));
  }, [data, theme]);

  return <svg ref={svgRef} className="w-full" />;
}

export default App;
