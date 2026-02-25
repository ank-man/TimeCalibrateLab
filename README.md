# TimeCalibrateLab

**An Interactive Bayesian Phylogenetic Time Calibration Simulator**

TimeCalibrateLab is a browser-based educational tool that makes Bayesian phylogenetic time calibration intuitive through interactive probabilistic visualization. It demonstrates the core equation:

```
Posterior ∝ Likelihood × Prior
```

## Features

- **Interactive Phylogenetic Tree** — Click nodes, view calibrations, see posterior updates in real time
- **Calibration Control** — Configure Uniform, Lognormal, and Exponential fossil calibrations with live density plots
- **Clock Models** — Toggle between strict and relaxed (uncorrelated lognormal) molecular clocks
- **Data Panel** — Adjust alignment length and substitution rate; see how data strength affects inference
- **Prior vs Posterior Visualization** — D3.js density overlays with HPD interval shading
- **MCMC Engine** — Browser-based Metropolis-Hastings sampler running in a Web Worker
- **Educational Modes** — Beginner mode with tooltips; Advanced mode with trace plots and full parameter editing
- **Export** — Download posterior samples (CSV), tree (JSON), figures (PNG), session configs
- **Example Scenarios** — Pre-built scenarios: wide calibration, strong data, conflicting calibrations, high relaxed variance
- **Dark/Light Theme** — Professional academic styling

## Quick Start

```bash
cd timecalibratelab
npm install
npm run dev
```

Open http://localhost:5173 in your browser.

## Tech Stack

- **React 18** + **TypeScript** — UI framework
- **Vite 5** — Build tool
- **D3.js** — Tree visualization, density plots, trace plots
- **Zustand** — Global state management
- **Tailwind CSS 3** — Styling
- **Lucide React** — Icons
- **Web Workers** — Off-thread MCMC sampling

## Project Structure

```
src/
├── components/
│   ├── TreeView.tsx          # D3 phylogenetic tree with interactive nodes
│   ├── CalibrationPanel.tsx  # Fossil calibration configuration
│   ├── ClockPanel.tsx        # Strict/relaxed clock controls + rate heatmap
│   ├── DataPanel.tsx         # Alignment length, likelihood model controls
│   └── DensityPlot.tsx       # D3 prior vs posterior density overlay
├── math/
│   ├── priors.ts             # Uniform, Lognormal, Exponential, Gamma distributions
│   ├── likelihood.ts         # JC69 model + distance approximation
│   ├── mcmc.ts               # Metropolis-Hastings MCMC engine + example tree
│   └── hpd.ts                # HPD intervals, KDE, ESS computation
├── store/
│   └── useStore.ts           # Zustand global state (tree, calibrations, MCMC)
├── utils/
│   └── sampling.ts           # CSV/JSON export, file download utilities
├── workers/
│   └── mcmcWorker.ts         # Web Worker for off-thread MCMC
├── App.tsx                   # Main 4-panel layout with header/footer
├── main.tsx                  # Entry point
└── index.css                 # Tailwind directives + custom styles
```

## Documentation

The `docs/` directory contains comprehensive resources:

- **`YouTube_Script_Bayesian_Time_Calibration.md`** — Full ~60-minute video script explaining every concept from first principles
- **`Pipeline_OrthoFinder_to_MCMCTree.md`** — Step-by-step guide: proteomes → OrthoFinder → MAFFT → MCMCTree → dated tree
- **`MCMCTree_Parameters_Complete_Reference.md`** — Every MCMCTree control file parameter explained with math, tables, and guidance
- **`pipeline_orthofinder_to_mcmctree.sh`** — Automated bash pipeline script

## Mathematical Model

### Tree Parameterization
- Node ages `T_i`, branch lengths `d = r × t`
- Strict clock: single global rate `r`
- Relaxed clock: `r_i ~ LogNormal(μ, σ²)` per branch

### Calibration Priors
- **Uniform:** `p(T) = 1/(max - min)`
- **Lognormal with offset:** `T = offset + LogNormal(μ, σ)`
- **Exponential with offset:** `T = offset + Exp(λ)`

### Likelihood (JC69)
- P(same base) = `1/4 + 3/4 × e^(-4d/3)`
- P(different base) = `1/4 - 1/4 × e^(-4d/3)`
- Log-likelihood summed over L alignment sites

### MCMC
- Metropolis-Hastings with random walk proposals on node ages and rates
- 2,000–5,000 iterations; 25% burn-in
- Posterior mean, median, 95% HPD from kernel density estimation

## Example Scenarios

| Scenario | What It Demonstrates |
|----------|---------------------|
| Wide Calibration | Uniform(50,150) — posterior ≈ prior when alignment is short |
| Strong Data | L=10,000 — data overwhelms the prior |
| Conflicting Calibrations | Two nodes with incompatible priors — conflict warning |
| High Relaxed Variance | σ=1.0 — posterior spreads dramatically |

## Learning Lab

Click the **Learning Lab** button in the header to access 6 structured modules with guided lessons:

| Module | Topic | Lessons |
|--------|-------|---------|
| 1. Foundations | Why time-calibrate? The Bayesian equation. | 2 lessons, interactive tree exploration |
| 2. Priors & Calibrations | Uniform, Lognormal, Exponential distributions. Conflicts. | 2 lessons with live parameter adjustment |
| 3. Molecular Clock | Strict vs relaxed, the σ effect, rate heatmaps. | 1 lesson with clock model comparison |
| 4. Power of Data | Alignment length effect, JC69 vs distance likelihood. | 2 lessons showing data-prior interplay |
| 5. MCMC Sampling | Metropolis-Hastings, acceptance rates, trace plots, HPD intervals. | 2 lessons with convergence diagnostics |
| 6. Putting It Together | Complete analysis workflow from calibrations to interpretation. | 1 capstone lesson |

Each lesson includes **knowledge check quizzes**, **hands-on actions** with the simulator, and **hints**.

## Deploy to GitHub Pages

### Automatic (GitHub Actions)

1. Push this repo to GitHub
2. Go to **Settings → Pages → Source → GitHub Actions**
3. Push to `main` — the `.github/workflows/deploy.yml` workflow builds and deploys automatically
4. Your site will be live at `https://<username>.github.io/timecalibratelab/`

**Important:** Edit `vite.config.ts` and change `base: '/timecalibratelab/'` to match your actual repo name.

### Manual

```bash
npm run build
# Upload the dist/ folder to any static host (Netlify, Vercel, GitHub Pages, etc.)
```

### Local Preview of Production Build

```bash
npm run build
npx vite preview
```

## License

MIT

