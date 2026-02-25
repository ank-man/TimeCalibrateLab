# Contributing to TimeCalibrateLab

Thank you for your interest in contributing! This project aims to make Bayesian phylogenetic time calibration accessible and intuitive.

## Development Setup

```bash
git clone <repo-url>
cd timecalibratelab
npm install
npm run dev
```

## Architecture Overview

- **`src/math/`** — Pure TypeScript math functions (no React dependencies). Unit-testable.
- **`src/store/`** — Zustand global state. All parameter changes invalidate MCMC samples and trigger re-computation.
- **`src/components/`** — React components. D3 rendering happens inside `useEffect` hooks on SVG refs.
- **`src/workers/`** — Web Worker for MCMC sampling. Communicates via `postMessage`.
- **`docs/`** — YouTube script, pipeline guide, MCMCTree parameter reference.

## Guidelines

### Code Style
- TypeScript strict mode
- Functional React components with hooks
- Tailwind CSS for styling (no CSS modules)
- D3 for all data visualization (density plots, trees, traces)

### Math Module Rules
- All distribution functions must have both PDF and log-PDF variants
- Log-space arithmetic for numerical stability
- Every sampling function must handle edge cases (negative values, zero variance)

### Adding a New Distribution
1. Add PDF, logPDF, and sample functions to `src/math/priors.ts`
2. Add the type to `CalibrationDistribution` union
3. Update `calibrationLogPDF`, `calibrationPDF`, `calibrationSample`, `calibrationDensityCurve`
4. Update the CalibrationPanel UI to expose the new distribution's parameters

### Adding a New Clock Model
1. Implement the rate model in `src/math/mcmc.ts`
2. Update `computeRateLogPrior` and `computeLogLikelihood`
3. Add the option to `ClockPanel.tsx`
4. Update the Zustand store

### Performance
- MCMC must run in a Web Worker (never block the main thread)
- D3 rendering should not exceed 100ms for the default 8-taxon tree
- Use `requestAnimationFrame` for animations if needed

## Reporting Issues

Please include:
- Browser and OS
- Steps to reproduce
- Console errors (if any)
- Screenshot of unexpected behavior

## Documentation Contributions

The `docs/` directory is as important as the code. Improvements to the YouTube script, pipeline guide, or parameter reference are very welcome. Keep explanations precise and jargon-free where possible.
