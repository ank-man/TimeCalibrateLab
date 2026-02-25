# MCMCTree Parameters: Complete De-Black-Boxed Reference

Every parameter in MCMCTree's control file, explained with the underlying math, practical guidance, and common mistakes.

---

## Time Units Convention

**MCMCTree uses 100 Myr as its default time unit.**

| Real Age | MCMCTree Value |
|----------|---------------|
| 5 Ma     | 0.05          |
| 65 Ma    | 0.65          |
| 100 Ma   | 1.00          |
| 250 Ma   | 2.50          |

This means substitution rates are in units of substitutions/site/100Myr.
A mammalian nuclear rate of ~0.005 subst/site/Myr = 0.5 in MCMCTree units.

---

## Parameter-by-Parameter Reference

### `seed`
```
seed = -1
```

**What it does:** Initializes the pseudo-random number generator.

**Values:**
- `-1` — seed from system clock (different each run)
- Any positive integer — reproducible run

**Why it matters:** You MUST use different seeds for independent convergence chains. If two runs give the same output, you haven't tested convergence — you've just run the same chain twice.

**Recommendation:** Use `-1` for all runs; the clock-based seed ensures uniqueness.

---

### `seqfile`
```
seqfile = alignment.phy
```

**What it does:** Path to the sequence alignment in sequential PHYLIP format.

**Format requirements:**
```
  N_taxa  N_sites
Taxon1                         ATCGATCGATCG...
Taxon2                         ATCGATCGATCG...
```

- First line: number of taxa and sites, space-separated
- Taxon names: left-justified, max ~30 characters, padded with spaces
- Sequences: immediately follow the padded name
- Names must **exactly match** those in `treefile`
- For partitioned data (`ndata > 1`), stack multiple PHYLIP blocks

**Common mistake:** Species names differ between tree and alignment (e.g., `Homo_sapiens` vs `HomoSapiens`). MCMCTree will crash or produce garbage.

---

### `treefile`
```
treefile = calibrated_tree.nwk
```

**What it does:** Path to the rooted tree with fossil calibration annotations.

**Format:**
```
N_species  N_trees
((A,B)'B(0.05,0.10)',C)'B(0.20,0.40)';
```

- First line: counts
- Tree: Newick format, **no branch lengths**, with calibration strings after relevant `)`
- **Must be rooted** — MCMCTree does not root trees
- Calibrations are enclosed in single quotes

**Calibration distribution syntax:**

| Code | Distribution | Parameters | Example |
|------|-------------|------------|---------|
| `B(tL, tU)` | Uniform (soft bounds) | min, max | `'B(0.05, 0.10)'` — 5-10 Ma |
| `B(tL, tU, pL, pU)` | Uniform (custom tails) | min, max, left tail prob, right tail prob | `'B(0.05, 0.10, 0.01, 0.05)'` |
| `L(tL, p, c, pL)` | Truncated Cauchy (lower) | minimum, offset, scale, tail prob | `'L(0.05, 0.02, 0.1, 0.01)'` |
| `U(tU, p, c, pU)` | Truncated Cauchy (upper) | maximum, offset, scale, tail prob | `'U(1.00, 0.05, 0.1, 0.01)'` |
| `G(alpha, beta)` | Gamma | shape, rate | `'G(2, 10)'` — mean=0.2 |
| `S2N(loc, scale, shape)` | Skew-normal | location, scale, shape | `'S2N(0.65, 0.1, 5)'` |
| `ST(loc, scale, shape, df)` | Skew-t | location, scale, shape, degrees of freedom | `'ST(0.65, 0.1, 5, 4)'` |

**Understanding the Cauchy calibration `L(tL, p, c, pL)`:**
- `tL` = hard minimum (fossil age in 100 Myr)
- `p` = how far above tL the mode sits. Mode ≈ tL + p
- `c` = scale — larger c = wider distribution = more uncertainty
- `pL` = probability that the true age is BELOW tL. Usually 0.01 (1% — acknowledging fossil misdating)

**Example thought process:**
> "I have a fossil at 50 Ma for this clade. I think the clade is probably 55-60 Ma old but could be up to 80 Ma."
> → `'L(0.50, 0.08, 0.10, 0.01)'`
> - tL=0.50 (fossil at 50 Ma)
> - p=0.08 (mode at ~58 Ma)
> - c=0.10 (moderate uncertainty)
> - pL=0.01 (1% chance younger than 50 Ma)

---

### `ndata`
```
ndata = 1
```

**What it does:** Number of data partitions in the alignment file.

**If ndata = 1:** Single concatenated alignment.

**If ndata > 1:** The `seqfile` contains multiple PHYLIP blocks stacked vertically. Each partition gets its own rate multiplier, allowing rate heterogeneity among genes/partitions.

**Example for 3 partitions:**
```
ndata = 3
```
And `seqfile` contains:
```
  8  5000
[alignment block 1]
  8  3000
[alignment block 2]
  8  7000
[alignment block 3]
```

**When to partition:** If you have genes with very different evolutionary rates (e.g., mitochondrial + nuclear), partition them. For sets of similar nuclear genes, concatenation is usually fine.

---

### `seqtype`
```
seqtype = 0
```

**What it does:** Tells MCMCTree what kind of sequences you have.

| Value | Type | When to Use |
|-------|------|-------------|
| `0` | Nucleotides | Default for most analyses. Works with model=0-7. |
| `1` | Codons | More parameter-rich. Estimates dN/dS. Slower. Use for coding sequences when you want codon models. |
| `2` | Amino acids | When you only have protein alignments. Uses empirical amino acid models. |

**Recommendation:** `seqtype = 0` with nucleotide alignments is the standard approach.

---

### `usedata` ⭐ MOST IMPORTANT PARAMETER

```
usedata = 2
```

**What it does:** Controls how MCMCTree uses your sequence data.

| Value | Action | When to Use |
|-------|--------|-------------|
| `0` | **Ignore data entirely.** Samples from the prior only. | **ALWAYS RUN THIS FIRST** as a diagnostic. Compare the effective prior to your intended calibrations. |
| `1` | **Exact likelihood.** Computes the full phylogenetic likelihood at every MCMC step. | Small datasets only (<10 taxa, <1000 sites). Very slow for genomic data. |
| `2` | **Approximate likelihood.** Uses the Hessian matrix from `in.BV` to approximate the likelihood as a multivariate normal. | **Standard for all large analyses.** Must run `usedata = 3` first. |
| `3` | **Compute Hessian.** Calculates the gradient and Hessian of the log-likelihood at the ML branch lengths, writes `out.BV`. | **Run this once**, then rename `out.BV` → `in.BV`, then switch to `usedata = 2`. |

**The two-step workflow:**
```
Step 1: usedata = 3  →  produces out.BV  →  mv out.BV in.BV
Step 2: usedata = 2  →  reads in.BV  →  runs the actual MCMC
```

**The math behind the approximation:**

The exact log-likelihood L(b) as a function of branch lengths b is approximated by a second-order Taylor expansion around the ML estimates b̂:

```
L(b) ≈ L(b̂) + g'(b - b̂) + ½(b - b̂)'H(b - b̂)
```

Where:
- g = gradient (first derivatives) at b̂
- H = Hessian (second derivative matrix) at b̂

This is a multivariate normal approximation to the likelihood surface. It's extremely accurate when the likelihood is unimodal and well-behaved (true for most real datasets).

---

### `clock` ⭐ MAJOR SCIENTIFIC DECISION

```
clock = 2
```

**What it does:** Specifies the molecular clock model.

#### `clock = 1` — Strict Clock

**Model:** All branches share a single substitution rate r.

```
branch_length_i = r × time_i
```

**Prior on r:** Specified by `rgene_gamma`.

**When to use:** Almost never for real data. Useful as a baseline or for very clock-like data (e.g., some viral datasets).

**Effect on results:** Produces the narrowest credible intervals (because there's no rate uncertainty). But if the strict clock is wrong (it usually is), these narrow intervals are misleadingly precise.

#### `clock = 2` — Independent Rates (Uncorrelated Lognormal)

**Model:** Each branch i gets its own rate r_i drawn independently from a lognormal distribution:

```
log(r_i) ~ Normal(μ, σ²)
```

Where:
- μ = log of the overall mean rate (related to `rgene_gamma`)
- σ² = variance of log-rates among branches (estimated; prior from `sigma2_gamma`)

**When to use:** **Default choice for most analyses.** Makes no assumption about rate inheritance between parent and child branches.

**Effect on results:** Moderate credible intervals. Each branch rate is estimated independently, so distantly related lineages can have very different rates.

#### `clock = 3` — Autocorrelated Rates (Geometric Brownian Motion)

**Model:** A child branch's rate is correlated with its parent's rate:

```
log(r_child) ~ Normal(log(r_parent), σ² × t)
```

Where t is the time duration of the parent branch. Rates drift gradually along the tree like a random walk.

**When to use:** When you believe rates change gradually (e.g., correlated with body size, generation time). Often gives narrower intervals than clock=2 because the autocorrelation provides additional constraint.

**Effect on results:** Narrower credible intervals than clock=2. Rates change smoothly across the tree rather than jumping randomly.

**Recommendation:** Run BOTH clock=2 and clock=3 and compare. Report both sets of results. If they agree, great. If they disagree, discuss the biological plausibility of each model.

---

### `model`
```
model = 4
```

**What it does:** Nucleotide substitution model (only relevant for `seqtype = 0`).

| Value | Model | Parameters | Notes |
|-------|-------|-----------|-------|
| `0` | JC69 | None | Equal rates, equal base frequencies. Simplest. |
| `1` | K80 | κ (ti/tv ratio) | Distinguishes transitions from transversions. |
| `2` | F81 | π (base freqs) | Unequal base frequencies, equal rates. |
| `3` | F84 | κ, π | Transition bias + unequal frequencies. |
| `4` | **HKY85** | κ, π | **Most commonly used.** Transition/transversion bias + empirical base frequencies. |
| `5` | T92 | θ (GC content) | Simplified; rarely used. |
| `6` | TN93 | κ1, κ2, π | Separate AG and CT transition rates. |
| `7` | **GTR (REV)** | 6 rate params, π | Most general reversible model. Use for publication-quality. |

**Recommendation:** `model = 4` (HKY) is standard. Use `model = 7` (GTR) if you want maximum accuracy and have enough data. The choice usually has minor effects on divergence time estimates.

**For codons (`seqtype = 1`):** model specifies the codon model instead.
**For amino acids (`seqtype = 2`):** model specifies the amino acid replacement matrix.

---

### `alpha`
```
alpha = 0.5
```

**What it does:** Shape parameter of the discrete gamma distribution for among-site rate variation (ASRV).

**The model:** Not all sites evolve at the same rate. Some are conserved (slow), others are variable (fast). The gamma distribution models this heterogeneity:

```
Site rates ~ Gamma(α, α)    [mean = 1, variance = 1/α]
```

| α value | Rate variation | Meaning |
|---------|---------------|---------|
| `0` | None | All sites evolve at the same rate (rarely true) |
| `0.1` | Extreme | Very heterogeneous — some sites are nearly invariant, others evolve extremely fast |
| `0.5` | Substantial | Common default; good for most protein-coding genes |
| `1.0` | Moderate | Less heterogeneity |
| `5.0` | Mild | Sites are fairly homogeneous |
| `∞` | None | Equivalent to no gamma correction |

**How to set it:** Estimate it from a preliminary phylogenetic analysis:
```bash
# In IQ-TREE:
iqtree -s alignment.phy -m HKY+G4 -pre prelim
# Check the .iqtree output for "Gamma shape alpha"
```

Use that estimated value in MCMCTree.

---

### `ncatG`
```
ncatG = 5
```

**What it does:** Number of discrete categories to approximate the continuous gamma distribution for ASRV.

**The math:** The continuous gamma distribution is approximated by ncatG discrete rate categories, each with equal probability 1/ncatG. Category rates are the mean of each quantile interval.

| ncatG | Accuracy | Speed |
|-------|----------|-------|
| `1` | No rate variation (same as alpha=0) | Fastest |
| `4` | Good (standard) | Standard |
| `5` | Better (slight improvement over 4) | Slightly slower |
| `8` | Very good | Slower |

**Recommendation:** `ncatG = 4` or `5`. Going above 5 rarely changes results.

---

### `cleandata`
```
cleandata = 0
```

**What it does:** How to handle alignment gaps and ambiguous characters.

| Value | Action |
|-------|--------|
| `0` | Keep all sites; use proper ambiguity likelihood calculations |
| `1` | Delete any site with gaps or ambiguous characters in ANY taxon |

**Recommendation:** `cleandata = 0` — you've already trimmed your alignment with TrimAl; no need to discard more data. Setting `cleandata = 1` on gappy alignments can remove a large fraction of sites.

---

### `BDparas` — Birth-Death Tree Prior
```
BDparas = 1 1 0.1
```

**What it does:** Parameters for the birth-death process prior on divergence times: `λ μ ρ`

**The model:** Before looking at data or fossils, what do we believe about the shape of the species tree? The birth-death process generates tree topologies and node ages:

- **λ (birth/speciation rate):** Rate at which lineages split
- **μ (death/extinction rate):** Rate at which lineages go extinct
- **ρ (sampling fraction):** Proportion of extant species in your tree

**The math:**
```
P(tree | λ, μ, ρ) = [birth-death likelihood of node times and tree shape]
```

This prior says: "The pattern of divergence times should be consistent with a birth-death process with these parameters."

**How to set ρ:**
```
ρ = (number of species in your tree) / (total extant species in the clade)
```
Examples:
- 8 primates out of ~500 → ρ ≈ 0.016
- 20 mammals out of ~6,500 → ρ ≈ 0.003
- 50 angiosperms out of ~350,000 → ρ ≈ 0.00014

**How to set λ and μ:**
- Usually `λ ≈ μ` (moderate turnover)
- `λ = 1, μ = 1` is a common default
- The absolute values matter less than the ratio λ/μ
- Higher λ and μ → more node-rich tree prior
- This prior has **relatively little effect** when calibrations are informative

**Sensitivity check:** Run with different BDparas (e.g., `1 1 0.01` vs `2 2 0.1`) and verify results are stable. If they change substantially, your calibrations aren't constraining things enough.

---

### `rgene_gamma` ⭐ CRITICAL FOR RATE ESTIMATION
```
rgene_gamma = 2 20 1
```

**What it does:** Gamma-Dirichlet prior on overall substitution rates across partitions.

**Parameters:** `α β [concentration]`

The prior on the mean rate r across all branches:
```
r ~ Gamma(α, β)
Mean rate = α/β
Variance = α/β²
```

**Worked examples:**

| Group | Typical Rate (subst/site/Myr) | In 100Myr Units | Suggested rgene_gamma |
|-------|-------------------------------|-----------------|----------------------|
| Mammal nuclear | 0.002–0.005 | 0.2–0.5 | `2 5 1` (mean=0.4) |
| Mammal mitochondrial | 0.01–0.02 | 1.0–2.0 | `2 1.5 1` (mean=1.33) |
| Plant nuclear | 0.001–0.003 | 0.1–0.3 | `2 10 1` (mean=0.2) |
| Plant chloroplast | 0.0005–0.002 | 0.05–0.2 | `2 20 1` (mean=0.1) |
| Insect nuclear | 0.003–0.01 | 0.3–1.0 | `2 3 1` (mean=0.67) |
| Fungal nuclear | 0.001–0.005 | 0.1–0.5 | `2 7 1` (mean=0.29) |

**The α parameter controls how informative the prior is:**
- `α = 1` → exponential distribution (vague, long tail)
- `α = 2` → mildly informative (peaked but wide; **recommended default**)
- `α = 5` → moderately informative
- `α = 20` → very informative (use only if you're very confident about the rate)

**How to estimate the rate for your group:**
1. Run a quick tree search with branch lengths: `iqtree -s alignment.phy -m HKY+G`
2. Sum all branch lengths and divide by total tree time (from literature)
3. Or: use the root-to-tip distance for a known-age clade

**The third parameter (1 or 2):**
- `1` = alpha parameter of the Dirichlet for rate proportions across partitions (when ndata > 1)
- Higher values → rates are more similar across partitions
- For single-partition analyses, this doesn't matter

---

### `sigma2_gamma` ⭐ CONTROLS RATE VARIATION
```
sigma2_gamma = 1 10 1
```

**What it does:** Gamma prior on σ², the variance of log-rates among branches. Only relevant for `clock = 2` or `clock = 3`.

**Parameters:** `α β [concentration]`

```
σ² ~ Gamma(α, β)
Mean σ² = α/β
```

**What σ² means biologically:**

σ² controls how much substitution rates vary among branches. It's the single most important parameter for determining the width of your credible intervals.

| σ² value | Rate variation | Effect on CIs |
|----------|---------------|---------------|
| 0.01 | Nearly clock-like | Very narrow CIs |
| 0.05 | Low variation | Narrow CIs |
| 0.1 | Moderate variation | Moderate CIs |
| 0.5 | High variation | Wide CIs |
| 1.0 | Extreme variation | Very wide CIs |
| 2.0+ | Pathological | CIs so wide they're uninformative |

**Recommended settings:**

| Scenario | sigma2_gamma | Mean σ² |
|----------|-------------|---------|
| Clock-like data (e.g., slow-evolving nuclear) | `2 40 1` | 0.05 |
| Typical data (most analyses) | `1 10 1` | 0.10 |
| Variable-rate data (e.g., mixed nuclear+mito) | `1 4.5 1` | 0.22 |
| Highly variable (e.g., parasites, viruses) | `1 2 1` | 0.50 |

**The relationship to credible interval width:**

Roughly: doubling σ² increases the 95% HPD width by ~40-60%. This is why this parameter is so consequential.

**Prior sensitivity:** This is one of the most important parameters to test sensitivity for. Run with different sigma2_gamma values and report whether conclusions change.

---

### `RootAge`
```
RootAge = 'B(0.55, 1.00)'
```

**What it does:** Prior distribution on the age of the root node.

**Syntax:** Same calibration distribution syntax as tree node calibrations.

**Important:** This is required even if the root node already has a calibration in the tree file. MCMCTree uses this to set the overall time scale for the MCMC.

**Recommendation:** Use the same distribution as your root calibration in the tree file. Make sure it's consistent with your other calibrations (root must be older than all descendant calibrations).

---

### `finetune`
```
finetune = 1: .1 .1 .1 .1 .1 .1
```

**What it does:** Controls the step sizes for MCMC proposals.

**Format:** `auto: s1 s2 s3 s4 s5 s6`

- `auto` = `1` means auto-adjust during burnin (let MCMCTree tune itself). **Almost always use 1.**
- `auto` = `0` means use fixed step sizes (only for debugging)

**The six step sizes correspond to:**
1. Node ages (t)
2. Overall rates (μ) across partitions
3. Mixing parameter
4. Other parameters (κ, α, etc.)
5. Rate-related proposal
6. Another rate-related proposal

**When step sizes matter:**

If you see acceptance rates printed during the run:
- **20-40% acceptance** → good (optimal is ~23% for high-dimensional problems)
- **<10% acceptance** → steps too large, reduce the corresponding value
- **>70% acceptance** → steps too small, increase the value

With `finetune = 1:`, MCMCTree auto-tunes during burnin to achieve ~30% acceptance. This works well in practice.

---

### `burnin`
```
burnin = 100000
```

**What it does:** Number of MCMC iterations to discard before recording samples.

**Why it's needed:** The chain starts from arbitrary initial values. It needs time to "find" the region of high posterior probability. Samples from this initial wandering phase are not representative of the posterior.

**How much is enough?** Look at the trace plot in Tracer. The burnin should cover the entire initial transient phase where the log-posterior is climbing.

| Dataset size | Recommended burnin |
|-------------|-------------------|
| Small (<10 taxa, <1000 sites) | 10,000–50,000 |
| Medium (~20 taxa, ~10,000 sites) | 50,000–200,000 |
| Large (>20 taxa, genomic data) | 100,000–500,000 |

**Better too much than too little.** Discarding extra good samples is wasteful but harmless. Including pre-convergence samples biases your results.

---

### `sampfreq`
```
sampfreq = 50
```

**What it does:** Record every k-th MCMC sample (thinning).

**Why thin?** Consecutive MCMC samples are autocorrelated — each sample is similar to the previous one because proposals are small. Thinning reduces autocorrelation in the stored samples, improving ESS per sample.

| sampfreq | Disk usage | Autocorrelation | ESS per stored sample |
|----------|-----------|----------------|---------------------|
| 1 | Huge file | Very high | Low |
| 10 | Large | Moderate | Moderate |
| 50 | Manageable | Low | **Good** |
| 100 | Small | Very low | Good |
| 500 | Tiny | Negligible | Good, but wastes computation |

**Recommendation:** `sampfreq = 50` is standard. Use 100 for very large datasets. Don't go below 10.

---

### `nsample`
```
nsample = 200000
```

**What it does:** Number of post-burnin samples to record.

**Total MCMC iterations = burnin + nsample × sampfreq**

Example: `burnin=100000, sampfreq=50, nsample=200000` → total = 100,000 + 200,000 × 50 = **10.1 million iterations**

| Purpose | nsample | With sampfreq=50 | Total iterations |
|---------|---------|-----------------|-----------------|
| Quick test | 5,000 | 250,000 post-burnin | 350,000 |
| Preliminary | 20,000 | 1,000,000 post-burnin | 1,100,000 |
| **Publication** | **200,000** | **10,000,000 post-burnin** | **10,100,000** |
| High confidence | 500,000 | 25,000,000 post-burnin | 25,100,000 |

**How to know if you have enough:** Check ESS in Tracer. If ESS > 200 for all parameters, you have enough samples. If not, increase nsample (or decrease sampfreq, or both).

---

### `print`
```
print = 1
```

**What it does:** Controls what gets printed in the output.

| Value | Output |
|-------|--------|
| `1` | Node ages, 95% CIs, dated tree |
| `-1` | All of above PLUS individual branch rates (useful for rate analysis) |
| `2` | Detailed output (large files) |

**Recommendation:** Use `print = 1` for standard analyses. Use `print = -1` if you want to analyze rate variation across branches (e.g., for rate heatmaps).

---

## Complete Diagnostic Workflow

```
1. Run usedata=0 (prior only)
   → Check: Does the effective prior match your intended calibrations?
   → If not: Adjust calibrations or BDparas

2. Run usedata=3 (compute Hessian)
   → mv out.BV in.BV

3. Run usedata=2 (MCMC) — Chain 1 (seed=-1)
4. Run usedata=2 (MCMC) — Chain 2 (seed=-42)

5. Check convergence in Tracer:
   → ESS > 200 for all parameters? If no → run longer
   → Chains agree? If no → run much longer or check model
   → Trace plots stationary? If no → increase burnin

6. Sensitivity analyses:
   → clock=2 vs clock=3
   → Different calibration strategies
   → Different rgene_gamma and sigma2_gamma

7. Report:
   → Posterior means + 95% HPD
   → Clock model, substitution model, calibrations
   → ESS values and convergence evidence
   → Sensitivity analysis results
```
