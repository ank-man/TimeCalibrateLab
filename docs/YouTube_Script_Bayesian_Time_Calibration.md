# YouTube Script: "De-Black-Boxing Bayesian Phylogenetic Time Calibration — From OrthoFinder to Dated Trees"

**Estimated Runtime:** 35–45 minutes
**Target Audience:** Graduate students, postdocs, and researchers who use MCMCTree but want to understand what it actually does.

---

## INTRO [0:00–2:00]

**[ON SCREEN: Animated phylogenetic tree with geological timescale]**

Hey everyone. Today we're going to do something that almost nobody does in comparative genomics tutorials — we're going to actually *understand* what happens when you time-calibrate a phylogenetic tree.

Most tutorials tell you: "Run OrthoFinder, get a tree, throw it into MCMCTree with some fossil calibrations, and out comes a dated phylogeny." But nobody tells you *what the math is doing*, *what each parameter means*, or *why your choices matter*.

Today, I'm going to:

1. Walk you through the **complete pipeline** from raw proteomes to a time-calibrated tree
2. **Open every black box** — every MCMCTree parameter, explained in plain English
3. Show you an **interactive simulator** where you can *see* how priors, data, and clock models shape your results

By the end, you'll never blindly copy-paste MCMCTree control files again.

Let's start.

---

## PART 1: THE BIG PICTURE [2:00–7:00]

**[ON SCREEN: Slide — "What is Time Calibration?"]**

### What Are We Actually Doing?

A phylogenetic tree from sequence data gives you **topology** — who is related to whom — and **branch lengths** in units of **substitutions per site**. That tells you *how much* evolution happened, but not *when*.

Time calibration converts those branch lengths into **absolute time** — millions of years. To do this, we need three ingredients:

**[ON SCREEN: Three-column diagram]**

1. **A tree with branch lengths** (from sequence data)
2. **A molecular clock model** (the relationship between substitutions and time)
3. **Fossil calibrations** (anchor points that connect the tree to geological time)

### The Core Equation

Everything in Bayesian time calibration comes back to one equation:

**[ON SCREEN: Large animated equation]**

```
Posterior ∝ Likelihood × Prior
```

Or more precisely:

```
P(times, rates | data, fossils) ∝ P(data | times, rates) × P(times | fossils) × P(rates)
```

Let me unpack this:

- **Posterior** — what we want: the probability distribution of divergence times given everything we know
- **Likelihood** — what the sequence data tells us: given these divergence times and substitution rates, how probable is the alignment we observed?
- **Prior on times** — what fossils tell us: before looking at any sequences, what do we believe about when nodes diverged?
- **Prior on rates** — our prior belief about how fast sequences evolve

The software uses **MCMC** (Markov Chain Monte Carlo) to sample from this posterior distribution, because we can't compute it analytically.

### Why Does This Matter?

**[ON SCREEN: Two trees side by side — different calibrations, wildly different dates]**

Because your choices — which fossils, which priors, which clock model — *dramatically* affect the output. This isn't just running a program. You're making scientific decisions at every step, and understanding the math helps you make better ones.

---

## PART 2: THE MOLECULAR CLOCK [7:00–14:00]

**[ON SCREEN: Slide — "The Molecular Clock"]**

### Strict Clock

The original idea from Zuckerkandl and Pauling in 1962: all lineages evolve at the same constant rate.

```
distance = rate × time
d = r × t
```

If we know the distance (from the alignment) and the rate, we can solve for time. Simple.

**But the strict clock is almost always wrong.** Different lineages evolve at different rates. Mice evolve faster than elephants. Parasites evolve faster than their hosts. Generation time, population size, metabolic rate — all affect substitution rate.

### Relaxed Clock

**[ON SCREEN: Animation — branches colored by rate, rates drawn from a distribution]**

Modern methods use a **relaxed molecular clock**. Each branch gets its own rate, but rates are not completely free — they're drawn from a shared distribution:

```
r_i ~ LogNormal(μ, σ²)
```

Where:
- **μ** = the log of the mean rate across the tree
- **σ²** = how much rates vary among branches

**Key insight:** σ controls everything.
- σ = 0 → strict clock (all rates identical)
- σ = 0.1 → mild rate variation
- σ = 1.0 → massive rate variation, node age uncertainty explodes

**[ON SCREEN: Interactive demo — slider for σ showing posterior widening]**

In MCMCTree, this is controlled by the `clock` parameter and the `rgene_gamma` prior. We'll come back to these.

### Independent Rates vs Autocorrelated Rates

Two philosophies:

1. **Independent rates (IR):** Each branch's rate is drawn independently. Used by BEAST's uncorrelated lognormal model. In MCMCTree: `clock = 2`.

2. **Autocorrelated rates (AR):** A branch's rate is correlated with its parent's rate. The idea: closely related lineages should have similar rates. In MCMCTree: `clock = 3`.

Which to use? It's an empirical question. The autocorrelated model often gives narrower credible intervals. Many researchers run both and compare.

---

## PART 3: FOSSIL CALIBRATIONS [14:00–21:00]

**[ON SCREEN: Slide — "Fossil Calibrations: The Anchors"]**

### Why We Need Fossils

Sequence data alone tells you **relative** divergence (branch lengths) but not **absolute** time. Fossils provide the absolute time anchor.

A fossil calibration says: "This node is *at least* X million years old, because we found a fossil of that age that belongs to this clade."

### The Critical Distinction: Minimum vs Maximum

**[ON SCREEN: Diagram — fossil at node with one-sided arrow]**

A fossil gives you a **minimum age** — the clade is *at least* as old as the oldest fossil in it. But the clade could be much older than its oldest fossil. We just haven't found older fossils.

**Maximum ages** are much harder to justify. Some approaches:
- First appearance of a containing group
- Biogeographic events
- Geological constraints (this island didn't exist before X Ma)

### Calibration Distributions in MCMCTree

MCMCTree offers several distributions. Let me explain each one and when to use it:

**[ON SCREEN: Animated distribution plots for each type]**

#### 1. Uniform: `B(tL, tU)`
```
B(0.5, 1.5)
```
- Flat probability between minimum (tL) and maximum (tU)
- Says: "I think the true age is somewhere between 50 and 150 Ma, and I have no preference within that range"
- **When to use:** When you have both a firm minimum and a firm maximum but no reason to prefer any age within the range
- **Problem:** The hard bounds can be unrealistic. A 0.025 tail probability is added by default.

#### 2. Skew-Normal / Skew-t: `S2N(location, scale, shape)`
```
S2N(1.0, 0.2, 5)
```
- Asymmetric bell curve
- **When to use:** When you have a best estimate and asymmetric uncertainty

#### 3. Gamma: `G(alpha, beta)`
```
G(2, 20)
```
- Right-skewed distribution
- **When to use:** As a soft maximum; most mass near the minimum, long tail toward older ages

#### 4. Truncated Cauchy: `L(tL, p, c, pL)`

The most commonly used in MCMCTree for node calibrations:
```
L(0.5, 0.1, 0.2, 0.02)
```
- `tL` = minimum bound (the fossil age)
- `p` = offset parameter (controls how far the mode is from the minimum)
- `c` = scale parameter (controls the width of the distribution)
- `pL` = tail probability below the minimum (default 0.01 = 1% chance node is younger than the fossil)

#### The "SoftBound" Approach

In MCMCTree, most calibration styles place a small probability (1–2.5%) outside the bounds. This acknowledges that fossils could be misdated or misplaced phylogenetically. The parameter controlling this is the tail probability, often 0.025 by default in uniform bounds.

### How Many Calibrations Do You Need?

**[ON SCREEN: Simulation results — 1 vs 3 vs 6 calibrations]**

More is generally better, but quality matters more than quantity:
- Minimum: 1 calibration (but results will be dominated by the prior)
- Recommended: 3–5 well-justified calibrations spread across the tree
- Each calibration constrains nearby nodes; distant nodes benefit less
- **Conflicting calibrations** are a serious problem — MCMCTree will warn you, but you need to resolve them biologically

### The Interaction Between Data and Calibrations

**[ON SCREEN: Interactive demo — alignment length slider showing posterior shifting]**

This is the key insight most people miss:

- **Short alignments** (few genes) → the likelihood is weak → **calibrations dominate**
- **Long alignments** (many genes, genome-scale) → the likelihood is strong → **data dominates**

With genome-scale datasets, your calibration priors might barely matter for well-supported nodes. But for poorly supported nodes, calibrations are still critical.

This is exactly what Bayes' theorem predicts:
- Weak data + strong prior = posterior ≈ prior
- Strong data + weak prior = posterior ≈ likelihood
- Strong data + strong prior = posterior is a compromise

---

## PART 4: THE LIKELIHOOD — WHAT SEQUENCE DATA CONTRIBUTES [21:00–25:00]

**[ON SCREEN: Slide — "The Likelihood Function"]**

### From Alignment to Likelihood

The likelihood asks: "Given these divergence times and rates, how probable is the sequence alignment we observed?"

For two sequences separated by branch length `d`:

Under the JC69 model (simplest substitution model):
```
P(same base after distance d) = 1/4 + 3/4 × e^(-4d/3)
P(different base after distance d) = 1/4 - 1/4 × e^(-4d/3)
```

Where `d = rate × time`.

For an alignment of L sites where S are identical and D are different:
```
log-Likelihood = S × log(P(same)) + D × log(P(different))
```

**[ON SCREEN: Animation showing likelihood curve narrowing as L increases]**

More sites = sharper likelihood peak = more precise time estimates.

### The Approximate Likelihood in MCMCTree

MCMCTree doesn't compute the full likelihood at each MCMC step — that would be too slow for genomic data. Instead, it uses Thorne et al.'s **approximate likelihood**:

1. First, compute the gradient (first derivative) and Hessian (second derivative) of the log-likelihood at the maximum likelihood branch lengths
2. Then approximate the log-likelihood as a Taylor expansion around those values

This is why MCMCTree has **two steps**: first `usedata = 3` to compute the Hessian, then `usedata = 2` to run the MCMC with the approximation.

This approximation is very accurate for most datasets and makes MCMCTree vastly faster than BEAST for large datasets.

---

## PART 5: MCMC — HOW WE SAMPLE THE POSTERIOR [25:00–29:00]

**[ON SCREEN: Slide — "MCMC Sampling"]**

### The Problem

We want to know the posterior distribution:
```
P(times, rates | data, fossils)
```

But we can't compute it directly. The integral in the denominator (the marginal likelihood) is intractable.

### The Solution: Metropolis-Hastings

MCMC generates **samples** from the posterior without computing it directly:

**[ON SCREEN: Step-by-step animation of MCMC proposals]**

1. Start with initial node ages and rates
2. **Propose** a small change (e.g., nudge one node age by a random amount)
3. Compute the **acceptance ratio**:
   ```
   α = min(1, P(proposed) / P(current))
   ```
   Where P includes both the likelihood and the prior
4. Accept with probability α, reject otherwise
5. Record the current state
6. Repeat thousands of times

After enough iterations, the samples approximate the posterior distribution.

### Key MCMC Concepts

**Burn-in:** The first N samples are discarded because the chain hasn't converged yet. MCMCTree parameter: `burnin`.

**Sampling interval (thinning):** Only record every k-th sample to reduce autocorrelation. MCMCTree parameter: `sampfreq`.

**Number of samples:** How many post-burnin samples to collect. MCMCTree parameter: `nsample`.

**Convergence:** How do you know the chain has converged? Run it multiple times from different starting points. If they give similar posteriors, you're probably converged. Use Tracer to check ESS (Effective Sample Size) — you want ESS > 200 for each parameter.

**Proposal tuning:** MCMCTree automatically tunes step sizes, but you can adjust `finetune` parameters. We'll cover these in the parameter section.

---

## PART 6: COMPLETE PIPELINE — ORTHOFINDER TO MCMCTREE [29:00–38:00]

**[ON SCREEN: Slide — "The Complete Pipeline"]**

Now let's do the whole thing, step by step. I'll show the code, the logic, and the parameters.

### Step 0: Input Data

You need:
- Proteomes (FASTA files) for each species
- Some knowledge of fossil calibrations for your group

### Step 1: OrthoFinder

OrthoFinder finds orthologs, infers gene trees, and infers a species tree.

```bash
orthofinder -f proteomes/ -t 16 -a 8 -S diamond
```

**Outputs we need:**
- `Species_Tree/SpeciesTree_rooted_node_labels.txt` — the species tree (topology)
- `Orthogroup_Sequences/` — orthogroup protein FASTA files
- `Single_Copy_Orthologue_Sequences/` — single-copy orthologs (cleanest data)

### Step 2: Select and Align Orthologs

We want single-copy orthologs for time calibration. These are genes present in exactly one copy in every species.

```bash
# Align each single-copy orthogroup with MAFFT
for f in Single_Copy_Orthologue_Sequences/*.fa; do
    name=$(basename "$f" .fa)
    mafft --auto "$f" > alignments/${name}_aligned.fa
done
```

### Step 3: Back-Translate to Codons (Optional but Recommended)

Protein alignments are more reliable, but nucleotide data contains more information for divergence estimation. Use PAL2NAL to back-translate:

```bash
for f in alignments/*_aligned.fa; do
    name=$(basename "$f" _aligned.fa)
    pal2nal.pl "$f" cds/${name}.fa -output fasta > codon_alignments/${name}_codon.fa
done
```

### Step 4: Trim Alignments

Remove poorly aligned regions:

```bash
for f in codon_alignments/*.fa; do
    name=$(basename "$f" .fa)
    trimal -in "$f" -out trimmed/${name}_trimmed.fa -automated1
done
```

### Step 5: Concatenate Alignments

MCMCTree can handle partitioned data, but a concatenated alignment is simplest:

```bash
# Using AMAS
python3 AMAS.py concat -i trimmed/*_trimmed.fa -f fasta -d dna -p partitions.txt
```

This gives you:
- `concatenated.fa` — the supermatrix
- `partitions.txt` — partition boundaries

### Step 6: Convert to PHYLIP Format

MCMCTree requires PHYLIP format:

```bash
# Convert FASTA to sequential PHYLIP
# You can use a script or goalign:
goalign reformat phylip -i concatenated.fa -o alignment.phy
```

The format must be:
```
 Nspecies Nsites
Species1   ATCGATCG...
Species2   ATCGATCG...
```

**Important:** Species names must match between the tree and alignment exactly.

### Step 7: Prepare the Tree

Take the OrthoFinder species tree and:
1. Root it correctly
2. Remove branch lengths (MCMCTree estimates these)
3. Add calibration labels

**Starting tree from OrthoFinder:**
```
((((Human,Chimp),Gorilla),(Orangutan,Gibbon)),(OWMonkey,(NWMonkey,Lemur)));
```

**Calibrated tree for MCMCTree:**
```
8 1
((((Human,Chimp)'B(0.05,0.10)',Gorilla)'B(0.08,0.14)',(Orangutan,Gibbon))'B(0.11,0.28)',(OWMonkey,(NWMonkey,Lemur)'B(0.34,0.65)'))'B(0.55,1.00)';
```

**Format explanation:**
- First line: `N_species  N_trees` (8 species, 1 tree)
- Calibrations go after the closing parenthesis of the node they calibrate
- Times are in units of **100 Myr** by default (so 0.05 = 5 Ma, 1.00 = 100 Ma)

**Calibration syntax:**
- `'B(tL, tU)'` — uniform with soft bounds between tL and tU (in 100 Myr units)
- `'L(tL, p, c, pL)'` — truncated Cauchy
- `'G(alpha, beta)'` — gamma distribution
- `'S2N(location, scale, shape)'` — skew-normal

### Step 8: MCMCTree — Step 1: Compute the Hessian

**Create the control file `mcmctree_hessian.ctl`:**

```
          seed = -1
       seqfile = alignment.phy
      treefile = calibrated_tree.nwk
      mcmcfile = mcmc_hessian.txt
       outfile = out_hessian.txt

         ndata = 1
       seqtype = 0
       usedata = 3
         clock = 2
       RootAge = 'B(0.55,1.00)'

         model = 4
         alpha = 0.5
         ncatG = 5

   cleandata = 0

     BDparas = 1 1 0.1
 rgene_gamma = 2 20 1
sigma2_gamma = 1 10 1

      finetune = 1: .1 .1 .1 .1 .1 .1

         print = 1
        burnin = 2000
      sampfreq = 10
       nsample = 20000
```

```bash
mcmctree mcmctree_hessian.ctl
```

This produces `out.BV` (the Hessian file — gradient and curvature of the likelihood surface).

**Then rename it:**
```bash
mv out.BV in.BV
```

### Step 9: MCMCTree — Step 2: Run the MCMC

**Create `mcmctree_run.ctl`:**

```
          seed = -1
       seqfile = alignment.phy
      treefile = calibrated_tree.nwk
      mcmcfile = mcmc_samples.txt
       outfile = out_mcmctree.txt

         ndata = 1
       seqtype = 0
       usedata = 2
         clock = 2
       RootAge = 'B(0.55,1.00)'

         model = 4
         alpha = 0.5
         ncatG = 5

   cleandata = 0

     BDparas = 1 1 0.1
 rgene_gamma = 2 20 1
sigma2_gamma = 1 10 1

      finetune = 1: .1 .1 .1 .1 .1 .1

         print = 1
        burnin = 50000
      sampfreq = 50
       nsample = 100000
```

```bash
mcmctree mcmctree_run.ctl
```

### Step 10: Check Convergence

Run at least 2 independent chains (different seeds):
```bash
# Run 1
mkdir run1 && cp mcmctree_run.ctl in.BV alignment.phy calibrated_tree.nwk run1/
cd run1 && mcmctree mcmctree_run.ctl

# Run 2
mkdir run2 && cp mcmctree_run.ctl in.BV alignment.phy calibrated_tree.nwk run2/
cd run2 && sed -i 's/seed = -1/seed = -2/' mcmctree_run.ctl && mcmctree mcmctree_run.ctl
```

**Check with Tracer:**
- Load both `mcmc_samples.txt` files
- Check ESS > 200 for all parameters
- Visually confirm trace plots look like "hairy caterpillars" (good mixing)
- Compare posterior means between runs (should be very similar)

### Step 11: Summarize Results

MCMCTree outputs:
- `FigTree.tre` — tree with posterior mean ages and 95% HPD intervals (open in FigTree)
- `out_mcmctree.txt` — text summary with node ages
- `mcmc_samples.txt` — raw MCMC samples for further analysis

---

## PART 7: EVERY MCMCTREE PARAMETER EXPLAINED [38:00–50:00]

**[ON SCREEN: Slide — "Opening the Black Box: Every MCMCTree Parameter"]**

Here is every parameter in the MCMCTree control file, explained in plain language with guidance on what to set.

### I/O Parameters

| Parameter | What It Does | Typical Value |
|-----------|-------------|---------------|
| `seed` | Random number seed. Use -1 for clock-based seed, or set a specific integer for reproducibility. Different seeds for independent runs. | `-1` |
| `seqfile` | Path to the PHYLIP-format sequence alignment. | `alignment.phy` |
| `treefile` | Path to the Newick tree with calibration annotations. | `tree.nwk` |
| `mcmcfile` | Output file for MCMC samples. Each row = one sample of all parameters. | `mcmc.txt` |
| `outfile` | Output file for summary statistics and the dated tree. | `out.txt` |
| `ndata` | Number of data partitions in the alignment file. If you have one concatenated alignment, use 1. For partitioned analysis, stack multiple PHYLIP blocks in the seqfile. | `1` |

### Data Handling

| Parameter | What It Does | Options |
|-----------|-------------|---------|
| `seqtype` | Type of sequence data. | `0` = nucleotides, `1` = codons, `2` = amino acids |
| `usedata` | How to use the alignment. **This is the most important parameter.** | `0` = ignore data (sample from prior only — use this for prior predictive checks!), `1` = exact likelihood (slow, for small datasets), `2` = approximate likelihood using Hessian from in.BV (fast, standard for large datasets), `3` = compute the Hessian and write out.BV (run this first, then switch to 2) |
| `cleandata` | How to handle ambiguous sites and alignment gaps. | `0` = keep all sites (use ambiguity codes), `1` = remove all sites with ambiguities or gaps |

### Clock Model

| Parameter | What It Does | Options |
|-----------|-------------|---------|
| `clock` | The molecular clock model. **Major scientific decision.** | `1` = strict clock (one rate for entire tree — rarely appropriate), `2` = independent rates (each branch draws its rate independently from a lognormal — most common choice), `3` = autocorrelated rates (child branch rate is correlated with parent — biologically plausible for body size/generation time effects) |

**How to choose:** Run both clock=2 and clock=3, compare results. If they agree, report both. If they disagree, think about which model is more biologically appropriate for your group.

### Substitution Model

| Parameter | What It Does | Options |
|-----------|-------------|---------|
| `model` | Nucleotide substitution model. | `0` = JC69 (equal rates, equal frequencies), `1` = K80 (transition/transversion), `2` = F81, `3` = F84, `4` = HKY85 (most commonly used — different transition/transversion rates + empirical base frequencies), `5` = T92, `6` = TN93, `7` = REV (GTR — most parameter-rich) |
| `alpha` | Shape parameter of the gamma distribution for among-site rate variation. Lower = more rate heterogeneity among sites. | `0` = no rate variation (all sites evolve at same rate), `0.5` = substantial variation (common default), `> 1` = mild variation. Typically estimated from preliminary phylogenetic analysis (RAxML, IQ-TREE). |
| `ncatG` | Number of discrete gamma categories for rate variation. More categories = better approximation but slower. | `4` or `5` (standard) |

### Prior on Divergence Times: Birth-Death Process

| Parameter | What It Does |
|-----------|-------------|
| `BDparas` | `lambda mu rho` — parameters of the birth-death process used as the prior on divergence times. `lambda` = birth (speciation) rate, `mu` = death (extinction) rate, `rho` = sampling fraction (proportion of extant species included). |

**Typical values:** `BDparas = 1 1 0.1`

**What this means biologically:** The birth-death prior says that the tree of divergence times was generated by a process where lineages speciate at rate λ and go extinct at rate μ, and we've sampled a fraction ρ of the extant diversity.

**Practical guidance:**
- `lambda ≈ mu` → many lineages are born and die, tree is "bushy"
- `rho` = (species in your tree) / (total species in the clade). If you have 8 primates out of ~500, rho ≈ 0.016
- This prior has **relatively little effect** when you have good calibrations. Don't stress too much about it, but don't set absurd values.

### Prior on Rates: rgene_gamma

| Parameter | What It Does |
|-----------|-------------|
| `rgene_gamma` | `alpha beta something` — Gamma-Dirichlet prior on substitution rates across loci. `alpha` and `beta` parameterize a gamma distribution: mean = alpha/beta, variance = alpha/beta². The third number (1 or 2) controls the concentration parameter. |

**This sets your prior belief about substitution rates.**

```
rgene_gamma = 2 20 1
```

**Mean rate:** alpha/beta = 2/20 = 0.1 substitutions per site per 100 Myr = 0.001 per site per Myr = 10⁻³

**How to set it:**
1. Look up substitution rates for your group in the literature
2. Or estimate from a preliminary analysis (e.g., r8s, BEAST)
3. Convert to the right time units (MCMCTree usually works in 100 Myr)
4. Set alpha/beta to give that mean, with alpha controlling how informative the prior is:
   - alpha = 1 → exponential (vague)
   - alpha = 2 → mildly informative (good default)
   - alpha = 10 → quite informative (only if you're confident)

**Example calculations:**
- Mammals, nuclear DNA: ~0.005 subst/site/Myr → in 100 Myr units: 0.5 → `rgene_gamma = 2 4 1`
- Plants, chloroplast DNA: ~0.001 subst/site/Myr → in 100 Myr units: 0.1 → `rgene_gamma = 2 20 1`
- Fast-evolving virus: ~0.001 subst/site/year → need different time scale entirely

### Prior on Rate Variation: sigma2_gamma

| Parameter | What It Does |
|-----------|-------------|
| `sigma2_gamma` | `alpha beta something` — Gamma prior on σ², the variance of the log-normal distribution of rates among branches. Only relevant for `clock = 2` or `clock = 3`. |

```
sigma2_gamma = 1 10 1
```

**Mean σ²:** alpha/beta = 1/10 = 0.1

**What σ² means:**
- σ² = 0 → strict clock (no rate variation)
- σ² = 0.01 → very low rate variation (near-clock-like)
- σ² = 0.1 → moderate rate variation (common for many datasets)
- σ² = 1.0 → extreme rate variation (divergence times become very uncertain)

**How to set it:**
- Start with `sigma2_gamma = 1 10 1` (mean 0.1, moderately informative)
- If your group has very clock-like behavior: `sigma2_gamma = 2 40 1` (mean 0.05)
- If rates are highly variable: `sigma2_gamma = 1 4.5 1` (mean ~0.22)

### MCMC Tuning Parameters

| Parameter | What It Does | Guidance |
|-----------|-------------|----------|
| `finetune` | `auto: t1 t2 t3 t4 t5 t6` — Step sizes for MCMC proposals. First value: 1 = auto-adjust during burnin, 0 = fixed. Following values are initial step sizes for: node ages, rates, mixing, paras, SPR?, expansion? | `1: .1 .1 .1 .1 .1 .1` (auto-tune is almost always fine) |
| `burnin` | Number of MCMC iterations to discard at the start. The chain needs time to find the high-probability region. | `50000` minimum for real analyses. `200000` for large datasets. Check with Tracer. |
| `sampfreq` | Record every k-th sample after burnin. Reduces autocorrelation in saved samples. | `50` to `100`. Higher = less autocorrelated but slower to accumulate samples. |
| `nsample` | Number of post-burnin samples to collect. Total iterations = burnin + nsample × sampfreq. | `100000` for preliminary; `200000–500000` for publication. |
| `print` | Controls output verbosity. | `1` = standard output (dates, HPD intervals). `-1` = also print rate estimates for each branch. |

### Root Age

| Parameter | What It Does |
|-----------|-------------|
| `RootAge` | Prior on the age of the root node. Can be a calibration distribution just like node calibrations. **Required even if the root node already has a calibration in the tree file** — MCMCTree uses this as an overall time scale. |

```
RootAge = 'B(0.55, 1.00)'
```
This means: root is between 55 and 100 Ma (in 100 Myr units: 0.55 to 1.00).

---

## PART 8: COMMON MISTAKES AND HOW TO AVOID THEM [50:00–55:00]

**[ON SCREEN: Slide — "Common Pitfalls"]**

### 1. Wrong Time Units
MCMCTree works in units of 100 Myr by default. If your fossil is 65 Ma old, enter 0.65, not 65. **This is the #1 mistake.**

### 2. Never Running Prior-Only Analysis
Set `usedata = 0` and run MCMCTree. Compare the effective prior to your intended prior. If they differ substantially, your calibrations may be interacting in unexpected ways.

### 3. Insufficient Burn-in / Samples
Don't trust results if ESS < 200. Increase burnin and nsample. For genomic data, burnin = 200,000 and nsample = 500,000 is reasonable.

### 4. Not Running Multiple Chains
Always run at least 2 independent chains. If posteriors differ, you haven't converged.

### 5. Ignoring the Effective Prior
The birth-death tree prior interacts with your calibrations to create an "effective prior" that may differ from what you specified. Always check by running with `usedata = 0`.

### 6. Using a Strict Clock When It's Inappropriate
Almost all real datasets violate the strict clock. Use `clock = 2` or `clock = 3` unless you have strong evidence for clock-like evolution.

### 7. Recycling Someone Else's Calibrations Without Checking
Fossil calibrations from other studies may not apply to your taxon sampling. Always verify that the fossil actually falls within the clade it's calibrating in your specific tree.

---

## PART 9: INTERPRETING RESULTS [55:00–58:00]

**[ON SCREEN: FigTree screenshot with annotated tree]**

### What MCMCTree Gives You

1. **Posterior mean ages** for each node
2. **95% HPD (Highest Posterior Density) intervals** — the shortest interval containing 95% of the posterior probability
3. **Rate estimates** for each branch (if `print = -1`)

### How to Report Results

In your paper:
- Report posterior means AND 95% HPD intervals
- State the clock model, substitution model, calibration distributions
- State the number of genes/sites, burnin, samples, and ESS values
- Ideally: report results under both IR and AR clock models

### What Wide HPD Intervals Mean

A wide interval (e.g., 45–120 Ma) means:
- Not enough data (short alignment)
- Conflicting signals
- Calibration is too vague
- Rate variation is too high

This is not a failure — it's honesty about uncertainty.

---

## OUTRO [58:00–60:00]

**[ON SCREEN: TimeCalibrateLab simulator demo]**

So that's the complete picture. Every step from proteomes to dated trees, every parameter explained.

The key takeaway: Bayesian time calibration is not a black box. It's the equation:

```
Posterior ∝ Likelihood × Prior
```

Your sequence data provides the likelihood. Your fossils and clock model provide the prior. MCMCTree combines them into the posterior.

If you want to play with these concepts interactively, check out TimeCalibrateLab — the simulator I built that lets you adjust priors, alignment length, and clock models and *see* how the posterior changes in real time. Link in the description.

Thanks for watching. If this helped you, like, subscribe, and drop a comment with your organism group — I'd love to know what you're dating.

---

## DESCRIPTION BOX TEXT

```
🧬 De-Black-Boxing Bayesian Phylogenetic Time Calibration

Complete pipeline from OrthoFinder to MCMCTree dated phylogeny.
Every parameter explained. No black boxes.

📋 TIMESTAMPS:
0:00 - Introduction
2:00 - The Big Picture: What is Time Calibration?
7:00 - The Molecular Clock (Strict vs Relaxed)
14:00 - Fossil Calibrations (Distributions, Strategy)
21:00 - The Likelihood Function
25:00 - MCMC Sampling Explained
29:00 - Complete Pipeline: OrthoFinder → MCMCTree
38:00 - Every MCMCTree Parameter Explained
50:00 - Common Mistakes
55:00 - Interpreting Results

🔗 LINKS:
- TimeCalibrateLab Interactive Simulator: [link]
- Pipeline Scripts: [GitHub link]
- MCMCTree manual: http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf
- OrthoFinder: https://github.com/davidemms/OrthoFinder
- Tracer: https://beast.community/tracer

📚 KEY REFERENCES:
- Yang (2007) MCMCTree — doi:10.1093/molbev/msl150
- dos Reis et al. (2016) Bayesian molecular clock dating — doi:10.1016/j.tig.2015.10.007
- Emms & Kelly (2019) OrthoFinder — doi:10.1186/s13059-019-1832-y

#phylogenetics #evolution #bioinformatics #timetree #bayesian #mcmctree
```
