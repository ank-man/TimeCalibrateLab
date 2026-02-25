# Complete Pipeline: OrthoFinder to MCMCTree Time-Calibrated Phylogeny

## Overview

```
Proteomes → OrthoFinder → Single-Copy Orthologs → MAFFT Alignment → 
PAL2NAL Codon Alignment → TrimAl → Concatenation → RAxML/IQ-TREE (optional) →
MCMCTree Hessian → MCMCTree MCMC → Tracer Convergence Check → Dated Tree
```

---

## Prerequisites

```bash
# Install required software
conda create -n timetree python=3.10
conda activate timetree

# Core tools
conda install -c bioconda orthofinder mafft trimal pal2nal iqtree raxml-ng
conda install -c bioconda paml    # contains mcmctree

# Utilities
pip install biopython
conda install -c bioconda amas    # for concatenation

# Visualization
# Install FigTree from: https://github.com/rambaut/figtree
# Install Tracer from: https://beast.community/tracer
```

---

## STEP 1: Run OrthoFinder

**Input:** Directory of proteome FASTA files (one per species)

```bash
mkdir -p project/proteomes
# Copy your proteome files into project/proteomes/
# Files should be named: Species1.fa, Species2.fa, etc.
# Headers should be clean (no special characters)

orthofinder -f project/proteomes/ \
    -t 16 \      # threads for sequence search
    -a 8 \       # threads for analysis
    -S diamond \  # faster than blast
    -M msa \      # multiple sequence alignment method for gene trees
    -o project/orthofinder_results
```

**Key outputs:**
```
orthofinder_results/
├── Species_Tree/
│   └── SpeciesTree_rooted_node_labels.txt    ← YOUR SPECIES TREE
├── Orthogroup_Sequences/                      ← all orthogroups
├── Single_Copy_Orthologue_Sequences/          ← single-copy (cleanest)
├── Orthogroups/
│   └── Orthogroups.tsv                        ← orthogroup membership
└── Comparative_Genomics_Statistics/
    └── Statistics_Overall.tsv                 ← summary stats
```

**Check:**
```bash
# How many single-copy orthologs?
ls project/orthofinder_results/Single_Copy_Orthologue_Sequences/*.fa | wc -l

# If too few (<50), consider relaxing to orthogroups present in all species
# but allow some missing data
```

---

## STEP 2: Prepare CDS Sequences (for codon alignments)

If you want codon-level alignments (recommended for more signal):

```bash
mkdir -p project/cds_by_orthogroup

# You need CDS FASTA files matching protein IDs
# Create a script to extract CDS sequences for each orthogroup

python3 << 'PYEOF'
import os
from Bio import SeqIO

# Build lookup: protein_id → CDS sequence
cds_lookup = {}
cds_dir = "project/cds_files"   # directory with per-species CDS files
for cds_file in os.listdir(cds_dir):
    if not cds_file.endswith(('.fa', '.fasta', '.fna')):
        continue
    for rec in SeqIO.parse(os.path.join(cds_dir, cds_file), "fasta"):
        cds_lookup[rec.id] = rec

# For each single-copy orthogroup, write matching CDS
sc_dir = "project/orthofinder_results/Single_Copy_Orthologue_Sequences"
out_dir = "project/cds_by_orthogroup"
os.makedirs(out_dir, exist_ok=True)

for og_file in os.listdir(sc_dir):
    if not og_file.endswith('.fa'):
        continue
    og_name = og_file.replace('.fa', '')
    cds_records = []
    for rec in SeqIO.parse(os.path.join(sc_dir, og_file), "fasta"):
        if rec.id in cds_lookup:
            cds_records.append(cds_lookup[rec.id])
        else:
            print(f"WARNING: No CDS for {rec.id} in {og_name}")
    if len(cds_records) == len(list(SeqIO.parse(os.path.join(sc_dir, og_file), "fasta"))):
        SeqIO.write(cds_records, os.path.join(out_dir, f"{og_name}_cds.fa"), "fasta")
PYEOF
```

---

## STEP 3: Align Protein Sequences

```bash
mkdir -p project/protein_alignments

for f in project/orthofinder_results/Single_Copy_Orthologue_Sequences/*.fa; do
    name=$(basename "$f" .fa)
    echo "Aligning $name..."
    mafft --auto --thread 4 "$f" > "project/protein_alignments/${name}_aligned.fa"
done

echo "Aligned $(ls project/protein_alignments/*.fa | wc -l) orthogroups"
```

---

## STEP 4: Back-Translate to Codon Alignments

```bash
mkdir -p project/codon_alignments

for f in project/protein_alignments/*_aligned.fa; do
    name=$(basename "$f" _aligned.fa)
    cds_file="project/cds_by_orthogroup/${name}_cds.fa"
    
    if [ -f "$cds_file" ]; then
        pal2nal.pl "$f" "$cds_file" -output fasta \
            > "project/codon_alignments/${name}_codon.fa" 2>/dev/null
        
        if [ $? -ne 0 ]; then
            echo "WARNING: PAL2NAL failed for $name (length mismatch?)"
        fi
    else
        echo "SKIP: No CDS file for $name, using protein alignment"
        cp "$f" "project/codon_alignments/${name}_codon.fa"
    fi
done
```

**If you don't have CDS:** Just use the protein alignments directly. Set `seqtype = 2` in MCMCTree.

---

## STEP 5: Trim Alignments

```bash
mkdir -p project/trimmed_alignments

for f in project/codon_alignments/*_codon.fa; do
    name=$(basename "$f" _codon.fa)
    trimal -in "$f" \
           -out "project/trimmed_alignments/${name}_trimmed.fa" \
           -automated1
done

# Remove any alignments that are too short after trimming
for f in project/trimmed_alignments/*.fa; do
    len=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    if [ "$len" -lt 100 ]; then
        echo "Removing short alignment: $f ($len bp)"
        rm "$f"
    fi
done

echo "$(ls project/trimmed_alignments/*.fa | wc -l) alignments passed trimming"
```

---

## STEP 6: Filter and Concatenate

```bash
mkdir -p project/concatenated

# Use AMAS to concatenate
python3 -m amas concat \
    -i project/trimmed_alignments/*_trimmed.fa \
    -f fasta \
    -d dna \         # use 'aa' for protein
    -p project/concatenated/partitions.txt \
    -t project/concatenated/concatenated.fa \
    --out-format fasta

# Check alignment dimensions
python3 -c "
from Bio import AlignIO
aln = AlignIO.read('project/concatenated/concatenated.fa', 'fasta')
print(f'Species: {len(aln)}')
print(f'Sites: {aln.get_alignment_length()}')
for rec in aln:
    print(f'  {rec.id}')
"
```

---

## STEP 7: Convert to PHYLIP Format

MCMCTree requires sequential PHYLIP format:

```bash
python3 << 'PYEOF'
from Bio import AlignIO

aln = AlignIO.read("project/concatenated/concatenated.fa", "fasta")
n_taxa = len(aln)
n_sites = aln.get_alignment_length()

with open("project/mcmctree/alignment.phy", "w") as f:
    f.write(f"  {n_taxa}  {n_sites}\n")
    for rec in aln:
        # MCMCTree needs names ≤30 chars, no spaces
        name = rec.id[:30].replace(" ", "_")
        # Write name padded to 31 chars, then sequence
        f.write(f"{name:<31}{str(rec.seq)}\n")

print(f"Written: {n_taxa} taxa, {n_sites} sites")
PYEOF
```

**Verify the format:**
```bash
head -2 project/mcmctree/alignment.phy
# Should look like:
#   8  125000
# Human                          ATCGATCG...
```

---

## STEP 8: Prepare the Calibrated Tree

### 8a. Get the species tree topology

```bash
mkdir -p project/mcmctree

# Copy the OrthoFinder species tree
cp project/orthofinder_results/Species_Tree/SpeciesTree_rooted_node_labels.txt \
   project/mcmctree/species_tree_raw.nwk

# View it
cat project/mcmctree/species_tree_raw.nwk
```

### 8b. Clean the tree and add calibrations

```python
#!/usr/bin/env python3
"""
prepare_calibrated_tree.py

Takes an OrthoFinder species tree, removes branch lengths,
and adds fossil calibration annotations for MCMCTree.

EDIT THE CALIBRATIONS SECTION BELOW FOR YOUR DATA.
"""

import re
import sys

# ============================================================
# EDIT THIS: Your input tree (Newick with or without branch lengths)
# ============================================================
input_tree_file = "project/mcmctree/species_tree_raw.nwk"
output_tree_file = "project/mcmctree/calibrated_tree.nwk"

# ============================================================
# EDIT THIS: Define your fossil calibrations
# 
# Format: { "node_defined_by_two_descendants": "calibration_string" }
#
# MCMCTree time units: 100 Myr (so 65 Ma = 0.65)
#
# Calibration types:
#   'B(tL, tU)'           — Uniform between tL and tU (soft bounds, 2.5% tails)
#   'B(tL, tU, pL, pU)'   — Uniform with custom tail probabilities
#   'L(tL, p, c, pL)'     — Truncated Cauchy (lower bound)
#   'U(tU, p, c, pU)'     — Truncated Cauchy (upper bound)  
#   'G(alpha, beta)'      — Gamma distribution
#   'S2N(loc, scale, shape)' — Skew-normal
#   'ST(loc, scale, shape, df)' — Skew-t
#
# Node is identified by the MRCA of two species names
# ============================================================

CALIBRATIONS = {
    # Example: root calibration (MRCA of all species)
    # Adjust species names to match YOUR data
    ("Human", "Lemur"):       "'B(0.55, 1.00)'",      # Root: 55-100 Ma
    ("Human", "Chimp"):       "'B(0.05, 0.10)'",      # Human-Chimp: 5-10 Ma
    ("Human", "OWMonkey"):    "'B(0.20, 0.35)'",      # Catarrhini: 20-35 Ma
    ("NWMonkey", "Lemur"):    "'B(0.55, 0.80)'",      # Strepsirrhini-Haplorhini: 55-80 Ma
}

# Root age (REQUIRED by MCMCTree separately)
ROOT_AGE = "'B(0.55, 1.00)'"

# ============================================================
# Processing logic (generally don't need to edit below)
# ============================================================

def read_tree(filepath):
    with open(filepath) as f:
        tree = f.read().strip()
    return tree

def strip_branch_lengths(tree):
    """Remove branch lengths and support values, keep topology."""
    # Remove branch lengths (:number)
    tree = re.sub(r':\d+\.?\d*([eE][+-]?\d+)?', '', tree)
    # Remove support values (numbers before closing parens or commas)
    # Be careful not to remove species names
    return tree

def get_species_list(tree):
    """Extract species names from Newick string."""
    # Remove everything except names and structural chars
    clean = re.sub(r':[^,\)]+', '', tree)  # remove branch lengths
    clean = re.sub(r'\)[^,\);]*', ')', clean)  # remove internal labels
    names = re.findall(r'[A-Za-z_][A-Za-z0-9_]*', clean)
    return names

def find_mrca_node(tree_str, sp1, sp2):
    """
    Find the closing parenthesis of the MRCA node of two species.
    Returns the position after the ')' where we should insert calibration.
    
    This is a simplified approach — for complex trees, use ete3.
    """
    # Parse parentheses to find MRCA
    # Build a nested structure
    depth = 0
    positions = []  # (depth, position_of_close_paren)
    taxa_at_depth = {}
    
    i = 0
    current_taxon = ""
    open_positions = []
    
    while i < len(tree_str):
        c = tree_str[i]
        if c == '(':
            depth += 1
            open_positions.append(i)
        elif c == ')':
            close_pos = i
            if depth not in taxa_at_depth:
                taxa_at_depth[depth] = set()
            positions.append((depth, close_pos))
            depth -= 1
        elif c == ',':
            if current_taxon.strip():
                if depth not in taxa_at_depth:
                    taxa_at_depth[depth] = set()
                taxa_at_depth[depth].add(current_taxon.strip())
            current_taxon = ""
        elif c not in ';\n\r':
            current_taxon += c
        i += 1
    
    return None  # Fallback — use ete3 for robust MRCA finding

# Simple approach: manually construct the tree string
tree = read_tree(input_tree_file)
tree = strip_branch_lengths(tree)

print(f"Input tree (cleaned): {tree}")
print(f"Species: {get_species_list(tree)}")
print(f"\nCalibrations to apply:")
for (sp1, sp2), cal in CALIBRATIONS.items():
    print(f"  MRCA({sp1}, {sp2}): {cal}")

# For robust MRCA finding, use ete3:
try:
    from ete3 import Tree
    
    t = Tree(tree, format=1)
    
    calibrated_tree = tree.rstrip(';')
    
    # Apply calibrations by finding MRCA nodes
    for (sp1, sp2), cal in CALIBRATIONS.items():
        try:
            mrca = t.get_common_ancestor(sp1, sp2)
            # Get the newick substring for this node
            node_newick = mrca.write(format=9).rstrip(';')
            # Insert calibration after the closing paren of this clade
            if node_newick in calibrated_tree:
                calibrated_tree = calibrated_tree.replace(
                    node_newick, 
                    node_newick + cal,
                    1  # only first occurrence
                )
                print(f"  Applied calibration for MRCA({sp1},{sp2})")
            else:
                print(f"  WARNING: Could not find subtree for MRCA({sp1},{sp2})")
        except Exception as e:
            print(f"  ERROR: {e}")
    
    calibrated_tree = calibrated_tree + ";"
    
except ImportError:
    print("\nWARNING: ete3 not installed. Writing tree without automatic calibration insertion.")
    print("You will need to manually add calibrations to the tree file.")
    print("Install ete3: pip install ete3")
    calibrated_tree = tree

# Count species
n_species = len(get_species_list(calibrated_tree))

# Write output
with open(output_tree_file, 'w') as f:
    f.write(f"{n_species} 1\n")
    f.write(calibrated_tree + "\n")

print(f"\nCalibrated tree written to: {output_tree_file}")
print(f"Number of species: {n_species}")
print(f"\nDon't forget to set RootAge = {ROOT_AGE} in the control file!")
```

### 8c. Manual alternative (most reliable)

If the script above is tricky, just manually edit the tree:

```
8 1
((((Human,Chimp)'B(0.05,0.10)',Gorilla),(Orangutan,Gibbon))'B(0.11,0.28)',(OWMonkey,(NWMonkey,Lemur)))'B(0.55,1.00)';
```

**Rules:**
- First line: `N_species N_trees`
- Calibration goes right after the `)` that closes the clade
- Times in units of 100 Myr
- No branch lengths in the tree
- Species names must exactly match the alignment

---

## STEP 9: MCMCTree — Compute the Hessian (Approximate Likelihood)

```bash
cd project/mcmctree

# Create the Hessian control file
cat > mcmctree_hessian.ctl << 'EOF'
          seed = -1
       seqfile = alignment.phy
      treefile = calibrated_tree.nwk
      mcmcfile = mcmc_hessian.txt
       outfile = out_hessian.txt

         ndata = 1
       seqtype = 0            * 0=nucleotides, 1=codons, 2=AAs
       usedata = 3            * 3 = compute Hessian (gradient + curvature)
         clock = 2            * 2 = independent rates (relaxed clock)
       RootAge = 'B(0.55,1.00)'

         model = 4            * 4 = HKY85
         alpha = 0.5          * gamma shape for among-site rate variation
         ncatG = 5            * number of gamma categories

     cleandata = 0            * 0 = keep gaps/ambiguities

       BDparas = 1 1 0.1     * birth rate, death rate, sampling fraction
   rgene_gamma = 2 20 1      * gamma prior on mean rate
  sigma2_gamma = 1 10 1      * gamma prior on rate variance (sigma^2)

      finetune = 1: .1 .1 .1 .1 .1 .1  * auto-tune step sizes

         print = 1
        burnin = 2000
      sampfreq = 2
       nsample = 5000
EOF

# Run it
mcmctree mcmctree_hessian.ctl

# This produces out.BV — the Hessian matrix
# Rename it for the next step
mv out.BV in.BV

echo "Hessian computed. in.BV is ready."
```

**What just happened:**
MCMCTree computed the gradient and Hessian (second derivative matrix) of the log-likelihood evaluated at the maximum likelihood branch lengths. This is a quadratic approximation of the likelihood surface that makes the subsequent MCMC orders of magnitude faster.

---

## STEP 10: MCMCTree — Run the MCMC (Main Analysis)

```bash
# Create the MCMC control file
cat > mcmctree_run.ctl << 'EOF'
          seed = -1
       seqfile = alignment.phy
      treefile = calibrated_tree.nwk
      mcmcfile = mcmc_samples.txt
       outfile = out_mcmctree.txt

         ndata = 1
       seqtype = 0            * 0=nucleotides
       usedata = 2            * 2 = approximate likelihood (uses in.BV)
         clock = 2            * 2 = independent rates relaxed clock
       RootAge = 'B(0.55,1.00)'

         model = 4            * HKY85
         alpha = 0.5          * gamma shape
         ncatG = 5            * gamma categories

     cleandata = 0

       BDparas = 1 1 0.1
   rgene_gamma = 2 20 1
  sigma2_gamma = 1 10 1

      finetune = 1: .1 .1 .1 .1 .1 .1

         print = 1
        burnin = 100000       * discard first 100k iterations
      sampfreq = 50           * save every 50th sample
       nsample = 200000       * collect 200k samples
                              * total iterations = 100000 + 200000*50 = 10.1 million
EOF

mcmctree mcmctree_run.ctl
```

---

## STEP 11: Run a Prior-Only Analysis (ESSENTIAL DIAGNOSTIC)

```bash
mkdir -p project/mcmctree/prior_only
cp calibrated_tree.nwk alignment.phy project/mcmctree/prior_only/

cat > project/mcmctree/prior_only/mcmctree_prior.ctl << 'EOF'
          seed = -1
       seqfile = alignment.phy
      treefile = calibrated_tree.nwk
      mcmcfile = mcmc_prior.txt
       outfile = out_prior.txt

         ndata = 1
       seqtype = 0
       usedata = 0            * 0 = IGNORE data, sample from prior only
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
EOF

cd project/mcmctree/prior_only
mcmctree mcmctree_prior.ctl
```

**Why this matters:** The "effective prior" (what MCMCTree actually uses) can differ from your specified calibrations because:
- The birth-death tree prior interacts with calibrations
- Node age ordering constraints truncate distributions
- Multiple calibrations constrain each other

Compare the prior-only node ages with your intended calibrations. If they differ substantially, adjust your calibrations or BDparas.

---

## STEP 12: Run Independent Chains for Convergence

```bash
# Chain 1 (already done above)
# Chain 2
mkdir -p project/mcmctree/run2
cp project/mcmctree/{in.BV,alignment.phy,calibrated_tree.nwk} project/mcmctree/run2/
cp project/mcmctree/mcmctree_run.ctl project/mcmctree/run2/
cd project/mcmctree/run2
sed -i 's/seed = -1/seed = -42/' mcmctree_run.ctl
sed -i 's/mcmc_samples/mcmc_samples_run2/' mcmctree_run.ctl
mcmctree mcmctree_run.ctl

# Chain 3 (optional but recommended)
mkdir -p project/mcmctree/run3
cp project/mcmctree/{in.BV,alignment.phy,calibrated_tree.nwk} project/mcmctree/run3/
cp project/mcmctree/mcmctree_run.ctl project/mcmctree/run3/
cd project/mcmctree/run3
sed -i 's/seed = -1/seed = -123/' mcmctree_run.ctl
sed -i 's/mcmc_samples/mcmc_samples_run3/' mcmctree_run.ctl
mcmctree mcmctree_run.ctl
```

---

## STEP 13: Check Convergence

### Using Tracer (GUI)

1. Open Tracer
2. File → Import Trace File → load all `mcmc_samples*.txt`
3. Check:
   - **ESS > 200** for all parameters (especially `t_n*` node ages)
   - **Trace plots** look like "hairy caterpillars" (good mixing)
   - **Posterior means agree** between independent runs (within 5%)

### Using R (command line)

```R
#!/usr/bin/env Rscript
# check_convergence.R

library(coda)

# Load chains
chain1 <- read.table("project/mcmctree/mcmc_samples.txt", header=TRUE)
chain2 <- read.table("project/mcmctree/run2/mcmc_samples_run2.txt", header=TRUE)

# Convert to mcmc objects
mc1 <- mcmc(chain1)
mc2 <- mcmc(chain2)

# ESS for each parameter
cat("=== ESS Chain 1 ===\n")
print(effectiveSize(mc1))
cat("\n=== ESS Chain 2 ===\n")
print(effectiveSize(mc2))

# Check if all ESS > 200
ess1 <- effectiveSize(mc1)
ess2 <- effectiveSize(mc2)
if (all(ess1 > 200) & all(ess2 > 200)) {
    cat("\n✓ All ESS values > 200. Chains appear converged.\n")
} else {
    cat("\n✗ Some ESS values < 200. Increase burnin or nsample.\n")
    cat("Low ESS parameters:\n")
    print(names(ess1[ess1 < 200]))
    print(names(ess2[ess2 < 200]))
}

# Compare posterior means
cat("\n=== Posterior Mean Comparison ===\n")
means1 <- colMeans(chain1)
means2 <- colMeans(chain2)
diff_pct <- abs(means1 - means2) / ((means1 + means2) / 2) * 100
cat("Max % difference between chains:", max(diff_pct, na.rm=TRUE), "%\n")
if (max(diff_pct, na.rm=TRUE) < 5) {
    cat("✓ Chains agree within 5%.\n")
} else {
    cat("✗ Chains disagree. Run longer or check calibrations.\n")
}

# Gelman-Rubin diagnostic
mc_list <- mcmc.list(mc1, mc2)
cat("\n=== Gelman-Rubin R-hat ===\n")
gr <- gelman.diag(mc_list, multivariate=FALSE)
print(gr$psrf)
if (all(gr$psrf[,1] < 1.05)) {
    cat("\n✓ All R-hat < 1.05. Good convergence.\n")
} else {
    cat("\n✗ Some R-hat > 1.05. Chains haven't converged.\n")
}
```

---

## STEP 14: Extract and Visualize Results

### The dated tree

MCMCTree writes `FigTree.tre` — a Nexus-format tree with:
- Posterior mean node ages
- 95% HPD intervals as node annotations

Open in FigTree:
```bash
figtree project/mcmctree/FigTree.tre
```

In FigTree:
1. Node Labels → Display: `height_95%_HPD`
2. Scale axis → check "Scale by factor" and enter 100 (to convert to Ma)
3. Time Scale → Offset = 0, Scale factor = 100

### Extract node ages programmatically

```python
#!/usr/bin/env python3
"""
extract_node_ages.py
Parse MCMCTree output to get node ages and HPD intervals.
"""

import re

def parse_mcmctree_output(filepath):
    """Parse out_mcmctree.txt for node age estimates."""
    with open(filepath) as f:
        content = f.read()
    
    # Find the section with node ages
    # MCMCTree prints: "Node  TimeUnit  95%_CI"
    results = []
    
    # Look for lines like: "t_n7     0.6523 (0.5512, 0.7891)"
    pattern = r'(t_n\d+)\s+([\d.]+)\s+\(([\d.]+),\s*([\d.]+)\)'
    matches = re.findall(pattern, content)
    
    for match in matches:
        node_id, mean, lower, upper = match
        results.append({
            'node': node_id,
            'mean_100myr': float(mean),
            'mean_ma': float(mean) * 100,
            'lower_95_ma': float(lower) * 100,
            'upper_95_ma': float(upper) * 100,
        })
    
    return results

results = parse_mcmctree_output("project/mcmctree/out_mcmctree.txt")

print(f"{'Node':<10} {'Mean (Ma)':<12} {'95% HPD (Ma)':<25}")
print("-" * 47)
for r in results:
    print(f"{r['node']:<10} {r['mean_ma']:<12.2f} [{r['lower_95_ma']:.2f}, {r['upper_95_ma']:.2f}]")
```

### Create a publication-quality figure

```R
#!/usr/bin/env Rscript
# plot_timetree.R

library(ape)
library(ggtree)
library(treeio)
library(ggplot2)

# Read MCMCTree output
tree <- read.beast("project/mcmctree/FigTree.tre")

# Scale to Ma (multiply by 100 if needed)
# Check your time units first!

p <- ggtree(tree, mrsd="2025-01-01") +   # adjust if needed
    theme_tree2() +
    geom_range(range='height_0.95_HPD', color='steelblue', alpha=0.4, size=2) +
    geom_nodelab(aes(label=round(height, 1)), size=2.5, hjust=-0.1) +
    geom_tiplab(size=3, fontface="italic") +
    scale_x_continuous(
        name = "Time (Ma)",
        breaks = seq(0, 200, 20)
    ) +
    coord_cartesian(xlim = c(-10, max(tree@data$height) * 1.3)) +
    theme(
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size=10)
    )

ggsave("project/timetree.pdf", p, width=10, height=6)
ggsave("project/timetree.png", p, width=10, height=6, dpi=300)

cat("Figure saved to project/timetree.pdf and .png\n")
```

---

## STEP 15: Sensitivity Analyses (For Publication)

### A. Compare clock models

Run with `clock = 2` (independent rates) AND `clock = 3` (autocorrelated rates):

```bash
mkdir -p project/mcmctree/clock3
cp project/mcmctree/{in.BV,alignment.phy,calibrated_tree.nwk,mcmctree_run.ctl} \
   project/mcmctree/clock3/
cd project/mcmctree/clock3
sed -i 's/clock = 2/clock = 3/' mcmctree_run.ctl
mcmctree mcmctree_run.ctl
```

Compare results between clock models. If they agree, report both. If they disagree, discuss which is more appropriate.

### B. Vary calibration strategies

1. Run with each calibration removed one at a time (leave-one-out)
2. Run with wider/narrower calibrations
3. Compare prior-only vs posterior to see which calibrations matter most

### C. Vary alignment data

1. Run with different subsets of genes (e.g., 50%, 75%, 100%)
2. Run with protein vs. nucleotide alignments
3. Assess whether results are stable across data subsets

---

## Quick Reference: MCMCTree Control File Template

```
*=============================================================
* MCMCTree Control File — Template with annotations
* Lines starting with * are comments
*=============================================================

     seed = -1              * random seed (-1 = use clock)
  seqfile = alignment.phy   * PHYLIP format alignment
 treefile = tree.nwk        * calibrated Newick tree
 mcmcfile = mcmc.txt        * MCMC sample output
  outfile = out.txt         * summary output

    ndata = 1               * number of partitions
  seqtype = 0               * 0:nucleotide 1:codon 2:amino acid
  usedata = 2               * 0:prior only 1:exact 2:approx 3:hessian
    clock = 2               * 1:strict 2:independent 3:autocorrelated
  RootAge = 'B(0.55,1.00)'  * root age prior (100Myr units)

    model = 4               * 0:JC 1:K80 4:HKY 7:GTR
    alpha = 0.5             * gamma shape (0=no rate var)
    ncatG = 5               * gamma categories

cleandata = 0               * 0:keep ambig 1:remove

  BDparas = 1 1 0.1         * birth death sampling
rgene_gamma = 2 20 1        * rate prior: mean=alpha/beta
sigma2_gamma = 1 10 1       * rate variance prior

 finetune = 1: .1 .1 .1 .1 .1 .1  * MCMC step sizes

    print = 1               * 1:ages -1:ages+rates
   burnin = 100000          * iterations to discard
 sampfreq = 50              * thinning interval
  nsample = 200000          * post-burnin samples
```

---

## Troubleshooting

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| MCMCTree crashes immediately | Mismatch between tree species names and alignment names | Make names identical |
| `Error: node age out of bounds` | Calibration conflicts or initial ages violate constraints | Check calibrations don't conflict; use wider bounds |
| Very low ESS (<50) | Poor mixing; chain hasn't converged | Increase burnin and nsample; adjust finetune |
| Prior-only and posterior identical | Data has no signal (alignment too short or too conserved) | Use more genes; check alignment quality |
| Posterior much wider than expected | High rate variation (sigma2 too large) or too few calibrations | Add calibrations; reduce sigma2_gamma prior mean |
| `in.BV` not found | Forgot to run usedata=3 first or forgot to rename out.BV | Run Hessian step; `mv out.BV in.BV` |
| Very different results between runs | Chain hasn't converged | Run much longer (10x burnin and nsample) |
