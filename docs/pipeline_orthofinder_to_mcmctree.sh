#!/usr/bin/env bash
#===============================================================================
# Pipeline: OrthoFinder → MCMCTree Time-Calibrated Phylogeny
#
# A complete, reproducible pipeline for building a time-calibrated species tree
# from proteome FASTA files using OrthoFinder for ortholog detection,
# MAFFT for alignment, TrimAl for trimming, and MCMCTree for divergence dating.
#
# Usage:
#   1. Edit the CONFIGURATION section below
#   2. chmod +x pipeline_orthofinder_to_mcmctree.sh
#   3. ./pipeline_orthofinder_to_mcmctree.sh
#
# Prerequisites:
#   conda install -c bioconda orthofinder mafft trimal pal2nal paml iqtree
#   pip install biopython ete3
#
# Author: TimeCalibrateLab
# License: MIT
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION — EDIT THIS SECTION
#===============================================================================

# Project name and directory
PROJECT="timetree_project"
WORKDIR="$(pwd)/${PROJECT}"

# Input data
PROTEOME_DIR="$(pwd)/proteomes"       # Directory containing *.fa proteome files
CDS_DIR=""                             # Optional: directory with CDS *.fa files
                                       # Leave empty if you only have proteins

# Threads
THREADS=16
ANALYSIS_THREADS=8

# Sequence type for MCMCTree
# "nucleotide" if you have CDS and want codon back-translation
# "protein" if you only have protein sequences
SEQ_TYPE="nucleotide"

# MCMCTree settings
CLOCK_MODEL=2                          # 1=strict, 2=independent rates, 3=autocorrelated
SUBST_MODEL=4                          # 0=JC, 4=HKY, 7=GTR (for nucleotides)
GAMMA_ALPHA=0.5                        # Among-site rate variation shape
GAMMA_CATS=5                           # Discrete gamma categories

# Rate prior: Gamma(alpha, beta) → mean rate = alpha/beta (in 100 Myr units)
# Typical mammalian nuclear: alpha=2, beta=5 (mean=0.4 ≈ 0.004 subst/site/Myr)
# Typical plant nuclear:     alpha=2, beta=20 (mean=0.1 ≈ 0.001 subst/site/Myr)
RGENE_ALPHA=2
RGENE_BETA=20

# Rate variance prior: Gamma(alpha, beta) → mean sigma² = alpha/beta
# 1 10 → mean 0.1 (moderate rate variation, good default)
SIGMA2_ALPHA=1
SIGMA2_BETA=10

# Birth-death tree prior
BD_LAMBDA=1                            # Speciation rate
BD_MU=1                                # Extinction rate
BD_RHO=0.1                             # Sampling fraction (species in tree / total in clade)

# MCMC chain settings
BURNIN=100000
SAMPFREQ=50
NSAMPLE=200000

# Root age prior (in 100 Myr units)
# Example: 'B(0.55, 1.00)' means uniform 55-100 Ma
ROOT_AGE="'B(0.55, 1.00)'"

# Minimum alignment length after trimming (bp)
MIN_ALN_LENGTH=100

# Number of independent MCMC chains for convergence check
N_CHAINS=2

#===============================================================================
# CALIBRATIONS — EDIT THIS
#
# Define your fossil calibrations here.
# These will be manually inserted into the tree file.
# After running Steps 1-7, you must MANUALLY edit the tree file
# (see STEP 8 below for instructions).
#===============================================================================

# Calibration notes (for reference — actual insertion is manual):
# Node: MRCA(Human,Chimp)    → 'B(0.05, 0.10)'   # 5-10 Ma
# Node: MRCA(Human,OWMonkey) → 'B(0.20, 0.35)'   # 20-35 Ma  
# Node: Root                 → 'B(0.55, 1.00)'   # 55-100 Ma

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

log_step() {
    echo ""
    echo "================================================================"
    echo "  STEP $1: $2"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "================================================================"
    echo ""
}

log_info() {
    echo "[INFO] $1"
}

log_warn() {
    echo "[WARN] $1" >&2
}

log_error() {
    echo "[ERROR] $1" >&2
}

check_tool() {
    if ! command -v "$1" &>/dev/null; then
        log_error "$1 is not installed or not in PATH"
        log_error "Install with: $2"
        exit 1
    fi
}

count_files() {
    ls "$1" 2>/dev/null | wc -l
}

#===============================================================================
# PREFLIGHT CHECKS
#===============================================================================

log_step "0" "Preflight checks"

check_tool "orthofinder" "conda install -c bioconda orthofinder"
check_tool "mafft" "conda install -c bioconda mafft"
check_tool "trimal" "conda install -c bioconda trimal"
check_tool "mcmctree" "conda install -c bioconda paml"

if [ "${SEQ_TYPE}" = "nucleotide" ] && [ -n "${CDS_DIR}" ]; then
    check_tool "pal2nal.pl" "conda install -c bioconda pal2nal"
fi

# Check Python modules
python3 -c "from Bio import SeqIO" 2>/dev/null || {
    log_error "BioPython not installed. Run: pip install biopython"
    exit 1
}

# Check input
if [ ! -d "${PROTEOME_DIR}" ]; then
    log_error "Proteome directory not found: ${PROTEOME_DIR}"
    exit 1
fi

N_PROTEOMES=$(count_files "${PROTEOME_DIR}/*.fa*")
if [ "${N_PROTEOMES}" -eq 0 ]; then
    log_error "No FASTA files (*.fa, *.fasta) found in ${PROTEOME_DIR}"
    exit 1
fi
log_info "Found ${N_PROTEOMES} proteome files in ${PROTEOME_DIR}"

# Create working directory
mkdir -p "${WORKDIR}"
log_info "Working directory: ${WORKDIR}"

#===============================================================================
# STEP 1: Run OrthoFinder
#===============================================================================

log_step "1" "Running OrthoFinder"

ORTHOFINDER_OUT="${WORKDIR}/orthofinder_results"

if [ -d "${ORTHOFINDER_OUT}" ]; then
    log_info "OrthoFinder output already exists, skipping..."
else
    orthofinder \
        -f "${PROTEOME_DIR}" \
        -t "${THREADS}" \
        -a "${ANALYSIS_THREADS}" \
        -S diamond \
        -M msa \
        -o "${ORTHOFINDER_OUT}"
fi

# Find the actual results directory (OrthoFinder nests under Results_*)
RESULTS_DIR=$(find "${ORTHOFINDER_OUT}" -maxdepth 1 -name "Results_*" -type d | head -1)
if [ -z "${RESULTS_DIR}" ]; then
    # Maybe it's directly in the output
    RESULTS_DIR="${ORTHOFINDER_OUT}"
fi

SC_DIR="${RESULTS_DIR}/Single_Copy_Orthologue_Sequences"
SPECIES_TREE="${RESULTS_DIR}/Species_Tree/SpeciesTree_rooted_node_labels.txt"

if [ ! -d "${SC_DIR}" ]; then
    log_error "Single_Copy_Orthologue_Sequences not found at: ${SC_DIR}"
    log_error "Check OrthoFinder output structure."
    exit 1
fi

N_SC=$(count_files "${SC_DIR}/*.fa")
log_info "Found ${N_SC} single-copy orthogroups"
log_info "Species tree: ${SPECIES_TREE}"

#===============================================================================
# STEP 2: Align protein sequences
#===============================================================================

log_step "2" "Aligning protein sequences with MAFFT"

PROTEIN_ALN_DIR="${WORKDIR}/protein_alignments"
mkdir -p "${PROTEIN_ALN_DIR}"

for f in "${SC_DIR}"/*.fa; do
    name=$(basename "$f" .fa)
    outfile="${PROTEIN_ALN_DIR}/${name}_aligned.fa"
    if [ -f "${outfile}" ]; then
        continue
    fi
    log_info "Aligning ${name}..."
    mafft --auto --thread "${ANALYSIS_THREADS}" "$f" > "${outfile}" 2>/dev/null
done

N_ALIGNED=$(count_files "${PROTEIN_ALN_DIR}/*_aligned.fa")
log_info "Aligned ${N_ALIGNED} orthogroups"

#===============================================================================
# STEP 3: Back-translate to codons (if CDS available)
#===============================================================================

FINAL_ALN_DIR="${PROTEIN_ALN_DIR}"

if [ "${SEQ_TYPE}" = "nucleotide" ] && [ -n "${CDS_DIR}" ] && [ -d "${CDS_DIR}" ]; then
    log_step "3" "Back-translating to codon alignments with PAL2NAL"

    CODON_ALN_DIR="${WORKDIR}/codon_alignments"
    mkdir -p "${CODON_ALN_DIR}"

    # Build CDS lookup by orthogroup
    CDS_OG_DIR="${WORKDIR}/cds_by_orthogroup"
    mkdir -p "${CDS_OG_DIR}"

    python3 << PYEOF
import os
from Bio import SeqIO

cds_lookup = {}
cds_dir = "${CDS_DIR}"
for cds_file in os.listdir(cds_dir):
    if not cds_file.endswith(('.fa', '.fasta', '.fna')):
        continue
    for rec in SeqIO.parse(os.path.join(cds_dir, cds_file), "fasta"):
        cds_lookup[rec.id] = rec

sc_dir = "${SC_DIR}"
out_dir = "${CDS_OG_DIR}"
count = 0
for og_file in os.listdir(sc_dir):
    if not og_file.endswith('.fa'):
        continue
    og_name = og_file.replace('.fa', '')
    prot_recs = list(SeqIO.parse(os.path.join(sc_dir, og_file), "fasta"))
    cds_recs = []
    ok = True
    for rec in prot_recs:
        if rec.id in cds_lookup:
            cds_recs.append(cds_lookup[rec.id])
        else:
            ok = False
            break
    if ok and len(cds_recs) == len(prot_recs):
        SeqIO.write(cds_recs, os.path.join(out_dir, f"{og_name}_cds.fa"), "fasta")
        count += 1
print(f"Matched CDS for {count} orthogroups")
PYEOF

    for f in "${PROTEIN_ALN_DIR}"/*_aligned.fa; do
        name=$(basename "$f" _aligned.fa)
        cds_file="${CDS_OG_DIR}/${name}_cds.fa"
        outfile="${CODON_ALN_DIR}/${name}_codon.fa"

        if [ -f "${outfile}" ]; then
            continue
        fi

        if [ -f "${cds_file}" ]; then
            pal2nal.pl "$f" "${cds_file}" -output fasta > "${outfile}" 2>/dev/null || {
                log_warn "PAL2NAL failed for ${name}, using protein alignment"
                cp "$f" "${outfile}"
            }
        else
            log_warn "No CDS for ${name}, using protein alignment"
            cp "$f" "${outfile}"
        fi
    done

    FINAL_ALN_DIR="${CODON_ALN_DIR}"
    log_info "Codon alignments in: ${CODON_ALN_DIR}"
else
    log_step "3" "Skipping codon back-translation (no CDS provided)"
fi

#===============================================================================
# STEP 4: Trim alignments
#===============================================================================

log_step "4" "Trimming alignments with TrimAl"

TRIMMED_DIR="${WORKDIR}/trimmed_alignments"
mkdir -p "${TRIMMED_DIR}"

for f in "${FINAL_ALN_DIR}"/*.fa; do
    name=$(basename "$f" .fa)
    outfile="${TRIMMED_DIR}/${name}_trimmed.fa"
    if [ -f "${outfile}" ]; then
        continue
    fi
    trimal -in "$f" -out "${outfile}" -automated1 2>/dev/null || {
        log_warn "TrimAl failed for ${name}, copying untrimmed"
        cp "$f" "${outfile}"
    }
done

# Remove short alignments
REMOVED=0
for f in "${TRIMMED_DIR}"/*_trimmed.fa; do
    len=$(grep -v ">" "$f" | tr -d '\n ' | wc -c)
    if [ "${len}" -lt "${MIN_ALN_LENGTH}" ]; then
        rm "$f"
        REMOVED=$((REMOVED + 1))
    fi
done

N_TRIMMED=$(count_files "${TRIMMED_DIR}/*_trimmed.fa")
log_info "${N_TRIMMED} alignments passed trimming (${REMOVED} removed as too short)"

#===============================================================================
# STEP 5: Concatenate alignments
#===============================================================================

log_step "5" "Concatenating alignments"

CONCAT_DIR="${WORKDIR}/concatenated"
mkdir -p "${CONCAT_DIR}"

python3 << PYEOF
import os
from Bio import SeqIO
from collections import OrderedDict

trimmed_dir = "${TRIMMED_DIR}"
concat_file = "${CONCAT_DIR}/concatenated.fa"
partition_file = "${CONCAT_DIR}/partitions.txt"

# Collect all alignments
alignments = []
for f in sorted(os.listdir(trimmed_dir)):
    if not f.endswith('.fa'):
        continue
    records = OrderedDict()
    for rec in SeqIO.parse(os.path.join(trimmed_dir, f), "fasta"):
        # Use species name (first part before underscore or whole ID)
        records[rec.id] = str(rec.seq)
    if records:
        alignments.append((f, records))

if not alignments:
    print("ERROR: No alignments found!")
    exit(1)

# Get all species names (use first alignment to determine species set)
all_species = set()
for name, records in alignments:
    all_species.update(records.keys())
all_species = sorted(all_species)

# Concatenate
concat = {sp: "" for sp in all_species}
partitions = []
pos = 1

for aln_name, records in alignments:
    aln_len = len(next(iter(records.values())))
    partitions.append(f"{aln_name}\t{pos}-{pos + aln_len - 1}")
    
    for sp in all_species:
        if sp in records:
            concat[sp] += records[sp]
        else:
            # Fill missing species with gaps
            concat[sp] += "-" * aln_len
    pos += aln_len

total_len = len(next(iter(concat.values())))

# Write concatenated FASTA
with open(concat_file, 'w') as f:
    for sp, seq in concat.items():
        f.write(f">{sp}\n{seq}\n")

# Write partitions
with open(partition_file, 'w') as f:
    for p in partitions:
        f.write(p + "\n")

print(f"Concatenated: {len(all_species)} species, {total_len} sites, {len(alignments)} genes")
PYEOF

#===============================================================================
# STEP 6: Convert to PHYLIP format
#===============================================================================

log_step "6" "Converting to PHYLIP format for MCMCTree"

MCMC_DIR="${WORKDIR}/mcmctree"
mkdir -p "${MCMC_DIR}"

python3 << PYEOF
from Bio import SeqIO

records = list(SeqIO.parse("${CONCAT_DIR}/concatenated.fa", "fasta"))
n_taxa = len(records)
n_sites = len(records[0].seq)

with open("${MCMC_DIR}/alignment.phy", "w") as f:
    f.write(f"  {n_taxa}  {n_sites}\n")
    for rec in records:
        name = rec.id[:30].replace(" ", "_")
        f.write(f"{name:<31}{str(rec.seq)}\n")

print(f"PHYLIP file: {n_taxa} taxa, {n_sites} sites")
PYEOF

#===============================================================================
# STEP 7: Prepare the species tree (topology only, no branch lengths)
#===============================================================================

log_step "7" "Preparing species tree topology"

python3 << PYEOF
import re

# Read OrthoFinder species tree
with open("${SPECIES_TREE}") as f:
    tree = f.read().strip()

# Remove branch lengths
tree = re.sub(r':\d+\.?\d*([eE][+-]?\d+)?', '', tree)
# Remove internal node labels (N0, N1, etc. from OrthoFinder)
tree = re.sub(r'\)N\d+', ')', tree)

# Count species
species = re.findall(r'[A-Za-z_][A-Za-z0-9_.]*', tree.replace('(', '').replace(')', '').replace(';', ''))
species = [s for s in species if s]
n_species = len(species)

# Write tree file (without calibrations — user must add these manually)
with open("${MCMC_DIR}/species_tree_uncalibrated.nwk", 'w') as f:
    f.write(f"{n_species} 1\n")
    f.write(tree + "\n")

print(f"Species tree: {n_species} taxa")
print(f"Topology: {tree}")
print(f"\nWritten to: ${MCMC_DIR}/species_tree_uncalibrated.nwk")
print(f"\n{'='*60}")
print(f"IMPORTANT: You must now MANUALLY add calibrations!")
print(f"{'='*60}")
print(f"1. Open ${MCMC_DIR}/species_tree_uncalibrated.nwk")
print(f"2. Add calibration annotations after relevant ) characters")
print(f"3. Save as ${MCMC_DIR}/calibrated_tree.nwk")
print(f"\nExample:")
print(f"  Before: ((A,B),C);")
print(f"  After:  ((A,B)'B(0.05,0.10)',C)'B(0.55,1.00)';")
print(f"\nRemember: times are in units of 100 Myr!")
PYEOF

#===============================================================================
# STEP 8: MANUAL STEP — Add calibrations
#===============================================================================

log_step "8" "ADD CALIBRATIONS TO TREE (MANUAL STEP)"

CALIBRATED_TREE="${MCMC_DIR}/calibrated_tree.nwk"

if [ ! -f "${CALIBRATED_TREE}" ]; then
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────┐"
    echo "  │  ACTION REQUIRED: Add fossil calibrations to the tree  │"
    echo "  │                                                         │"
    echo "  │  1. Open: ${MCMC_DIR}/species_tree_uncalibrated.nwk    │"
    echo "  │  2. Add calibrations (see docs for syntax)             │"
    echo "  │  3. Save as: ${CALIBRATED_TREE}                        │"
    echo "  │  4. Re-run this script                                 │"
    echo "  └─────────────────────────────────────────────────────────┘"
    echo ""
    echo "  Calibration syntax examples:"
    echo "    Uniform:    'B(0.05, 0.10)'       — 5-10 Ma"
    echo "    Cauchy:     'L(0.05, 0.02, 0.1, 0.01)'"
    echo "    Gamma:      'G(2, 10)'"
    echo ""
    echo "  Time units: 100 Myr (so 65 Ma = 0.65)"
    echo ""
    exit 0
fi

log_info "Calibrated tree found: ${CALIBRATED_TREE}"

#===============================================================================
# STEP 9: MCMCTree — Compute Hessian (approximate likelihood)
#===============================================================================

log_step "9" "Computing Hessian matrix (usedata=3)"

HESSIAN_DIR="${MCMC_DIR}/hessian"
mkdir -p "${HESSIAN_DIR}"

# Determine seqtype value
if [ "${SEQ_TYPE}" = "nucleotide" ]; then
    SEQTYPE_VAL=0
else
    SEQTYPE_VAL=2
fi

# Write Hessian control file
cat > "${HESSIAN_DIR}/mcmctree.ctl" << CTLEOF
          seed = -1
       seqfile = ../alignment.phy
      treefile = ../calibrated_tree.nwk
      mcmcfile = mcmc_hessian.txt
       outfile = out_hessian.txt

         ndata = 1
       seqtype = ${SEQTYPE_VAL}
       usedata = 3
         clock = ${CLOCK_MODEL}
       RootAge = ${ROOT_AGE}

         model = ${SUBST_MODEL}
         alpha = ${GAMMA_ALPHA}
         ncatG = ${GAMMA_CATS}

     cleandata = 0

       BDparas = ${BD_LAMBDA} ${BD_MU} ${BD_RHO}
   rgene_gamma = ${RGENE_ALPHA} ${RGENE_BETA} 1
  sigma2_gamma = ${SIGMA2_ALPHA} ${SIGMA2_BETA} 1

      finetune = 1: .1 .1 .1 .1 .1 .1

         print = 1
        burnin = 2000
      sampfreq = 2
       nsample = 5000
CTLEOF

if [ -f "${MCMC_DIR}/in.BV" ]; then
    log_info "in.BV already exists, skipping Hessian computation"
else
    cd "${HESSIAN_DIR}"
    mcmctree mcmctree.ctl
    
    if [ -f "out.BV" ]; then
        cp out.BV "${MCMC_DIR}/in.BV"
        log_info "Hessian computed. in.BV ready."
    else
        log_error "Hessian computation failed — out.BV not produced"
        exit 1
    fi
    cd "${WORKDIR}"
fi

#===============================================================================
# STEP 10: MCMCTree — Prior-only analysis (DIAGNOSTIC)
#===============================================================================

log_step "10" "Running prior-only analysis (usedata=0)"

PRIOR_DIR="${MCMC_DIR}/prior_only"
mkdir -p "${PRIOR_DIR}"

cat > "${PRIOR_DIR}/mcmctree.ctl" << CTLEOF
          seed = -1
       seqfile = ../alignment.phy
      treefile = ../calibrated_tree.nwk
      mcmcfile = mcmc_prior.txt
       outfile = out_prior.txt

         ndata = 1
       seqtype = ${SEQTYPE_VAL}
       usedata = 0
         clock = ${CLOCK_MODEL}
       RootAge = ${ROOT_AGE}

         model = ${SUBST_MODEL}
         alpha = ${GAMMA_ALPHA}
         ncatG = ${GAMMA_CATS}

     cleandata = 0

       BDparas = ${BD_LAMBDA} ${BD_MU} ${BD_RHO}
   rgene_gamma = ${RGENE_ALPHA} ${RGENE_BETA} 1
  sigma2_gamma = ${SIGMA2_ALPHA} ${SIGMA2_BETA} 1

      finetune = 1: .1 .1 .1 .1 .1 .1

         print = 1
        burnin = 50000
      sampfreq = 50
       nsample = 100000
CTLEOF

cd "${PRIOR_DIR}"
mcmctree mcmctree.ctl
cd "${WORKDIR}"

log_info "Prior-only analysis complete. Check out_prior.txt for effective priors."
log_info "Compare effective priors to your intended calibrations!"

#===============================================================================
# STEP 11: MCMCTree — Run MCMC chains
#===============================================================================

log_step "11" "Running ${N_CHAINS} independent MCMC chains"

for CHAIN in $(seq 1 ${N_CHAINS}); do
    CHAIN_DIR="${MCMC_DIR}/run${CHAIN}"
    mkdir -p "${CHAIN_DIR}"
    cp "${MCMC_DIR}/in.BV" "${CHAIN_DIR}/"

    SEED=$((CHAIN * -1))

    cat > "${CHAIN_DIR}/mcmctree.ctl" << CTLEOF
          seed = ${SEED}
       seqfile = ../alignment.phy
      treefile = ../calibrated_tree.nwk
      mcmcfile = mcmc_samples.txt
       outfile = out_mcmctree.txt

         ndata = 1
       seqtype = ${SEQTYPE_VAL}
       usedata = 2
         clock = ${CLOCK_MODEL}
       RootAge = ${ROOT_AGE}

         model = ${SUBST_MODEL}
         alpha = ${GAMMA_ALPHA}
         ncatG = ${GAMMA_CATS}

     cleandata = 0

       BDparas = ${BD_LAMBDA} ${BD_MU} ${BD_RHO}
   rgene_gamma = ${RGENE_ALPHA} ${RGENE_BETA} 1
  sigma2_gamma = ${SIGMA2_ALPHA} ${SIGMA2_BETA} 1

      finetune = 1: .1 .1 .1 .1 .1 .1

         print = -1
        burnin = ${BURNIN}
      sampfreq = ${SAMPFREQ}
       nsample = ${NSAMPLE}
CTLEOF

    log_info "Starting chain ${CHAIN} (seed=${SEED})..."
    cd "${CHAIN_DIR}"
    mcmctree mcmctree.ctl &
    cd "${WORKDIR}"
done

# Wait for all chains to finish
wait
log_info "All ${N_CHAINS} chains completed."

#===============================================================================
# STEP 12: Check convergence
#===============================================================================

log_step "12" "Checking convergence between chains"

python3 << PYEOF
import numpy as np
import os

mcmc_dir = "${MCMC_DIR}"
n_chains = ${N_CHAINS}

chains = []
for i in range(1, n_chains + 1):
    fpath = os.path.join(mcmc_dir, f"run{i}", "mcmc_samples.txt")
    if os.path.exists(fpath):
        data = np.genfromtxt(fpath, names=True, deletechars='')
        chains.append(data)
        print(f"Chain {i}: {len(data)} samples loaded")
    else:
        print(f"WARNING: Chain {i} output not found at {fpath}")

if len(chains) < 2:
    print("Cannot check convergence with fewer than 2 chains")
    exit(0)

# Compare posterior means for node ages (columns starting with 't_')
print(f"\n{'Parameter':<15} {'Chain 1 Mean':>14} {'Chain 2 Mean':>14} {'Diff %':>10}")
print("-" * 55)

all_ok = True
for col in chains[0].dtype.names:
    if col.startswith('t_n') or col == 'mu' or col == 'sigma2':
        m1 = np.mean(chains[0][col])
        m2 = np.mean(chains[1][col])
        avg = (m1 + m2) / 2
        if avg != 0:
            diff_pct = abs(m1 - m2) / avg * 100
        else:
            diff_pct = 0
        
        flag = "" if diff_pct < 5 else " ⚠"
        if diff_pct >= 5:
            all_ok = False
        print(f"{col:<15} {m1:>14.4f} {m2:>14.4f} {diff_pct:>9.1f}%{flag}")

print()
if all_ok:
    print("✓ All parameters agree within 5% between chains.")
else:
    print("⚠ Some parameters differ by >5%. Consider running longer chains.")

# Simple ESS calculation
def compute_ess(x):
    n = len(x)
    if n < 10:
        return n
    mean = np.mean(x)
    var = np.var(x)
    if var == 0:
        return n
    sum_acf = 0
    for lag in range(1, min(n // 2, 200)):
        acf = np.mean((x[:-lag] - mean) * (x[lag:] - mean)) / var
        if acf < 0.05:
            break
        sum_acf += acf
    return max(1, n / (1 + 2 * sum_acf))

print(f"\n{'Parameter':<15} {'ESS Chain 1':>14} {'ESS Chain 2':>14}")
print("-" * 45)
low_ess = False
for col in chains[0].dtype.names:
    if col.startswith('t_n') or col == 'mu' or col == 'sigma2':
        ess1 = compute_ess(chains[0][col])
        ess2 = compute_ess(chains[1][col])
        flag = "" if (ess1 > 200 and ess2 > 200) else " ⚠"
        if ess1 < 200 or ess2 < 200:
            low_ess = True
        print(f"{col:<15} {ess1:>14.0f} {ess2:>14.0f}{flag}")

print()
if low_ess:
    print("⚠ Some ESS values < 200. Increase burnin and/or nsample.")
else:
    print("✓ All ESS values > 200. Good convergence.")
PYEOF

#===============================================================================
# STEP 13: Extract results
#===============================================================================

log_step "13" "Extracting results"

# Use chain 1 results as the primary output
RESULTS_FILE="${MCMC_DIR}/run1/out_mcmctree.txt"
FIGTREE_FILE="${MCMC_DIR}/run1/FigTree.tre"

if [ -f "${RESULTS_FILE}" ]; then
    echo ""
    echo "  ┌──────────────────────────────────────────┐"
    echo "  │        TIME-CALIBRATED PHYLOGENY          │"
    echo "  │           RESULTS SUMMARY                 │"
    echo "  └──────────────────────────────────────────┘"
    echo ""
    
    # Extract node age lines
    grep -E "^[[:space:]]*(t_n|Node)" "${RESULTS_FILE}" | head -30 || true
    
    echo ""
    echo "Full results: ${RESULTS_FILE}"
    
    if [ -f "${FIGTREE_FILE}" ]; then
        echo "FigTree tree:  ${FIGTREE_FILE}"
        echo ""
        echo "Open with: figtree ${FIGTREE_FILE}"
    fi
else
    log_error "Results file not found: ${RESULTS_FILE}"
fi

# Copy final outputs to a results directory
FINAL_DIR="${WORKDIR}/RESULTS"
mkdir -p "${FINAL_DIR}"
cp "${MCMC_DIR}/run1/out_mcmctree.txt" "${FINAL_DIR}/" 2>/dev/null || true
cp "${MCMC_DIR}/run1/FigTree.tre" "${FINAL_DIR}/" 2>/dev/null || true
cp "${MCMC_DIR}/run1/mcmc_samples.txt" "${FINAL_DIR}/" 2>/dev/null || true
cp "${MCMC_DIR}/calibrated_tree.nwk" "${FINAL_DIR}/" 2>/dev/null || true

#===============================================================================
# DONE
#===============================================================================

echo ""
echo "================================================================"
echo "  PIPELINE COMPLETE"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================"
echo ""
echo "  Output directory: ${FINAL_DIR}/"
echo ""
echo "  Files:"
echo "    out_mcmctree.txt   — Node ages with 95% HPD intervals"
echo "    FigTree.tre        — Dated tree for FigTree visualization"
echo "    mcmc_samples.txt   — Raw MCMC samples (load in Tracer)"
echo "    calibrated_tree.nwk — Your calibrated input tree"
echo ""
echo "  Next steps:"
echo "    1. Check convergence in Tracer (ESS > 200 for all params)"
echo "    2. Compare prior-only vs posterior (${MCMC_DIR}/prior_only/)"
echo "    3. Run sensitivity analyses (different clock, calibrations)"
echo "    4. Visualize with FigTree or ggtree in R"
echo ""
