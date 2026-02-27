import pandas as pd
import os

# ── Paths ──────────────────────────────────────────────────────────────────
STAR_DIR   = os.path.expanduser("~/rnaseq-pipeline/results/star")
COUNTS_IN  = os.path.expanduser("~/rnaseq-pipeline/results/counts/counts_clean.txt")
DESEQ2_IN  = os.path.expanduser("~/rnaseq-pipeline/results/deseq2/DE_results.csv")
SAMPLES_IN = os.path.expanduser("~/rnaseq-pipeline/config/samples.csv")
OUT_DIR    = os.path.expanduser("~/rnaseq-pipeline/shiny_inputs")
os.makedirs(OUT_DIR, exist_ok=True)

samples = ["EV_15C_rep1","EV_15C_rep2","EV_15C_rep3",
           "EV_20C_rep1","EV_20C_rep2","EV_20C_rep3"]

# ── File 1: Counts Matrix ──────────────────────────────────────────────────
counts = pd.read_csv(COUNTS_IN, sep="\t")
counts.rename(columns={"Geneid": "gene_id"}, inplace=True)
counts.to_csv(f"{OUT_DIR}/counts_matrix.csv", index=False)
print(f"[1] counts_matrix.csv  — {counts.shape[0]} genes x {counts.shape[1]-1} samples")

# ── File 2: Sample Metadata ────────────────────────────────────────────────
meta = pd.read_csv(SAMPLES_IN)
meta.rename(columns={"sample": "sample_id"}, inplace=True)
meta.to_csv(f"{OUT_DIR}/sample_metadata.csv", index=False)
print(f"[2] sample_metadata.csv — {len(meta)} samples")

# ── File 3: QC Metrics from STAR logs + samtools markdup metrics ───────────
def parse_star_log(path):
    metrics = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "Number of input reads" in line:
                metrics["input_reads"] = int(line.split("|")[1].strip())
            elif "Uniquely mapped reads %" in line:
                metrics["pct_mapped"] = float(line.split("|")[1].strip().replace("%",""))
    return metrics

def parse_markdup_metrics(path):
    """Parse samtools markdup metrics file for true PCR duplication rate.
    duplication_rate = DUPLICATE PRIMARY TOTAL / EXAMINED
    """
    examined = None
    dup_total = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("EXAMINED:"):
                examined = int(line.split(":")[1].strip())
            elif line.startswith("DUPLICATE PRIMARY TOTAL:"):
                dup_total = int(line.split(":")[1].strip())
    if examined and examined > 0:
        return round(dup_total / examined, 4)
    return None

qc_rows = []
for s in samples:
    star_log    = os.path.join(STAR_DIR, s, f"{s}.Log.final.out")
    markdup_log = os.path.join(STAR_DIR, s, f"{s}.markdup.metrics.txt")
    star     = parse_star_log(star_log)
    dup_rate = parse_markdup_metrics(markdup_log)
    qc_rows.append({
        "sample_id":             s,
        "condition":             "EV_15C" if "15C" in s else "EV_20C",
        "library_size_millions": round(star["input_reads"] / 1e6, 2),
        "pct_mapped":            star["pct_mapped"],
        "duplication_rate":      dup_rate,
    })

qc = pd.DataFrame(qc_rows)
qc.to_csv(f"{OUT_DIR}/qc_metrics.csv", index=False)
print(f"[3] qc_metrics.csv — {len(qc)} samples")
print(qc.to_string(index=False))

# ── File 4: DESeq2 Results ─────────────────────────────────────────────────
de = pd.read_csv(DESEQ2_IN)
de.rename(columns={"gene": "gene_id"}, inplace=True)
de.to_csv(f"{OUT_DIR}/DE_results.csv", index=False)
print(f"\n[4] DE_results.csv — {len(de)} genes, columns: {list(de.columns)}")

print(f"\nAll files written to {OUT_DIR}/")
print("  counts_matrix.csv")
print("  sample_metadata.csv")
print("  qc_metrics.csv")
print("  DE_results.csv")
