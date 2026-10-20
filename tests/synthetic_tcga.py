"""Generate a small synthetic TCGA-BRCA snapshot with the layout analyze_tcga.py reads.

The data are random with a few planted signals used by the tests:

* ``PHGDH`` expression is higher in clinical TNBC than in receptor-positive tumours;
* ``WT1`` expression increases the PFI hazard (HR about 1.6 per SD);
* ``LAMTOR2`` is present only under its previous symbol ``ROBLD3`` in RNA and CNA;
* ``LINC03007`` is absent from every platform;
* ``DUPGENE`` (not a candidate) is duplicated in the RNA matrix.

Usage as a script::

    python tests/synthetic_tcga.py OUTPUT_DIR [--n 320] [--seed 0]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

KO = ["LAMTOR2", "WDR24", "DERL2", "NME6", "TP53BP1", "PARG"]
ACT = ["TP53I11", "FOXQ1", "CERK", "PHF19", "SLC35D1", "PHLDA2", "WT1", "HEPHL1",
       "LINC03007", "SLC47A1", "ELL2", "INHBA", "PHGDH", "ABCB1"]
GENES = KO + ACT
PATHS = ["MTORC1_SIGNALING", "UNFOLDED_PROTEIN_RESPONSE", "OXIDATIVE_PHOSPHORYLATION",
         "GLYCOLYSIS", "DNA_REPAIR", "G2M_CHECKPOINT", "TGF_BETA_SIGNALING",
         "EPITHELIAL_MESENCHYMAL_TRANSITION", "APOPTOSIS", "HYPOXIA", "E2F_TARGETS"]
RPPA_READOUTS = ["RPS6KB1|p70S6K_pT389", "RPS6|S6_pS235_S236", "RPS6|S6_pS240_S244",
                 "EIF4EBP1|4E-BP1_pT37T46", "MAPK1 MAPK3|MAPK_pT202_Y204",
                 "AKT1 AKT2 AKT3|Akt_pS473", "TP53BP1|53BP1", "RAD51|Rad51", "XBP1|XBP1"]
ENTREZ = {g: 100000 + i for i, g in enumerate(GENES)}
ENTREZ["LAMTOR2"] = 28956
SUBTYPES = ["BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Basal", "BRCA_Normal", None]
STAGES = ["Stage I", "Stage IA", "Stage IB", "Stage IIA", "Stage IIB", "Stage IIIA",
          "Stage IIIB", "Stage IIIC", "Stage IV", "Stage X", "[Discrepancy]", None]
CBIO_HEADER = "#{0}\n#{0}\n#STRING\n#1\n"


def hallmark_sets(rng: np.random.Generator, genes_per_set: int = 24) -> dict[str, list[str]]:
    """Synthetic Hallmark-like sets, with some candidate genes as members."""
    sets = {}
    for p in PATHS:
        members = [f"{p[:6]}_G{i:02d}" for i in range(genes_per_set)]
        if p == "DNA_REPAIR":
            members += ["TP53BP1", "PARG"]
        if p == "GLYCOLYSIS":
            members += ["PHGDH"]
        if p == "MTORC1_SIGNALING":
            members += ["WDR24"]
        sets[p] = members
    shared = [f"SHARED_G{i:02d}" for i in range(6)]
    sets["E2F_TARGETS"] += shared
    sets["G2M_CHECKPOINT"] += shared
    return sets


def _write_cbio(path: Path, frame: pd.DataFrame, header: str) -> None:
    path.write_text(CBIO_HEADER.format(header) + frame.to_csv(sep="\t", index=False))


def make_dataset(raw: Path, n: int = 320, seed: int = 0,
                 duplicate_candidate: str | None = None) -> dict:
    """Write every input file into ``raw`` and return the planted truth."""
    raw = Path(raw)
    raw.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    patients = [f"TCGA-AA-{i:04d}" for i in range(1, n + 1)]
    samples = [p + "-01" for p in patients]

    # Receptor status: roughly one fifth clinical TNBC; a few indeterminate/missing.
    tnbc = rng.random(n) < 0.22
    er = np.where(tnbc, "Negative", rng.choice(["Positive", "Negative"], n, p=[.85, .15]))
    pr = np.where(tnbc, "Negative", rng.choice(["Positive", "Negative"], n, p=[.7, .3]))
    ihc = np.where(tnbc, "Negative", rng.choice(["Positive", "Negative", "Equivocal"], n, p=[.2, .6, .2]))
    fish = np.where(tnbc, rng.choice(["Negative", None], n, p=[.6, .4]),
                    rng.choice(["Positive", "Negative", None], n, p=[.15, .45, .4]))
    noise = rng.random(n) < 0.04
    er = np.where(noise, rng.choice(["Indeterminate", None], n), er)
    firehose = pd.DataFrame({"PATIENT_ID": patients, "ER_STATUS_BY_IHC": er,
                             "PR_STATUS_BY_IHC": pr, "IHC_HER2": ihc, "HER2_FISH_STATUS": fish})

    subtype = np.where(tnbc, rng.choice(["BRCA_Basal", "BRCA_Normal", None], n, p=[.75, .15, .1]),
                       rng.choice(SUBTYPES, n, p=[.45, .2, .12, .05, .08, .1]))
    patient = pd.DataFrame({"PATIENT_ID": patients, "SUBTYPE": subtype,
                            "AGE": rng.integers(28, 90, n)})

    sample = pd.DataFrame({"SAMPLE_ID": samples, "PATIENT_ID": patients, "SAMPLE_TYPE": "Primary"})
    extra = pd.DataFrame({"SAMPLE_ID": [patients[0] + "-06", patients[1] + "-01B"],
                          "PATIENT_ID": patients[:2], "SAMPLE_TYPE": ["Metastatic", "Primary"]})
    sample = pd.concat([sample, extra], ignore_index=True)

    # Expression: log-normal RSEM with planted contrast (PHGDH in TNBC).
    sets = hallmark_sets(rng)
    set_genes = sorted({g for members in sets.values() for g in members} - set(GENES))
    rna_genes = [g for g in GENES if g != "LINC03007"] + set_genes + ["EXTRA_G01", "EXTRA_G02", "DUPGENE"]
    log_expr = pd.DataFrame(rng.normal(8, 1.5, (len(rna_genes), n)), index=rna_genes, columns=samples)
    log_expr.loc["PHGDH", tnbc] += 1.2
    log_expr.loc[set_genes[:3]] = 0.0  # constant genes give an undefined z-score
    rsem = np.round(2 ** log_expr - 1, 4).clip(lower=0)
    rna = rsem.reset_index().rename(columns={"index": "Hugo_Symbol"})
    rna.insert(1, "Entrez_Gene_Id", [ENTREZ.get(g, 200000 + i) for i, g in enumerate(rna_genes)])
    rna.loc[rna.Hugo_Symbol.eq("LAMTOR2"), "Hugo_Symbol"] = "ROBLD3"
    dup = rna.loc[rna.Hugo_Symbol.eq("DUPGENE")].copy()
    dup["Entrez_Gene_Id"] = 999999
    missing_symbol = rna.iloc[[-1]].copy()
    missing_symbol["Hugo_Symbol"] = np.nan
    missing_symbol["Entrez_Gene_Id"] = 999998
    rna = pd.concat([rna, dup, missing_symbol], ignore_index=True)
    if duplicate_candidate:
        twin = rna.loc[rna.Hugo_Symbol.eq(duplicate_candidate)].copy()
        twin["Entrez_Gene_Id"] = 999997
        twin.iloc[:, 2:] = 0.0
        rna = pd.concat([rna, twin], ignore_index=True)
    rna.to_csv(raw / "data_mrna_seq_v2_rsem.txt", sep="\t", index=False)

    (raw / "h.all.v2025.1.Hs.symbols.gmt").write_text("".join(
        f"HALLMARK_{p}\thttps://example.org/{p}\t" + "\t".join(members) + "\n"
        for p, members in sets.items()))

    # Copy number: GISTIC calls; LAMTOR2 reachable only through Entrez 28956 -> ROBLD3.
    cna_genes = [g for g in GENES if g != "LINC03007"] + ["EXTRA_G01"]
    calls = rng.choice([-2, -1, 0, 1, 2], (len(cna_genes), n), p=[.03, .22, .5, .2, .05]).astype(float)
    calls[:, rng.random(n) < 0.05] = np.nan
    cna = pd.DataFrame(calls, index=cna_genes, columns=samples).reset_index().rename(columns={"index": "Hugo_Symbol"})
    cna.insert(1, "Entrez_Gene_Id", [ENTREZ.get(g, 300000) for g in cna_genes])
    cna.loc[cna.Hugo_Symbol.eq("LAMTOR2"), "Hugo_Symbol"] = "ROBLD3"
    cna.to_csv(raw / "data_cna.txt", sep="\t", index=False)

    # RPPA for a subset of samples.
    rppa_samples = [s for s in samples if rng.random() < 0.7]
    rppa = pd.DataFrame(rng.normal(0, 1, (len(RPPA_READOUTS) + 1, len(rppa_samples))),
                        index=RPPA_READOUTS + ["EXTRA|extra"], columns=rppa_samples)
    rppa.index.name = "Composite.Element.REF"
    rppa.reset_index().to_csv(raw / "data_rppa.txt", sep="\t", index=False)

    # Survival from TCGA-CDR with a planted WT1 effect on PFI.
    wt1_z = (log_expr.loc["WT1"] - log_expr.loc["WT1"].mean()) / log_expr.loc["WT1"].std(ddof=1)
    hazard = 0.0006 * np.exp(0.47 * wt1_z.to_numpy())
    pfi_time = rng.exponential(1 / hazard)
    censor = rng.uniform(200, 3000, n)
    pfi = (pfi_time <= censor).astype(float)
    pfi_time = np.minimum(pfi_time, censor)
    os_time = pfi_time + rng.exponential(400, n)
    os_event = ((pfi == 1) & (rng.random(n) < .6)).astype(float)
    dfi = pfi.copy()
    dfi_time = pfi_time.copy()
    dfi[rng.random(n) < .3] = np.nan
    cdr = pd.DataFrame({"bcr_patient_barcode": patients, "type": "BRCA",
                        "age_at_initial_pathologic_diagnosis": patient.AGE.astype(float),
                        "ajcc_pathologic_tumor_stage": rng.choice(STAGES, n, p=[.08, .1, .04, .25, .2, .12, .07, .04, .03, .02, .02, .03]),
                        "PFI": pfi, "PFI.time": np.round(pfi_time), "DFI": dfi, "DFI.time": np.round(dfi_time),
                        "OS": os_event, "OS.time": np.round(os_time), "Redaction": None})
    cdr.loc[rng.random(n) < .03, "age_at_initial_pathologic_diagnosis"] = np.nan
    cdr.loc[rng.random(n) < .02, "PFI.time"] = 0.0
    other = cdr.iloc[:5].copy()
    other["type"] = "LUAD"
    other["bcr_patient_barcode"] = ["TCGA-ZZ-%04d" % i for i in range(5)]
    cdr = pd.concat([cdr, other], ignore_index=True)
    with pd.ExcelWriter(raw / "TCGA-CDR-SupplementalTableS1.xlsx") as xl:
        pd.DataFrame({"note": ["synthetic"]}).to_excel(xl, sheet_name="ReadMe", index=False)
        cdr.to_excel(xl, sheet_name="TCGA-CDR", index=False)

    purity = pd.DataFrame({"Sample": [s + "A" for s in samples], "Purity": np.round(rng.uniform(.2, .95, n), 2)})
    purity.loc[rng.random(n) < .06, "Purity"] = np.nan
    purity = purity.iloc[:-7]  # a few samples without purity
    purity.to_csv(raw / "purity.tsv", sep="\t", index=False)

    _write_cbio(raw / "data_clinical_sample.txt", sample, "Sample Identifier")
    _write_cbio(raw / "data_clinical_patient.txt", patient, "Patient Identifier")
    _write_cbio(raw / "data_clinical_firehose_patient.txt", firehose, "Patient Identifier")

    return {"n": n, "tnbc": pd.Series(tnbc, index=samples), "sets": sets,
            "constant_genes": set_genes[:3], "rppa_samples": rppa_samples}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("output", type=Path)
    parser.add_argument("--n", type=int, default=320)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    truth = make_dataset(args.output, args.n, args.seed)
    print(f"wrote synthetic snapshot for {truth['n']} patients to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
