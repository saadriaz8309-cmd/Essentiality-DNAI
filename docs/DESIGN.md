# DNAI — design (locked for v1)

## Goal

A literature-unbiased gene-essentiality **discovery** pipeline for *Drosophila melanogaster*.
It is not a leaderboard classifier. The output is a ranked list of **candidate essential genes**:
genes the model confidently predicts essential but that FlyBase marks viable or has never annotated.

## Why it is novel

Existing fly-essentiality models use features that track how well studied a gene is
(GO terms, text-mined interactions, lethal-allele counts), and *Drosophila* is very
annotation-biased. DNAI excludes every study-biased feature and uses only
annotation-independent signals. It treats research intensity as a **confounder** to correct for.
The candidate list is the essential genes that the biased models systematically miss.

Constraints: CPU-only, local, cheap. No foundation models, GPU, or API. All-fly (no human DepMap in v1).

## Labels

- Source: FlyBase `genotype_phenotype_data` (lethal vs viable phenotypes), joined to genes via `fbal_to_fbgn`.
- Essentiality score: `ES = lethal / total` over a gene's lethal + viable phenotype annotations.
  - The original handoff wrote `ES = lethal² / total²`. That is ES squared, which is monotone in
    ES, so the ~0.9 / ~0.1 thresholds would be ~0.95 / ~0.32 on the plain ratio. To be confirmed before Stage 2.
- **Confident positive**: ES > ~0.9 with ≥ N annotations.
- **Confident negative**: ES < ~0.1 with explicit viable evidence.
- **Unlabeled**: everything else, including all unannotated genes.
- **Critical rule:** a gene with no lethal annotation is *unlabeled*, never negative.
- Consider restricting to single-allele genotypes for a cleaner signal.

## Features (annotation-independent only)

| Feature | Source |
|---|---|
| protein length (longest isoform) | `dmel-all-translation.fasta` |
| gene length, exon count, GC | `dmel-all-r6.69.gff.gz` (pending) |
| mean expression, expression breadth / entropy | `gene_rpkm_report` |
| paralog count | `dmel_paralogs` |
| network degree | STRING v12 physical links, **experimental + database channels only** (drop textmining + coexpression) (pending) |

## Model

LightGBM with bagging-PU: K learners, each trained on all positives plus a random unlabeled
subsample used as pseudo-negatives. Scores are averaged across learners.

## CV / evaluation

- Group-aware CV, grouped by **paralog family** (prevents near-duplicate leakage).
- **AUPRC** is the primary metric.
- Baselines include a deliberately biased model (with GO/annotation features). The AUPRC gap is
  the study-bias effect, which is a finding, not a failure.
- External validation: eLife Viswanatha et al. (article 36333) ~1,235 cell-fitness genes, held out from training.
- Research-intensity covariate: per-gene reference/allele count. It is **kept out of the features**
  and used only to stratify evaluation (v1) or regress out (v2).

## Build order

1. Stage 1: gene universe (canonical protein-coding FBgn table + FBpp map)
2. Stage 2: labels (PU setup)
3. Stage 3: feature matrix (network feature added when STRING lands)
4. Stage 5: PU LightGBM + group-aware CV
5. Stage 6: evaluation (AUPRC, biased baseline, external validation)
6. Stage 7: candidate list + Streamlit app

## Raw file parsing notes

All files are FB2026_03 / genome r6.69. Files go in `data/raw/` under the clean names below
(`.tsv`, `.tsv.gz`, `.txt`, `.gz` and no extension are all accepted).

| Clean name | Role | Layout |
|---|---|---|
| `fbgn_annotation_ID` | ID backbone (17,872 genes) | no header, `#`/`##` comments; cols: symbol, organism (Dmel), primary_FBgn, secondary_FBgns (comma-sep), annotation_ID (CG), secondary_annotation_IDs |
| `gene_map_table` | FBgn + cytological location | `#` comments; cols: organism, symbol, primary_FBid, recombination_loc, cytogenetic_loc, sequence_loc |
| `fbgn_fbtr_fbpp_expanded` | FBgn↔FBpp + gene_type | cols: organism, gene_type, FBgn, gene_symbol, gene_fullname, annotation_ID, transcript_type, FBtr, transcript_symbol, FBpp, polypeptide_symbol. Filter `gene_type == "protein_coding_gene"` |
| `gene_rpkm_report` | expression (2.19M rows) | cols: Release, FBgn, GeneSymbol, Parent_FBlc, Parent_name, RNASource_FBlc, RNASource_name, RPKM (index 7), ...; one row per gene × sample |
| `fbal_to_fbgn` | allele → gene (308,912) | cols: AlleleID (FBal), AlleleSymbol, GeneID (FBgn), GeneSymbol |
| `genotype_phenotype_data` | labels (400,940 rows; ~59.7k lethal / ~48k viable) | header: genotype_symbols, genotype_FBids, phenotype_name, phenotype_id, qualifier_names, qualifier_ids, reference. Parse phenotype_name for "lethal"/"viable"; extract `FBal\d+` tokens from genotype_FBids |
| `dmel_paralogs` | paralog pairs (223,056) | col0 = FBgn, col5 = paralog FBgn; paralog count = distinct partners per FBgn |
| `dmel-all-translation.fasta` | protein length (30,837 proteins) | FBpp headers with `parent=FBgn...,FBtr...` and `length=NNN`; take the longest isoform per gene |

Still to download: `dmel-all-r6.69.gff.gz`, `7227.protein.physical.links.full.v12.0.txt.gz`,
eLife 36333 Supplementary file 1.
