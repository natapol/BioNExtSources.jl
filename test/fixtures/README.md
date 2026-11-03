# Test fixtures

Small extracts of real database exports, downloaded on 2026-11-03, used to test the readers
against the actual file formats:

| File | Source |
|---|---|
| `kegg_hsa04115.kgml`, `kegg_hsa00010.kgml` | KEGG REST `get/<map>/kgml` (p53 signalling, glycolysis) |
| `kegg_conv_uniprot.txt` | KEGG REST `conv/uniprot/hsa:7157+hsa:4193+hsa:1026` |
| `uniprot_p53_mdm2.tsv` | `fetch_uniprot("accession:P04637 OR accession:Q00987")` (all `UNIPROT_FIELDS`) |
| `intact_p53.mitab27.txt` | IntAct PSICQUIC, query `P04637`, MITAB 2.7, first 20 rows |
| `biogrid_sample.tab3.txt`, `biogrid_sample.chemtab.txt` | BioGRID 5.0.263, Fanconi anemia project, first rows |
| `reactome_R-HSA-69541.owl.gz` | Reactome BioPAX Level 3 export of R-HSA-69541 (Stabilization of p53) |

`signor_handmade.tsv` and `signor_handmade_complexes.csv` are **hand-written** in SIGNOR's
documented export layout (SIGNOR's server could not be reached when the tests were written);
their SIGNOR ids (`SIGNOR-TEST-*`) are not real. Replace them with a real export when possible.

KEGG content is free for academic use only (https://www.kegg.jp/kegg/legal.html); the other
databases publish under open licences. Check each source's terms before redistributing data.
