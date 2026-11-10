# BioNExtSources.jl

[![CI](https://github.com/natapol/BioNExtSources.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/natapol/BioNExtSources.jl/actions/workflows/CI.yml)

Populate [BioNExt](https://github.com/natapol/BioNExt.jl) knowledge graphs from the major biological databases.
Every reader returns an ordinary `KnowledgeGraph` that follows the BioNExt ontology and
passes `validate`, and `populate!` merges graphs from several databases into one network,
recognising the same molecule under different identifiers.

```julia
using BioNExt, BioNExtSources

kg = KnowledgeGraph()
populate!(kg,
    read_uniprot(fetch_uniprot("organism_id:9606 AND reviewed:true AND gene:TP53")),
    read_intact(fetch_intact("P04637 AND species:9606")),
    read_kgml(fetch_kgml("hsa04115"); uniprot = fetch_kegg_uniprot("hsa"), taxon = 9606),
    read_reactome(fetch_reactome_biopax("R-HSA-69541")),
    read_signor(fetch_signor()),
    read_biogrid("BIOGRID-ORGANISM-Homo_sapiens-5.0.263.tab3.txt"))
validate(kg)          # -> Issue[]
```

`examples/p53_network.jl` builds a complete network this way.

## Installation

Julia 1.10 or newer. BioNExt is registered in Julia's General registry; BioNExtSources isn't
registered yet, so add it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/natapol/BioNExtSources.jl")   # also installs BioNExt
```

To work on both packages, clone them side by side and develop BioNExt into this one:

```sh
git clone https://github.com/natapol/BioNExt.jl
git clone https://github.com/natapol/BioNExtSources.jl
cd BioNExtSources.jl
julia --project=. -e 'using Pkg; Pkg.develop(path = "../BioNExt.jl")'
```

## Databases

| Database | Read | Download | Format | Tested against |
|---|---|---|---|---|
| **UniProt** | `read_uniprot` | `fetch_uniprot(query)` | TSV | live API and stored export |
| **KEGG** | `read_kgml`, `read_kegg_conv` | `fetch_kgml(map)`, `fetch_kegg_uniprot(org)`, `fetch_kegg(op)` | KGML, REST text | live API; p53-signalling and glycolysis maps |
| **IntAct** and any PSICQUIC service (MINT, BioGRID, InnateDB, …) | `read_intact`, `read_mitab` | `fetch_intact(query)`, `fetch_psicquic(service, query)` | PSI-MITAB 2.5–2.8 | live API and stored export |
| **BioGRID** | `read_biogrid`, `read_biogrid_chemicals`; also `read_mitab` | `fetch_biogrid(genes; accesskey)` | TAB3, chemtab, `.zip` releases | release 5.0.263 files (156,669 interactions); the REST download needs an access key and was not tested |
| **SIGNOR** | `read_signor`, `read_signor_complexes`; causalTab via `read_mitab` | `fetch_signor(; organism)` | TSV, MITAB 2.8 | **hand-written file in SIGNOR's documented layout only**: SIGNOR's server was unreachable during development |
| **Reactome** | `read_reactome` | `fetch_reactome_biopax(id)` | BioPAX Level 3 | live API and stored export |
| **BioCyc** (EcoCyc, YeastCyc, HumanCyc, …) | `read_biocyc` | none: download manually (see below) | BioPAX Level 3 | through the shared BioPAX reader; no BioCyc file tested |
| **Pathway Commons** (Reactome, KEGG, PANTHER, HumanCyc, …) | `read_biopax` | `fetch_pathwaycommons(uri)` | BioPAX Level 3 | live API |

Downloads are cached in `~/.cache/bionext`, or in `$BIONEXT_CACHE`, and reused. Pass
`refresh = true` to download again, or `cache = dir` to use another folder. Transient server
errors are retried, after 2, 4 and 8 s. Every reader also accepts a local file, `.gz` and `.zip`
files, or an `IO`.

### Before you start: access and licences

- **KEGG** data are free for **academic use only**; other uses need a licence
  ([terms](https://www.kegg.jp/kegg/legal.html)).
- **BioCyc** requires a (free or paid) account and blocks scripted downloads. Export a pathway
  or a whole database from the BioCyc website as *BioPAX Level 3*, then call `read_biocyc(file)`.
- **BioGRID**'s REST service needs a free access key (https://webservice.thebiogrid.org). For a
  whole organism, it's simpler to download a release file from https://downloads.thebiogrid.org,
  such as `BIOGRID-ORGANISM-…tab3.zip`. `read_biogrid` reads the zip directly.
- UniProt, IntAct, Reactome and Pathway Commons are openly licensed. SIGNOR and BioGRID
  publish under their own open terms. Check each source's licence before redistributing a
  graph built from it. Every resource records its source in `provenance`.

## How data become BioNExt

### Identifiers

The same molecule must get the same id in every source, or the sources can't be merged:

| Kind | Id |
|---|---|
| proteins | `UniProt:<accession>`; isoforms merge into the canonical accession unless `isoforms = :keep`, which keeps `P04637-2` with `isoform_of` |
| small molecules | `CHEBI:<id>` when the source gives a ChEBI id, otherwise the source's id (`kegg.compound:C00002`, `DrugBank:DB00126`, `pubchem.compound:…`) |
| genes | `NCBIGene:<id>` (UniProt, BioGRID), `kegg.genes:hsa:7157` (KEGG) |
| gene products without a known UniProt accession | `gene_product/<gene id>`, with the gene id in `db_xrefs` |
| complexes, reactions, pathways, edges | the source's stable id: `Reactome:R-HSA-…`, `IntAct:EBI-…`, `SIGNOR:SIGNOR-…`, `kegg.reaction:R00710`, `kegg.pathway:hsa04115`, `biogrid_interaction/<id>` |
| organisms | `NCBITaxon:<id>` |

Identifiers that a source asserts to be the same molecule, such as secondary accessions and
STRING or IntAct ids, go to `external_ids`. `populate!` and `integrate!` use them to merge.
Related identifiers that are not equivalent, such as a protein's gene or transcripts, go to
`db_xrefs`. To make one source's ids meet another's, supply the mapping: UniProt's `STRING`
column, or `fetch_kegg_uniprot` for KEGG.

### Relationship types that fit the ontology

Several relationship types restrict their participants. `:Ubiquitination` needs
protein → protein, and `:RegulateExpression` needs a gene as its target. The databases don't
always fit those rules: an E3 ligase *complex* ubiquitinates, and a transcription factor
"up-regulates the quantity" of a *protein*. Every edge therefore passes through `fit_type`.
When the participants don't fit, it falls back to the nearest broader type that does
(`FALLBACK`): `:RegulateExpression` → `:Increase`, `:Ubiquitination` by a complex →
`:StateTransition`, `:Catalysis` by a small molecule → `:Activation`, and finally
`:Association`. This is why every imported graph passes `validate`. Pass `strict = false` to
any reader to keep the source's own type instead; `validate` then reports the mismatches.

Edges are `Interaction`s for symmetric associations (binding, association, co-localization)
and directed `Influence`s otherwise. Reactions are `Reaction`s.

### Per database

**UniProt:** `Protein` with its name, synonyms, organism, sequence, cellular locations (from
the subcellular-location text and GO cellular-component terms), modified residues (`"Ser15"`)
and cofactors (`SmallMolecule`s). With `genes = true` (the default), each protein also gets its
`Gene`, carrying the HGNC id and the GO molecular functions, plus a gene → protein `:Expression`
edge.

**KEGG (`read_kgml`):**
- **Gene boxes** become `Protein`s (KGML defines them as gene products), one per gene in the box.
- **Compounds** become `SmallMolecule`s, **ortholog boxes** `OrthologGroup`s, **groups** `MolecularComplex`es, and **linked maps** `PathwayModule`s.
- **Relations** give one edge per subtype for every pair of products: activation, inhibition, expression, phosphorylation, binding/association, …
- **ECrel relations** become `:SuccessiveReaction`.
- **Reactions** become `Reaction`s, with `:Catalysis` from the enzymes.
- **`genes = true`** adds `Gene`s, so "expression" can target genes as the ontology intends.

**PSI-MITAB (IntAct, …):**
- **Interactors** become nodes according to their PSI-MI type: protein, small molecule, gene, DNA, RNA (with its RNA type, such as miRNA or mRNA), or complex.
- **Interaction types** map to relationship types (`MI_INTERACTION`): direct interaction → `:Binding`, physical association → `:Association`, colocalization → `:CoOccurence`, phosphorylation/ubiquitination/… reactions → those types, other enzymatic reactions → `:StateTransition`.
- **MITAB 2.8 causal statements** (`MI_CAUSAL`) give directed `Influence`s. Biological roles (enzyme / enzyme target, …) set the direction.
- **Kept:** publications, confidence scores (such as the IntAct MI-score) and interaction ids.
- **Skipped:** negative interactions. `organisms = [9606]` keeps only human–human interactions.

**BioGRID:**
- **Physical interactions** connect proteins. BioGRID's method-specific scores go to `weight`, not `confidence`. The experimental system sets the type (`BIOGRID_SYSTEM`): two-hybrid, co-crystal structure and other direct methods give `:Binding`, affinity capture and co-fractionation give `:Association`, and co-localization gives `:CoOccurence`.
- **"Biochemical Activity"** rows that report a modification become directed edges, such as `:Phosphorylation`.
- **Genetic interactions** (`genetic = true`) connect `Gene`s.
- **`read_biogrid_chemicals`** links drugs to their targets: inhibitors and antagonists give `:Inhibition`, agonists give `:Activation`.

**SIGNOR:**
- **Edge type:** each row becomes a directed `Influence`. A mechanism that names a modification (phosphorylation, ubiquitination, …) sets the type; otherwise the effect does (`SIGNOR_EFFECT`).
- **Other fields:** the modified residue is added to the target's `modification_sites`, the SIGNOR id becomes the edge id, and the score becomes `confidence`.
- **Not imported:** protein families, phenotypes and stimuli, which have no BioNExt class. They are counted in the log.
- **Complexes:** pass `complexes = read_signor_complexes(file)` to fill in complex subunits.

**BioPAX (Reactome, BioCyc, Pathway Commons):** BioNExt was modelled on BioPAX, so the
mapping is direct. See `?read_biopax` for the full table. Proteins and small molecules are named
by their entity reference (UniProt, ChEBI). With `states = :collapse` (the default), the
modified and relocated forms of a molecule become one node, holding all its locations and
modification sites. `states = :keep` keeps one node per BioPAX entity instead. Entity sets are
expanded to their members. Catalyses and controls become `Influence`s on the reactions they
control.

### Provenance and skipped records

Each reader adds a `Provenance` record (`prov/<source>`), giving the source and file and the
date collected, and attaches it to every entity it creates. After `populate!`, an entity's
`provenance` lists every database that described it. Records a reader can't represent are
counted and reported in one log message per import, for example:

```
[ Info: SIGNOR import skipped: 1 × entity of type "stimulus"; 1 × relation with a skipped entity
```

## Combining sources

`populate!(kg, graphs...)` calls `integrate!` from BioNExt for each graph in turn:

- **Same entity:** entities of the same class that share an id or an `external_ids` entry
  become one. The surviving id is chosen by `prefer`, which defaults to `DEFAULT_PREFER`:
  UniProt, then ChEBI, NCBI Gene, HGNC, Ensembl, …
- **Same edge:** edges that become identical (same class, type and participants) become one,
  with the publications, provenance and highest confidence of all their sources (`keep_best`).
- **Result:** a `Resolution` listing what was merged and every conflicting value.

Read each source with the same `isoforms` setting, so that protein ids line up.

## Size and speed

Readers stream their input. BioGRID's Fanconi-anaemia project file, with 156,669 interactions
among 16,695 proteins, imports in about 45 s. Whole-organism BioGRID or IntAct releases are
several times larger; I haven't measured them, but expect minutes and several GB of memory. Call `index!(kg)` on the result before running many queries.

## Adding a database

A reader is a function that returns a `KnowledgeGraph`. The helpers in `src/common.jl` do
most of the work:

- `curie(db, id)` and `protein_id(acc)`: consistent identifiers
- `upsert!(kg, resource)`: add, or merge into an existing resource
- `edge!(kg, type, left, right; source, kwargs...)`: an edge whose type is fitted to the ontology and whose class follows from its type
- `taxon!(kg, taxid)`, `location_terms(x)`, `residue_site(res, pos)`
- `source_provenance!` and `stamp!`: provenance
- `with_source(f, src)`: plain, gzip and zip input
- `fetch_url(url; cache)`: cached downloads with retries

Add a test with a small real export in `test/fixtures/`, as the existing readers have.

## Tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'                                   # offline, stored exports
BIONEXT_ONLINE_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'         # also the live APIs
```

The stored exports are small, real extracts of each database. `test/fixtures/README.md` lists
where each one came from. The one exception is the SIGNOR file, which was written by hand.
