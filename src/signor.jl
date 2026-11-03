# SIGNOR: causal signalling relationships (native TSV export) and complexes.
# SIGNOR's causalTab export is PSI-MITAB 2.8: read it with `read_mitab(src; source = "SIGNOR")`.

# Column order of SIGNOR's TSV export (used when the file has no header line).
const SIGNOR_COLUMNS = ["ENTITYA", "TYPEA", "IDA", "DATABASEA", "ENTITYB", "TYPEB", "IDB", "DATABASEB",
    "EFFECT", "MECHANISM", "RESIDUE", "SEQUENCE", "TAX_ID", "CELL_DATA", "TISSUE_DATA", "MODULATOR_COMPLEX",
    "TARGET_COMPLEX", "MODIFICATIONA", "MODASEQ", "MODIFICATIONB", "MODBSEQ", "PMID", "DIRECT", "NOTES",
    "ANNOTATOR", "SENTENCE", "SIGNOR_ID", "SCORE"]

# Mechanisms that name a modification => relationship type (they take precedence over the effect).
const SIGNOR_MECHANISM = Dict(
    "phosphorylation" => :Phosphorylation, "dephosphorylation" => :Dephosphorylation,
    "ubiquitination" => :Ubiquitination, "methylation" => :Methylation, "glycosylation" => :Glycosylation,
    "transcriptional regulation" => nothing, "binding" => nothing,
)
# Effects => relationship type.
const SIGNOR_EFFECT = Dict(
    "up-regulates" => :Increase, "up-regulates activity" => :Activation, "up-regulates quantity" => :Increase,
    "up-regulates quantity by expression" => :RegulateExpression, "up-regulates quantity by stabilization" => :Increase,
    "down-regulates" => :Decrease, "down-regulates activity" => :Inhibition, "down-regulates quantity" => :Decrease,
    "down-regulates quantity by repression" => :Repression,
    "down-regulates quantity by destabilization" => :Degradation, "form complex" => :Binding,
)
# SIGNOR entity types => class; others (protein families, phenotypes, stimuli, ...) are skipped.
const SIGNOR_TYPE = Dict("protein" => Protein, "chemical" => SmallMolecule, "smallmolecule" => SmallMolecule,
                         "complex" => MolecularComplex, "mirna" => Rna, "rna" => Rna, "antibody" => Protein)

"""
    read_signor_complexes(src) -> Dict{String,Tuple{String,Vector{String}}}

Read SIGNOR's complex list (`SIGNOR_complexes.csv`: `SIGNOR ID;COMPLEX NAME;LIST OF ENTITIES`)
into `"SIGNOR:SIGNOR-C1" => (name, ["UniProt:P…", …])`, for `read_signor(...; complexes)`.
"""
function read_signor_complexes(src)
    out = Dict{String,Tuple{String,Vector{String}}}()
    with_source(src) do io
        for (n, line) in enumerate(eachline(io))
            f = split(line, ';')
            length(f) < 3 && continue
            n == 1 && occursin("ID", uppercase(f[1])) && continue          # header line
            id = curie("signor", strip(f[1], ['"', ' ']))
            out[id] = (String(strip(f[2], ['"', ' '])), [protein_id(m) for m in split_list(strip(f[3], '"'), ',')])
        end
    end
    return out
end

function signor_entity!(kg, row, side, complexes, isoforms, skipped)
    typ = lowercase(something(cell(row("TYPE" * side)), ""))
    T = get(SIGNOR_TYPE, typ, nothing)
    T === nothing && (skip!(skipped, "entity of type \"$typ\""); return nothing)
    raw, db = cell(row("ID" * side)), lowercase(something(cell(row("DATABASE" * side)), ""))
    raw === nothing && return nothing
    id = T === Protein && db in ("uniprot", "") ? protein_id(raw; isoforms) :
         db == "pubchem" ? curie("pubchem compound", replace(raw, r"^CID:" => "")) :
         db == "signor" ? curie("signor", raw) : curie(db, raw)
    r = T(; id, name = cell(row("ENTITY" * side)))
    if T === MolecularComplex && haskey(complexes, id)
        name, members = complexes[id]
        r.name = something(r.name, name)
        for m in members
            upsert!(kg, Protein(id = m))
            push!(r.subunits, m)
        end
    end
    upsert!(kg, r)
    return id
end

"""
    read_signor(src; complexes = Dict(), isoforms = :merge, strict = true) -> KnowledgeGraph

Read SIGNOR's TSV export (with or without its header line). Each row becomes a directed
`Influence` from entity A to entity B:

- the **type** comes from the mechanism when it names a modification (phosphorylation,
  ubiquitination, ...), otherwise from the effect (`SIGNOR_EFFECT`: "up-regulates
  activity" → `:Activation`, "down-regulates quantity by repression" → `:Repression`,
  "form complex" → `:Binding`, ...);
- the modified **residue** is added to the target's `modification_sites`;
- `SIGNOR_ID` is the edge id, `SCORE` its `confidence`, `PMID` its publication.

Proteins, chemicals, complexes and miRNAs are imported; protein families, phenotypes and
stimuli have no BioNExt class and are skipped (counted in the log). Pass `complexes =
read_signor_complexes(file)` to fill in complex subunits.
"""
function read_signor(src; complexes = Dict{String,Tuple{String,Vector{String}}}(), isoforms::Symbol = :merge,
                     strict::Bool = true)
    kg = KnowledgeGraph()
    skipped = Skipped()
    with_source(src) do io
        lines = Iterators.Stateful(eachline(io))
        firstline = peek(lines)
        header = occursin("ENTITYA", firstline) ? (popfirst!(lines); String.(strip.(split(firstline, '\t')))) : SIGNOR_COLUMNS
        col = Dict(h => i for (i, h) in enumerate(header))
        for line in lines
            isempty(strip(line)) && continue
            f = split(line, '\t')
            row(name) = (i = get(col, name, nothing); i === nothing || i > length(f) ? nothing : f[i])
            a = signor_entity!(kg, row, "A", complexes, isoforms, skipped)
            b = signor_entity!(kg, row, "B", complexes, isoforms, skipped)
            (a === nothing || b === nothing) && (skip!(skipped, "relation with a skipped entity"); continue)
            mech = lowercase(something(cell(row("MECHANISM")), ""))
            effect = lowercase(something(cell(row("EFFECT")), ""))
            type = something(get(SIGNOR_MECHANISM, mech, nothing), get(SIGNOR_EFFECT, effect, :Association))
            res = cell(row("RESIDUE"))                       # already "Ser15"-style
            if res !== nothing && kg[b] isa Protein && !(res in kg[b].modification_sites)
                push!(kg[b].modification_sites, res)
            end
            sid = cell(row("SIGNOR_ID"))
            pubs = String[p for p in (pmid(something(cell(row("PMID")), "")),) if p !== nothing]
            t = taxon!(kg, something(cell(row("TAX_ID")), "-"))
            if t !== nothing
                for x in (a, b)
                    kg[x] isa Protein && kg[x].organism === nothing && (kg[x].organism = t)
                end
            end
            edge!(kg, type, [a], [b]; id = sid === nothing ? nothing : curie("signor", sid), source = "signor",
                  strict, confidence = tryparse(Float64, something(cell(row("SCORE")), "")), publications = pubs)
        end
    end
    report(skipped, "SIGNOR")
    stamp!(kg, source_provenance!(kg, "SIGNOR"; url = "https://signor.uniroma2.it", file = src isa IO ? nothing : src))
    return kg
end

"""
    fetch_signor(; organism = 9606, cache = cache_dir(), refresh = false) -> path

SIGNOR's full TSV export for one organism (9606 human, 10090 mouse, 10116 rat), from
`https://signor.uniroma2.it/getData.php`.
"""
fetch_signor(; organism::Integer = 9606, kwargs...) =
    fetch_url("https://signor.uniroma2.it/getData.php?organism=$organism"; name = "signor_$organism.tsv", kwargs...)
