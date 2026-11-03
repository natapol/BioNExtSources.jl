# PSI-MITAB 2.5–2.8: IntAct, MINT, BioGRID (MITAB release), SIGNOR causalTab,
# Reactome interactions and every other PSICQUIC service.

struct MIValue
    db::String
    id::String
    text::String
end

# Split a MITAB field on '|' outside double quotes.
function split_mitab(s::AbstractString)
    out, buf, inq = String[], IOBuffer(), false
    for c in s
        c == '"' && (inq = !inq)
        (c == '|' && !inq) ? push!(out, String(take!(buf))) : print(buf, c)
    end
    push!(out, String(take!(buf)))
    return out
end

"""    mi_values(field) -> Vector{MIValue}

Parse a MITAB field like `psi-mi:"MI:0915"(physical association)|uniprotkb:TP53(gene name)`."""
function mi_values(s::AbstractString)
    s = strip(s)
    (isempty(s) || s == "-") && return MIValue[]
    out = MIValue[]
    for p in split_mitab(s)
        m = match(r"^([^:]+):(\"(?:[^\"\\]|\\.)*\"|[^(]*)(?:\((.*)\))?$", strip(p))
        m === nothing && continue
        push!(out, MIValue(m[1], strip(m[2], '"'), something(m[3], "")))
    end
    return out
end

mi_label(vs) = isempty(vs) ? "" : lowercase(strip(vs[1].text, '"'))
mi_code(vs) = isempty(vs) ? "" : vs[1].id

# PSI-MI interactor types => (class, RNA region term or nothing). Keyed by label and code.
const MI_INTERACTOR = Dict{String,Tuple{Type,Union{Symbol,Nothing}}}(
    "protein" => (Protein, nothing), "MI:0326" => (Protein, nothing),
    "peptide" => (Protein, nothing), "MI:0327" => (Protein, nothing),
    "small molecule" => (SmallMolecule, nothing), "MI:0328" => (SmallMolecule, nothing),
    "polysaccharide" => (SmallMolecule, nothing), "MI:0904" => (SmallMolecule, nothing),
    "xenobiotic" => (SmallMolecule, nothing), "chemical" => (SmallMolecule, nothing),
    "gene" => (Gene, nothing), "MI:0250" => (Gene, nothing),
    "deoxyribonucleic acid" => (Dna, nothing), "MI:0319" => (Dna, nothing),
    "double stranded deoxyribonucleic acid" => (Dna, nothing), "MI:0681" => (Dna, nothing),
    "single stranded deoxyribonucleic acid" => (Dna, nothing), "MI:0680" => (Dna, nothing),
    "nucleic acid" => (Dna, nothing), "MI:0318" => (Dna, nothing),
    "ribonucleic acid" => (Rna, nothing), "MI:0320" => (Rna, nothing),
    "messenger rna" => (Rna, :MessengerRNA), "MI:0324" => (Rna, :MessengerRNA),
    "transfer rna" => (Rna, :TransferRNA), "MI:0325" => (Rna, :TransferRNA),
    "ribosomal rna" => (Rna, :RibosomalRNA), "MI:0608" => (Rna, :RibosomalRNA),
    "small interfering rna" => (Rna, :SmallInterferingRNA), "MI:0610" => (Rna, :SmallInterferingRNA),
    "small nuclear rna" => (Rna, :SmallNuclearRNA), "MI:0607" => (Rna, :SmallNuclearRNA),
    "micro rna" => (Rna, :MicroRNA), "mirna" => (Rna, :MicroRNA), "MI:2204" => (Rna, :MicroRNA),
    "long non-coding ribonucleic acid" => (Rna, :LongNonCodingRNA), "MI:2190" => (Rna, :LongNonCodingRNA),
    "complex" => (MolecularComplex, nothing), "MI:0314" => (MolecularComplex, nothing),
    "protein complex" => (MolecularComplex, nothing), "MI:0315" => (MolecularComplex, nothing),
)

# PSI-MI interaction types => relationship type. Enzymatic ones are directed.
const MI_INTERACTION = Dict{String,Symbol}(
    "direct interaction" => :Binding, "MI:0407" => :Binding,
    "covalent binding" => :Binding, "MI:0195" => :Binding, "disulfide bond" => :Binding, "MI:0408" => :Binding,
    "physical association" => :Association, "MI:0915" => :Association,
    "association" => :Association, "MI:0914" => :Association,
    "colocalization" => :CoOccurence, "MI:0403" => :CoOccurence,
    "genetic interaction" => :Association, "MI:0208" => :Association,
    "phosphorylation reaction" => :Phosphorylation, "MI:0217" => :Phosphorylation,
    "dephosphorylation reaction" => :Dephosphorylation, "MI:0203" => :Dephosphorylation,
    "ubiquitination reaction" => :Ubiquitination, "MI:0220" => :Ubiquitination,
    "methylation reaction" => :Methylation, "MI:0213" => :Methylation,
    "glycosylation reaction" => :Glycosylation, "MI:0559" => :Glycosylation,
)
# Other enzymatic reactions (acetylation, cleavage, sumoylation, ...) change the target's state.
const ENZYMATIC = r"reaction$|cleavage|acetylation|sumoylation|neddylation|palmitoylation|hydroxylation|deamination|demethylation|deacetylation|deubiquitination"

# PSI-MI causal statements (MITAB 2.8, used by SIGNOR) => relationship type.
const MI_CAUSAL = Dict{String,Symbol}(
    "up-regulates" => :Increase, "MI:2235" => :Increase,
    "up-regulates activity" => :Activation, "MI:2236" => :Activation,
    "up-regulates quantity" => :Increase, "MI:2237" => :Increase,
    "up-regulates quantity by expression" => :RegulateExpression, "MI:2238" => :RegulateExpression,
    "up-regulates quantity by stabilization" => :Increase, "MI:2239" => :Increase,
    "down-regulates" => :Decrease, "MI:2240" => :Decrease,
    "down-regulates activity" => :Inhibition, "MI:2241" => :Inhibition,
    "down-regulates quantity" => :Decrease, "MI:2242" => :Decrease,
    "down-regulates quantity by repression" => :Repression, "MI:2243" => :Repression,
    "down-regulates quantity by destabilization" => :Degradation, "MI:2244" => :Degradation,
)

# Biological roles marking the target of a directed interaction.
const TARGET_ROLES = Set(["enzyme target", "MI:0502", "inhibited", "MI:0587", "stimulated", "MI:0844",
                          "regulator target", "MI:2275"])

"""
    mi_relationship(type_values; causal = MIValue[]) -> (type::Symbol, directed::Bool)

Relationship type for a MITAB interaction-type field (and causal-statement field).
"""
function mi_relationship(types::Vector{MIValue}; causal::Vector{MIValue} = MIValue[])
    if !isempty(causal)
        t = get(MI_CAUSAL, mi_label(causal), get(MI_CAUSAL, mi_code(causal), nothing))
        t === nothing || return (t, true)
    end
    label, code = mi_label(types), mi_code(types)
    t = get(MI_INTERACTION, label, get(MI_INTERACTION, code, nothing))
    t !== nothing && return (t, !(t in UNDIRECTED))
    occursin(ENZYMATIC, label) && return (:StateTransition, true)
    return (:Association, false)
end

const SAME_ENTITY_DBS = Set(["uniprotkb", "uniprot", "uniprot/swiss-prot", "uniprot/trembl", "intact",
                             "chebi", "complex portal", "refseq", "ensembl", "signor"])

# Create (or update) the node for one MITAB interactor; return its id or nothing.
function mitab_interactor!(kg, f, side, isoforms, skipped)
    idc, altc, aliasc, taxc, typec = side == 1 ? (1, 3, 5, 10, 21) : (2, 4, 6, 11, 22)
    ids = mi_values(get(f, idc, "-"))
    isempty(ids) && return nothing
    alts = mi_values(get(f, altc, "-"))
    typ = mi_values(get(f, typec, "-"))
    class = get(MI_INTERACTOR, mi_label(typ), get(MI_INTERACTOR, mi_code(typ), nothing))
    all_ids = [ids; alts]
    isuniprot(v) = lowercase(v.db) in ("uniprotkb", "uniprot", "uniprot/swiss-prot", "uniprot/trembl")
    if class === nothing
        any(isuniprot, all_ids) || (skip!(skipped, "interactor of type \"$(mi_label(typ))\""); return nothing)
        class = (Protein, nothing)
    end
    T, region = class
    chosen = T === Protein ? something(findfirst(isuniprot, all_ids), 1) :
             T === SmallMolecule ? something(findfirst(v -> lowercase(v.db) == "chebi", all_ids), 1) : 1
    v = all_ids[chosen]
    id = T === Protein && isuniprot(v) ? protein_id(v.id; isoforms) : curie(v.db, v.id)
    r = T(; id)
    for (k, a) in enumerate(all_ids)
        k == chosen && continue
        db = lowercase(a.db)
        x = isuniprot(a) ? protein_id(a.id; isoforms) : curie(a.db, a.id)
        x == id && continue
        if db in SAME_ENTITY_DBS
            x in r.external_ids || push!(r.external_ids, x)
        elseif !(startswith(x, "NCBIGene:") && !occursin(r"^\d+$", a.id))   # BioGRID lists gene symbols under Entrez
            x in r.db_xrefs || push!(r.db_xrefs, x)
        end
    end
    for a in mi_values(get(f, aliasc, "-"))
        n = strip(a.id)
        if r.name === nothing && a.text == "gene name"
            r.name = n
        else
            n in r.synonyms || push!(r.synonyms, n)
        end
    end
    r.name === nothing && !isempty(r.synonyms) && (r.name = popfirst!(r.synonyms))
    filter!(!=(r.name), r.synonyms)
    tax = mi_values(get(f, taxc, "-"))
    isempty(tax) || (r.organism = taxon!(kg, tax[1].id; name = length(tax) > 1 ? tax[2].text : tax[1].text))
    region === nothing || push!(r.regions, region)
    if T === Protein && isoforms === :keep && (parent = isoform_parent(id)) !== nothing
        upsert!(kg, Protein(id = parent, organism = r.organism))
        push!(r.isoform_of, parent)
    end
    upsert!(kg, r)
    return id
end

"""
    read_mitab(src; source = "IntAct", isoforms = :merge, organisms = nothing, strict = true)
        -> KnowledgeGraph

Read PSI-MITAB 2.5, 2.6, 2.7 or 2.8 (as from IntAct, MINT, BioGRID's MITAB release,
SIGNOR's causalTab, Reactome's interaction export, or `fetch_psicquic`).

- **Interactors** become nodes by their PSI-MI type: proteins are named by UniProt
  accession (from the primary or an alternative id), small molecules by ChEBI, complexes,
  genes, DNA and RNA (with RNA region types) by their primary id. Alternative ids that
  denote the same molecule go to `external_ids`, others to `db_xrefs`; gene-name aliases
  give the `name`; the taxon gives the `organism`.
- **Interactions** become edges: PSI-MI interaction types and, in MITAB 2.8, causal
  statements map to relationship types (`MI_INTERACTION`, `MI_CAUSAL`); physical
  associations become `Interaction`s, enzymatic reactions and causal statements directed
  `Influence`s (A → B, turned around when biological roles say A is the target).
  Publications, confidence scores and interaction ids are kept; spoke-expanded
  interactions get one edge per pair.
- `organisms` (vector of NCBI taxonomy ids) keeps only interactions whose interactors
  are all from those organisms. Negative interactions are skipped.
"""
function read_mitab(src; source::AbstractString = "IntAct", isoforms::Symbol = :merge,
                    organisms = nothing, strict::Bool = true)
    kg = KnowledgeGraph()
    skipped = Skipped()
    keep_tax = organisms === nothing ? nothing : Set(string.(organisms))
    with_source(src) do io
        for line in eachline(io)
            (isempty(strip(line)) || startswith(line, '#') || startswith(line, "ID(s) interactor A")) && continue
            f = split(line, '\t')
            length(f) < 15 && (skip!(skipped, "line with fewer than 15 columns"); continue)
            length(f) >= 36 && lowercase(strip(f[36])) == "true" && (skip!(skipped, "negative interaction"); continue)
            if keep_tax !== nothing
                taxa = [v.id for c in (10, 11) for v in mi_values(f[c])[1:min(end, 1)]]
                all(t -> t in keep_tax, taxa) || (skip!(skipped, "interaction outside `organisms`"); continue)
            end
            a = mitab_interactor!(kg, f, 1, isoforms, skipped)
            b = strip(f[2]) == "-" ? a : mitab_interactor!(kg, f, 2, isoforms, skipped)
            (a === nothing || b === nothing) && (skip!(skipped, "interaction with an unsupported interactor"); continue)
            causal = length(f) >= 46 ? mi_values(f[46]) : MIValue[]
            type, directed = mi_relationship(mi_values(f[12]); causal)
            if directed && length(f) >= 18
                roleA, roleB = mi_values(f[17]), mi_values(f[18])
                (mi_label(roleA) in TARGET_ROLES || mi_code(roleA) in TARGET_ROLES) &&
                    !(mi_label(roleB) in TARGET_ROLES || mi_code(roleB) in TARGET_ROLES) && ((a, b) = (b, a))
            end
            ids = mi_values(f[14])
            primary = isempty(ids) ? nothing :
                      lowercase(ids[1].db) == "biogrid" ? "biogrid_interaction/" * ids[1].id : curie(ids[1].db, ids[1].id)
            expanded = length(f) >= 16 && strip(f[16]) != "-"
            # a spoke-expanded pair is not itself an IntAct record: give it a minted id
            id = primary === nothing ? nothing :
                 expanded ? string(lowercase(source), "/", ids[1].id, "/", a, "+", b) : primary
            pubs = String[]
            for v in mi_values(f[9])
                lowercase(v.db) == "pubmed" && (p = pmid(v.id)) !== nothing && push!(pubs, p)
            end
            conf = nothing
            for v in mi_values(f[15])
                x = tryparse(Float64, v.id)
                x === nothing || (conf = x; break)
            end
            ext = String[curie(v.db, v.id) for v in ids[2:end] if lowercase(v.db) != "biogrid"]
            edge!(kg, type, [a], [b]; id, source, strict, confidence = conf, publications = pubs,
                  external_ids = isempty(ext) ? nothing : ext)
        end
    end
    report(skipped, source)
    stamp!(kg, source_provenance!(kg, source; file = src isa IO ? nothing : src))
    return kg
end

# --- PSICQUIC -------------------------------------------------------------------------------

const PSICQUIC_REGISTRY = "https://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=xml"

"""    psicquic_services() -> Dict{String,String}

Active PSICQUIC services (name ⇒ REST base URL) from the PSICQUIC registry."""
function psicquic_services(; kwargs...)
    out = Dict{String,String}()
    root = EzXML.root(EzXML.readxml(fetch_url(PSICQUIC_REGISTRY; name = "psicquic_registry.xml", kwargs...)))
    for s in EzXML.eachelement(root)
        f = Dict(EzXML.nodename(c) => EzXML.nodecontent(c) for c in EzXML.eachelement(s))
        haskey(f, "name") && haskey(f, "restUrl") &&
            (out[f["name"]] = replace(f["restUrl"], r"^http://www\.ebi\.ac\.uk" => "https://www.ebi.ac.uk"))
    end
    return out
end

"""
    fetch_psicquic(service, query; format = "tab27", page = 2500, limit = 100_000,
                   cache = cache_dir(), refresh = false) -> path

Run a MIQL query (e.g. `"P04637"`, `"species:human AND type:\\"direct interaction\\""`)
against a PSICQUIC service (`"IntAct"`, `"MINT"`, ... see `psicquic_services()`), page
through all results up to `limit` rows, and save them as one MITAB file for `read_mitab`.
"""
function fetch_psicquic(service::AbstractString, query::AbstractString; format::AbstractString = "tab27",
                        page::Integer = 2500, limit::Integer = 100_000, cache::AbstractString = cache_dir(),
                        refresh::Bool = false)
    path = joinpath(cache, "psicquic_" * lowercase(service) * "_" * string(hash((query, format, limit)); base = 16) * ".txt")
    (isfile(path) && !refresh) && return path
    services = psicquic_services(; cache)
    key = something(findfirst(k -> lowercase(k) == lowercase(service), collect(keys(services))), 0)
    key == 0 && throw(ArgumentError("unknown PSICQUIC service \"$service\"; available: $(join(sort!(collect(keys(services))), ", "))"))
    base = rstrip(services[collect(keys(services))[key]], '/') * "/query/"
    mkpath(cache)
    tmp = path * ".part"
    open(tmp, "w") do out
        first = 0
        while first < limit
            n = min(page, limit - first)
            url = base * escape_uri(query) * "?format=$format&firstResult=$first&maxResults=$n"
            buf = IOBuffer()
            download_retry(url, buf)
            chunk = String(take!(buf))
            rows = count(==('\n'), chunk)
            write(out, chunk)
            rows < n && break
            first += n
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""    fetch_intact(query; kwargs...) -> path

`fetch_psicquic("IntAct", query; kwargs...)`."""
fetch_intact(query::AbstractString; kwargs...) = fetch_psicquic("IntAct", query; kwargs...)

"""    read_intact(src; kwargs...) -> KnowledgeGraph

`read_mitab(src; source = "IntAct", kwargs...)`."""
read_intact(src; kwargs...) = read_mitab(src; source = "IntAct", kwargs...)
