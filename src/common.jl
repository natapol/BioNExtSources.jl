# Shared machinery for all importers: identifiers, vocabulary mapping, building
# resources and edges that satisfy the ontology, provenance, file access and downloads.

# --- identifiers ------------------------------------------------------------------

# Database names as written by the sources (lower-cased) => CURIE prefix used in the graph.
# Prefixes follow identifiers.org, so ids resolve and match across sources.
const DB_PREFIX = Dict(
    "uniprotkb" => "UniProt", "uniprot" => "UniProt", "uniprot knowledgebase" => "UniProt",
    "swiss-prot" => "UniProt", "uniprot/swiss-prot" => "UniProt", "uniprot/trembl" => "UniProt",
    "trembl" => "UniProt", "sp" => "UniProt", "up" => "UniProt",
    "chebi" => "CHEBI",
    "entrez gene/locuslink" => "NCBIGene", "entrezgene" => "NCBIGene", "entrez gene" => "NCBIGene",
    "ncbi gene" => "NCBIGene", "geneid" => "NCBIGene", "ncbigene" => "NCBIGene", "entrez" => "NCBIGene",
    "ncbi-geneid" => "NCBIGene",
    "ensembl" => "Ensembl", "refseq" => "RefSeq", "intact" => "IntAct", "reactome" => "Reactome",
    "pubmed" => "PMID", "go" => "GO", "gene ontology" => "GO",
    "complex portal" => "ComplexPortal", "complexportal" => "ComplexPortal",
    "signor" => "SIGNOR", "pubchem" => "pubchem.compound", "pubchem compound" => "pubchem.compound",
    "pubchem-compound" => "pubchem.compound", "cid" => "pubchem.compound",
    "drugbank" => "DrugBank", "hgnc" => "HGNC", "psi-mi" => "MI", "mi" => "MI",
    "psi-mod" => "MOD", "mod" => "MOD", "ec" => "EC", "ec-code" => "EC", "enzyme consortium" => "EC",
    "kegg compound" => "kegg.compound", "kegg.compound" => "kegg.compound", "cpd" => "kegg.compound",
    "kegg glycan" => "kegg.glycan", "gl" => "kegg.glycan", "kegg drug" => "kegg.drug", "dr" => "kegg.drug",
    "kegg.genes" => "kegg.genes", "kegg pathway" => "kegg.pathway", "path" => "kegg.pathway",
    "kegg reaction" => "kegg.reaction", "rn" => "kegg.reaction", "ko" => "kegg.orthology",
    "biogrid" => "BioGRID", "mirbase" => "miRBase", "mirbase mature" => "miRBase",
    "taxonomy" => "NCBITaxon", "taxid" => "NCBITaxon", "ncbi taxonomy" => "NCBITaxon",
    "rnacentral" => "RNAcentral", "pdb" => "PDB", "doi" => "DOI", "cas" => "CAS",
    "sgd" => "SGD", "mgi" => "MGI", "wormbase" => "WormBase", "flybase" => "FlyBase",
)

"""
    curie(db, id) -> String

Identifier as a CURIE with the prefix used throughout BioNExtSources, e.g.
`curie("uniprotkb", "P04637") == "UniProt:P04637"`, `curie("ChEBI", "CHEBI:15422") == "CHEBI:15422"`.
Unknown database names are kept as given.
"""
function curie(db::AbstractString, id::AbstractString)
    d = lowercase(strip(db))
    prefix = get(DB_PREFIX, d, String(strip(db)))
    local_ = String(strip(id, ['"', ' ']))
    m = match(r"^([A-Za-z][\w.\-]*):(.+)$", local_)
    m !== nothing && lowercase(m[1]) in (lowercase(prefix), d) && (local_ = m[2])
    return prefix * ":" * local_
end

"""
    protein_id(accession; isoforms = :merge) -> String

`"UniProt:P04637"` for a UniProt accession. With `isoforms = :merge` (default), isoform
suffixes are removed (`P04637-2` → `P04637`), so all isoforms become one node; with
`:keep` they stay separate (see `isoform_parent`).
"""
function protein_id(acc::AbstractString; isoforms::Symbol = :merge)
    a = replace(String(strip(acc)), r"^(?i:uniprot(kb)?|sp|up):" => "")
    isoforms === :merge && (a = replace(a, r"-\d+$" => ""))
    return "UniProt:" * a
end

"""    isoform_parent(id) -> Union{String,Nothing}

Canonical protein of an isoform id (`"UniProt:P04637-2"` → `"UniProt:P04637"`)."""
isoform_parent(id::AbstractString) = (m = match(r"^(UniProt:[^-]+)-\d+$", id)) === nothing ? nothing : String(m[1])

"""
    gene_product_id(gene_id) -> String

Id for "the product of this gene" when a source names gene products by their gene
(KEGG, BioGRID) and no UniProt accession is known: `"gene_product/NCBIGene:7157"`.
It is minted in the graph's own namespace and carries the gene id in `db_xrefs`.
"""
gene_product_id(gene_id::AbstractString) = "gene_product/" * gene_id

pmid(x) = (m = match(r"(\d+)", string(x))) === nothing ? nothing : "PMID:" * m[1]

# --- vocabulary mapping ---------------------------------------------------------------

go_id(iri::AbstractString) = (m = match(r"GO[_:](\d{7})", iri)) === nothing ? nothing : "GO:" * m[1]

# GO cellular component => location term: the ontology's own seeAlso links, plus common
# sub-compartments mapped to the enclosing compartment the vocabulary has.
const GO_LOCATION = let d = Dict{String,Symbol}()
    for t in terms(:CellularLocation), u in term(t).see_also
        (g = go_id(u)) === nothing || (d[g] = t)
    end
    merge!(d, Dict(
        "GO:0005654" => :Nucleus, "GO:0005635" => :Nucleus, "GO:0031965" => :Nucleus,
        "GO:0005829" => :Cytoplasm, "GO:0005759" => :Mitochondria, "GO:0005743" => :Mitochondria,
        "GO:0005741" => :Mitochondria, "GO:0031966" => :Mitochondria, "GO:0005758" => :Mitochondria,
        "GO:0000139" => :GolgiBody, "GO:0005796" => :GolgiBody, "GO:0010008" => :Endosome,
        "GO:0031901" => :Endosome, "GO:0005769" => :Endosome, "GO:0005770" => :Endosome,
        "GO:0070062" => :ExtracellularVesicle, "GO:0005778" => :Peroxisome, "GO:0005782" => :Peroxisome,
        "GO:0005774" => :Vacuole, "GO:0000329" => :Vacuole, "GO:0000324" => :Vacuole,
        "GO:0030867" => :RoughER, "GO:0030868" => :SmoothER, "GO:0031410" => :Vesicle,
        "GO:0030659" => :Vesicle, "GO:0009986" => :PlasmaMembrane, "GO:0005887" => :PlasmaMembrane,
    ))
    d
end

# Location names as used by UniProt, Reactome, ... => term. "endoplasmic reticulum" is
# deliberately absent: the vocabulary only has the rough and smooth ER.
const LOCATION_NAMES = let d = Dict{String,Symbol}()
    for t in terms(:CellularLocation)
        d[lowercase(term(t).label)] = t
        d[lowercase(string(t))] = t
    end
    merge!(d, Dict(
        "nucleoplasm" => :Nucleus, "nuclear envelope" => :Nucleus, "nuclear membrane" => :Nucleus,
        "chromatin" => :Nucleus, "cytosol" => :Cytoplasm, "mitochondrion" => :Mitochondria,
        "mitochondrial matrix" => :Mitochondria, "mitochondrial inner membrane" => :Mitochondria,
        "mitochondrial outer membrane" => :Mitochondria, "cell membrane" => :PlasmaMembrane,
        "golgi apparatus" => :GolgiBody, "golgi apparatus membrane" => :GolgiBody, "golgi membrane" => :GolgiBody,
        "lysosome membrane" => :LysosomalMembrane, "lysosomal membrane" => :LysosomalMembrane,
        "peroxisome" => :Peroxisome, "early endosome" => :Endosome, "late endosome" => :Endosome,
        "recycling endosome" => :Endosome, "endosome membrane" => :Endosome,
        "extracellular exosome" => :ExtracellularVesicle, "exosome" => :ExtracellularVesicle,
        "cytoplasmic vesicle" => :Vesicle, "multivesicular body" => :MultivesicularBody,
        "vacuole membrane" => :Vacuole, "nuclear pore" => :NuclearPoreComplex,
        "rough endoplasmic reticulum" => :RoughER, "smooth endoplasmic reticulum" => :SmoothER,
        "ergic" => :ERGogiIntermediateCompartment,
    ))
    d
end

"""
    location_terms(x) -> Vector{Symbol}

Cellular-location terms for a GO id (`"GO:0005654"`), a location name (`"nucleoplasm"`),
or a UniProt "Subcellular location" text. Unmappable locations are dropped.
"""
function location_terms(x::AbstractString)
    g = go_id(x)
    g !== nothing && return haskey(GO_LOCATION, g) ? [GO_LOCATION[g]] : Symbol[]
    text = replace(x, r"SUBCELLULAR LOCATION:\s*" => "", r"\{[^}]*\}" => "", r"Note=.*$"s => "",
                   r"\[[^\]]*\]:?" => "")
    out = Symbol[]
    for part in split(text, r"[.;]")
        # the whole phrase, else its first or last comma-separated piece ("Nucleus, nucleoplasm")
        for cand in (strip(part), strip(first(split(part, ','))), strip(last(split(part, ','))))
            t = get(LOCATION_NAMES, lowercase(cand), nothing)
            t === nothing && continue
            t in out || push!(out, t)
            break
        end
    end
    return out
end

const THREE_LETTER = Dict('A' => "Ala", 'R' => "Arg", 'N' => "Asn", 'D' => "Asp", 'C' => "Cys",
    'E' => "Glu", 'Q' => "Gln", 'G' => "Gly", 'H' => "His", 'I' => "Ile", 'L' => "Leu", 'K' => "Lys",
    'M' => "Met", 'F' => "Phe", 'P' => "Pro", 'S' => "Ser", 'T' => "Thr", 'W' => "Trp", 'Y' => "Tyr",
    'V' => "Val", 'U' => "Sec", 'O' => "Pyl")

const AMINO_ACIDS = Dict("alanine" => "Ala", "arginine" => "Arg", "asparagine" => "Asn",
    "aspartic acid" => "Asp", "aspartate" => "Asp", "cysteine" => "Cys", "glutamic acid" => "Glu",
    "glutamate" => "Glu", "glutamine" => "Gln", "glycine" => "Gly", "histidine" => "His",
    "isoleucine" => "Ile", "leucine" => "Leu", "lysine" => "Lys", "methionine" => "Met",
    "phenylalanine" => "Phe", "proline" => "Pro", "serine" => "Ser", "threonine" => "Thr",
    "tryptophan" => "Trp", "tyrosine" => "Tyr", "valine" => "Val", "selenocysteine" => "Sec")

"""
    residue_site(residue, position) -> String

`"Ser15"` from a one-letter code (`'S'`), or from any name containing an amino acid
(`"O-phospho-L-serine"`, UniProt's `"Phosphoserine"`). Unrecognised residues are kept as given.
"""
function residue_site(res, pos)
    res isa AbstractChar && return get(THREE_LETTER, uppercase(res), string(res)) * string(pos)
    name = lowercase(string(res))
    hits = [k for k in keys(AMINO_ACIDS) if occursin(k, name)]
    isempty(hits) && return string(res, pos)
    return AMINO_ACIDS[argmax(length, hits)] * string(pos)
end

# --- relationship types that respect the ontology -------------------------------------------

# Types that describe a symmetric association: emitted as `Interaction`, others as `Influence`.
const UNDIRECTED = Set([:Binding, :Association, :CoOccurence, :Dissociation, :Dimerization, :CoExpression,
                        :PositiveCorrelation, :NegativeCorrelation, :NoInteraction, :Fusion,
                        :Enhancement, :Suppression])

# When the participants do not match a type's expected classes, fall back step by step
# to a broader type; :Association accepts anything.
const FALLBACK = Dict(
    :Ubiquitination => :StateTransition, :Dephosphorylation => :StateTransition,
    :Dimerization => :Binding, :Binding => :Association,
    :RegulateExpression => :Increase, :Expression => :Increase, :TranscriptionTo => :Increase,
    :Translation => :Increase, :Repression => :Decrease, :Degradation => :Decrease,
    :Catalysis => :Activation, :Activation => :Increase, :Inhibition => :Decrease,
    :Increase => :Association, :Decrease => :Association, :CoOccurence => :Association,
    :CoExpression => :Association, :PositiveCorrelation => :Association,
    :NegativeCorrelation => :Association, :NoInteraction => :Association, :Enhancement => :Association,
    :Suppression => :Association, :Fusion => :Association, :Dissociation => :Association,
    :SignalCrosstalk => :Association,
)

function fits(kg, ids, expected)
    isempty(expected) && return true
    allowed = [BioNExt.CLASS_TYPES[string(c)] for c in expected]
    all(id -> !haskey(kg, id) || any(T -> kg[id] isa T, allowed), ids)
end

"""
    fit_type(kg, type, left, right; strict = true) -> Symbol

`type` if the participants (looked up in `kg`) match its expected classes; otherwise
the nearest broader type that fits (see `FALLBACK`). With `strict = false`, `type`
unchanged.
"""
function fit_type(kg::KnowledgeGraph, t::Symbol, left, right; strict::Bool = true)
    strict || return t
    while true
        x = term(t)
        fits(kg, left, x.expected_left) && fits(kg, right, x.expected_right) && return t
        t === :Association && return t
        t = get(FALLBACK, t, :Association)
    end
end

# --- building the graph -----------------------------------------------------------------------

"""
    upsert!(kg, r) -> resource in kg

Add `r`, or merge it into the resource with the same id (vector fields united, empty
fields filled; existing single values win).
"""
function upsert!(kg::KnowledgeGraph, r::Resource)
    haskey(kg, r.id) || (add!(kg, r); return r)
    existing = kg[r.id]
    typeof(existing) == typeof(r) || return existing          # first class seen wins
    BioNExt.merge_into!(existing, r, Conflict[], keep_first)
    return existing
end

"""    taxon!(kg, taxid; name = nothing) -> Union{String,Nothing}

Ensure a `Taxon` for an NCBI taxonomy id and return its id; negative or missing ids
(PSI-MI uses -1 for "in vitro") give `nothing`."""
function taxon!(kg::KnowledgeGraph, taxid; name = nothing)
    t = tryparse(Int, string(taxid))
    (t === nothing || t <= 0) && return nothing
    id = "NCBITaxon:$t"
    haskey(kg, id) ? (name !== nothing && kg[id].name === nothing && (kg[id].name = String(name))) :
        add!(kg, Taxon(id = id, name = name === nothing ? nothing : String(name)))
    return id
end

"""
    edge!(kg, type, left, right; id = nothing, source, strict = true, kwargs...) -> EdgeEntity

Add an edge whose class follows from its (fitted) type: `Interaction` for symmetric
associations, `Influence` otherwise. `kwargs` set other fields (`confidence`,
`publications`, ...). A missing `id` is minted from source, type and participants, so
the same statement from the same source always gets the same id.
"""
function edge!(kg::KnowledgeGraph, type::Symbol, left::Vector{String}, right::Vector{String};
               id = nothing, source::AbstractString = "edge", strict::Bool = true, kwargs...)
    t = fit_type(kg, type, left, right; strict)
    T = t in UNDIRECTED && all(x -> !(get(kg, x, nothing) isa Union{EdgeEntity,PathwayModule}), right) ?
        Interaction : Influence
    id = something(id, string(lowercase(source), "/", t, "/", join(left, "+"), "/", join(right, "+")))
    e = T(; id, type = t, left = unique(left), right = unique(right))
    for (k, v) in kwargs
        v === nothing || BioNExt.assign!(e, k, v)
    end
    return upsert!(kg, e)
end

"""
    source_provenance!(kg, source; url = nothing, file = nothing, release = nothing) -> String

Add a `Provenance` record for one import and return its id (`"prov/<source>"`, plus the
release if given).
"""
function source_provenance!(kg::KnowledgeGraph, source::AbstractString; url = nothing, file = nothing, release = nothing)
    id = "prov/" * lowercase(source) * (release === nothing ? "" : "/" * string(release))
    srcs = String[x for x in (url, file === nothing ? nothing : "file://" * abspath(file)) if x !== nothing]
    upsert!(kg, Provenance(; id, sources = srcs, collected_date = now()))
    return id
end

# Attach the import's provenance to every entity that came from it.
function stamp!(kg::KnowledgeGraph, prov::String)
    for r in kg
        r isa Entity && !(prov in r.provenance) && push!(r.provenance, prov)
    end
    return kg
end

# Count things an importer had to skip, and report them once at the end.
mutable struct Skipped
    counts::Dict{String,Int}
end
Skipped() = Skipped(Dict{String,Int}())
skip!(s::Skipped, why) = (s.counts[why] = get(s.counts, why, 0) + 1; nothing)
function report(s::Skipped, source)
    isempty(s.counts) && return
    @info "$source import skipped: " * join(("$v × $k" for (k, v) in sort!(collect(s.counts); by = last, rev = true)), "; ")
end

# --- files and downloads -------------------------------------------------------------------

"""
    with_source(f, src; member = nothing)

Call `f(io)` on `src`: an `IO`, a path, a gzip file (`.gz`) or a zip archive (`.zip`, the
first entry whose name contains `member`, or the first entry).
"""
function with_source(f, src; member = nothing)
    src isa IO && return f(src)
    path = String(src)
    if endswith(lowercase(path), ".gz")
        return open(io -> f(GzipDecompressorStream(io)), path)
    elseif endswith(lowercase(path), ".zip")
        z = ZipFile.Reader(path)
        try
            files = filter(e -> !endswith(e.name, "/"), z.files)
            i = member === nothing ? 1 : findfirst(e -> occursin(member, e.name), files)
            i === nothing && throw(ArgumentError("no entry matching \"$member\" in $path"))
            return f(files[i])
        finally
            close(z)
        end
    end
    return open(f, path)
end

"""    cache_dir() -> String

Where downloads are kept: `ENV["BIONEXT_CACHE"]`, or `~/.cache/bionext`."""
cache_dir() = get(ENV, "BIONEXT_CACHE", joinpath(homedir(), ".cache", "bionext"))

"""
    download_retry(url, dest; headers = Dict(), attempts = 4)

`Downloads.download` with retries (after 2, 4 and 8 s): public bioinformatics services
regularly answer a request with a transient error and the next one normally.
"""
function download_retry(url, dest; headers = Dict{String,String}(), attempts::Integer = 4)
    for k in 1:attempts
        try
            return Downloads.download(url, dest; headers)
        catch err
            k == attempts && rethrow()
            dest isa IO || rm(dest; force = true)
            dest isa IOBuffer && take!(dest)
            @debug "retrying $url" exception = err
            sleep(2.0^k)
        end
    end
end

"""
    fetch_url(url; cache = cache_dir(), name = nothing, refresh = false, headers = Dict()) -> path

Download `url` into the cache (once; again only with `refresh = true`) and return the
file path. `name` sets the file name; by default it is derived from the URL. Transient
failures are retried (see `download_retry`).
"""
function fetch_url(url::AbstractString; cache::AbstractString = cache_dir(), name = nothing,
                   refresh::Bool = false, headers = Dict{String,String}())
    mkpath(cache)
    fname = something(name, replace(replace(url, r"^https?://" => ""), r"[^A-Za-z0-9._-]+" => "_"))
    length(fname) > 200 && (fname = first(fname, 150) * "_" * string(hash(url); base = 16))
    path = joinpath(cache, fname)
    (isfile(path) && !refresh) && return path
    tmp = path * ".part"
    try
        download_retry(url, tmp; headers)
    catch err
        rm(tmp; force = true)
        throw(ErrorException("download failed: $url\n$(sprint(showerror, err))"))
    end
    mv(tmp, path; force = true)
    return path
end

split_list(s, sep = '|') = (s === nothing || strip(s) in ("", "-")) ? String[] :
    String[strip(x) for x in split(s, sep) if !isempty(strip(x))]
cell(s) = (s === nothing || strip(s) in ("", "-")) ? nothing : String(strip(s))
