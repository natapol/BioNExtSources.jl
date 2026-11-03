# BioPAX Level 3 (RDF/XML): Reactome, BioCyc, Pathway Commons and other pathway databases.
# BioNExt was modelled on BioPAX, so most classes map one to one.

struct BPRef
    id::String
end

struct BPObject
    class::String
    props::Vector{Pair{String,Any}}         # property => literal String or BPRef
end

attr(e, name) = (i = findfirst(a -> EzXML.nodename(a) == name, collect(EzXML.attributes(e)));
                 i === nothing ? nothing : EzXML.nodecontent(collect(EzXML.attributes(e))[i]))

# Read every BioPAX object (top-level or nested) into id => BPObject.
function biopax_objects(root)
    objs = Dict{String,BPObject}()
    counter = Ref(0)
    function object!(e)
        id = something(attr(e, "about"), (x = attr(e, "ID")) === nothing ? nothing : "#" * x,
                       "_:n" * string(counter[] += 1))
        props = Pair{String,Any}[]
        for p in EzXML.eachelement(e)
            name = EzXML.nodename(p)
            if (r = attr(p, "resource")) !== nothing
                push!(props, name => BPRef(r))
            elseif EzXML.haselement(p)
                push!(props, name => BPRef(object!(first(EzXML.eachelement(p)))))
            else
                push!(props, name => String(strip(EzXML.nodecontent(p))))
            end
        end
        objs[id] = BPObject(EzXML.nodename(e), props)
        return id
    end
    for e in EzXML.eachelement(root)
        EzXML.nodename(e) in ("Ontology",) || object!(e)
    end
    return objs
end

const BP_CONVERSIONS = Set(["BiochemicalReaction", "Conversion", "Transport", "TransportWithBiochemicalReaction",
                            "ComplexAssembly", "Degradation"])
const BP_CONTROLS = Set(["Catalysis", "Control", "Modulation", "TemplateReactionRegulation"])

"""
    read_biopax(src; source = "BioPAX", states = :collapse, isoforms = :merge, strict = true)
        -> KnowledgeGraph

Read a BioPAX Level 3 file (RDF/XML), as exported by Reactome, BioCyc and Pathway Commons.

| BioPAX | BioNExt |
|---|---|
| `Protein`, `SmallMolecule`, `Dna`(`Region`), `Rna`(`Region`), `Gene` | the same classes, named by their entity reference's UniProt / ChEBI / … id |
| `Complex` | `MolecularComplex` with `subunits` and `subcomplexes` |
| entity sets (`memberPhysicalEntity`) | expanded to their members |
| `BiochemicalReaction`, `Transport`, `ComplexAssembly`, other conversions | `Reaction` (`:Reversible` or `:Irreversible`, with EC numbers) |
| `TemplateReaction` | `Reaction` from template to product |
| `Catalysis` | `Influence` `:Catalysis` from controller to reaction |
| `Control`, `Modulation`, `TemplateReactionRegulation` | `Influence` `:Activation` / `:Inhibition` |
| `MolecularInteraction` / `GeneticInteraction` | `Interaction` `:Binding` / `:Association` |
| `Pathway` | `PathwayModule` with its components |
| cellular location (GO), modification features, PubMed xrefs, organisms | `cellular_location`, `modification_sites`, `publications`, `organism` |

Ids are stable database ids when the file has them (Reactome `R-HSA-…`, Pathway Commons
URIs). With `states = :collapse` (default), all forms of a molecule (locations,
modifications) become one node named by its reference id, which is what network analysis
and merging across sources need; `states = :keep` gives one node per BioPAX entity,
linked to the reference id through `db_xrefs`.
"""
function read_biopax(src; source::AbstractString = "BioPAX", states::Symbol = :collapse,
                     isoforms::Symbol = :merge, strict::Bool = true)
    kg = KnowledgeGraph()
    skipped = Skipped()
    root = EzXML.root(with_source(io -> EzXML.parsexml(read(io, String)), src))
    objs = biopax_objects(root)
    fileid = src isa IO ? "stream" : splitext(basename(String(src)))[1]
    get_(o, p) = [v for (k, v) in o.props if k == p]
    lit(o, p) = (v = findfirst(x -> x.first == p && x.second isa String, o.props)) === nothing ? nothing : o.props[v].second
    refs(o, p) = [objs[v.id] for v in get_(o, p) if v isa BPRef && haskey(objs, v.id)]
    refids(o, p) = [v.id for v in get_(o, p) if v isa BPRef && haskey(objs, v.id)]

    xrefs(o, class = nothing) = [(lit(x, "db"), lit(x, "id")) for x in refs(o, "xref")
                                 if (class === nothing || x.class == class) && lit(x, "db") !== nothing && lit(x, "id") !== nothing]
    function uri_curie(u)
        m = match(r"identifiers\.org/([^/]+)/(.+)$", u)
        m === nothing && (m = match(r"bioregistry\.io/([^:/]+):(.+)$", u))
        return m === nothing ? nothing : curie(m[1], m[2])
    end
    function stable_id(oid, o)
        for (db, id) in xrefs(o, "UnificationXref")
            lowercase(db) == "reactome" && return curie("reactome", id)
        end
        startswith(oid, "http") && return oid
        return lowercase(source) * "/" * fileid * "/" * lstrip(oid, ['#', '_', ':'])
    end
    function reference_id(er)
        er === nothing && return nothing
        u = xrefs(er, "UnificationXref")
        for (db, id) in u
            lowercase(db) in ("uniprot", "uniprot knowledgebase", "uniprotkb") && return protein_id(id; isoforms)
        end
        for (db, id) in u
            lowercase(db) == "chebi" && return curie("chebi", id)
        end
        isempty(u) || return curie(u[1]...)
        return nothing
    end
    pubs(o) = String[p for (db, id) in xrefs(o, "PublicationXref") if lowercase(db) == "pubmed" for p in (pmid(id),) if p !== nothing]
    # Names of an entity and of its reference; when states are collapsed the node stands for
    # the reference molecule, so the reference's names come first ("ATM", not "p-S1981-ATM").
    function names(o, er)
        own = [get_(o, "displayName"); get_(o, "standardName"); get_(o, "name")]
        ref = er === nothing ? [] : [get_(er, "displayName"); get_(er, "standardName"); get_(er, "name")]
        all_ = states === :collapse ? [ref; own] : [own; ref]
        return unique!(String[x for x in all_ if x isa String && !isempty(x) && !occursin(r"^\w+:\S+ ", x)])
    end
    function organism(o)
        for b in [refs(o, "organism"); [x for er in refs(o, "entityReference") for x in refs(er, "organism")]]
            for (db, id) in xrefs(b)
                lowercase(db) in ("taxonomy", "ncbi taxonomy", "ncbitaxon") && return taxon!(kg, id; name = lit(b, "displayName"))
            end
        end
        return nothing
    end
    function locations(o)
        out = Symbol[]
        for v in refs(o, "cellularLocation")
            ts = [t for (db, id) in xrefs(v) for t in location_terms(curie(db, id))]
            isempty(ts) && (ts = location_terms(something(lit(v, "term"), "")))
            isempty(ts) && skip!(skipped, "cellular location \"$(something(lit(v, "term"), "?"))\" (no vocabulary term)")
            foreach(t -> t in out || push!(out, t), ts)
        end
        return out
    end
    function sites(o)
        out = String[]
        for f in refs(o, "feature")
            f.class == "ModificationFeature" || continue
            pos = [lit(s, "sequencePosition") for s in refs(f, "featureLocation")]
            mod = [something(lit(m, "term"), "") for m in refs(f, "modificationType")]
            isempty(pos) || pos[1] === nothing || push!(out, residue_site(isempty(mod) ? "" : mod[1], pos[1]))
        end
        return unique!(out)
    end

    entity_cache = Dict{String,Vector{String}}()
    const_types = Dict("Protein" => Protein, "SmallMolecule" => SmallMolecule, "Dna" => Dna, "DnaRegion" => Dna,
                       "Rna" => Rna, "RnaRegion" => Rna, "Gene" => Gene)
    function entity_ids(oid)
        haskey(entity_cache, oid) && return entity_cache[oid]
        entity_cache[oid] = String[]                       # guards against cycles
        o = objs[oid]
        members = refids(o, "memberPhysicalEntity")
        ids = if o.class == "Complex" && isempty(members)
            c = MolecularComplex(id = stable_id(oid, o), name = (n = names(o, nothing); isempty(n) ? nothing : n[1]),
                                 cellular_location = locations(o), organism = organism(o), publications = pubs(o))
            for comp in refids(o, "component")
                target = objs[comp].class == "Complex" ? c.subcomplexes : c.subunits
                foreach(x -> x in target || push!(target, x), entity_ids(comp))
            end
            upsert!(kg, c)
            [c.id]
        elseif !isempty(members) && isempty(refids(o, "entityReference"))
            unique!(vcat(entity_ids.(members)...))
        elseif haskey(const_types, o.class)
            T = const_types[o.class]
            er = (e = refs(o, "entityReference"); isempty(e) ? nothing : e[1])
            ref = something(reference_id(er), er === nothing ? nothing : startswith(refids(o, "entityReference")[1], "http") ?
                            something(uri_curie(refids(o, "entityReference")[1]), refids(o, "entityReference")[1]) : nothing,
                            Some(nothing))
            canonical = ref === nothing ? stable_id(oid, o) : ref
            id = states === :keep ? stable_id(oid, o) : canonical
            n = names(o, er)
            r = T(; id, name = isempty(n) ? nothing : n[1], synonyms = n[2:end], organism = organism(o),
                  cellular_location = locations(o))
            states === :keep && canonical != id && push!(r.db_xrefs, canonical)
            T === Protein && append!(r.modification_sites, sites(o))
            upsert!(kg, r)
            [id]
        else
            skip!(skipped, "physical entity of class $(o.class)")
            String[]
        end
        return entity_cache[oid] = ids
    end

    interaction_cache = Dict{String,Vector{String}}()
    function interaction_ids(oid)
        haskey(interaction_cache, oid) && return interaction_cache[oid]
        interaction_cache[oid] = String[]
        o = objs[oid]
        n = names(o, nothing)
        name = isempty(n) ? nothing : n[1]
        ids = if o.class in BP_CONVERSIONS || o.class == "TemplateReaction"
            left = unique(vcat(entity_ids.(o.class == "TemplateReaction" ? refids(o, "template") : refids(o, "left"))...))
            right = unique(vcat(entity_ids.(o.class == "TemplateReaction" ? refids(o, "product") : refids(o, "right"))...))
            dir = something(lit(o, "conversionDirection"), "")
            dir == "RIGHT-TO-LEFT" && ((left, right) = (right, left))
            if isempty(left) || isempty(right)
                skip!(skipped, "$(o.class) without participants on both sides")
                String[]
            else
                r = Reaction(id = stable_id(oid, o), name = name, type = dir == "REVERSIBLE" ? :Reversible : :Irreversible,
                             left = left, right = right, publications = pubs(o),
                             ec_numbers = ["EC:" * x for x in get_(o, "eCNumber") if x isa String])
                upsert!(kg, r)
                [r.id]
            end
        elseif o.class in BP_CONTROLS
            controllers = String[]
            for c in refids(o, "controller")
                append!(controllers, objs[c].class == "Pathway" ? pathway_ids(c) : entity_ids(c))
            end
            controlled = String[]
            for c in refids(o, "controlled")
                append!(controlled, objs[c].class == "Pathway" ? pathway_ids(c) : interaction_ids(c))
            end
            ct = uppercase(something(lit(o, "controlType"), "ACTIVATION"))
            type = startswith(ct, "INHIBITION") ? :Inhibition : o.class == "Catalysis" ? :Catalysis : :Activation
            base = stable_id(oid, o)
            out = String[]
            for c in controllers, t in controlled
                id = length(controllers) * length(controlled) == 1 ? base : string(base, "/", c, "/", t)
                push!(out, edge!(kg, type, [c], [t]; id, source = lowercase(source), strict, publications = pubs(o)).id)
            end
            isempty(out) && skip!(skipped, "$(o.class) without controller or controlled")
            out
        elseif o.class in ("MolecularInteraction", "GeneticInteraction")
            ps = unique(vcat(entity_ids.(refids(o, "participant"))...))
            length(ps) < 2 ? (skip!(skipped, "$(o.class) with fewer than two participants"); String[]) :
                [edge!(kg, o.class == "GeneticInteraction" ? :Association : :Binding, ps[1:1], ps[2:end];
                       id = stable_id(oid, o), source = lowercase(source), strict, publications = pubs(o)).id]
        elseif o.class == "Pathway"
            pathway_ids(oid)
        else
            skip!(skipped, "interaction of class $(o.class)")
            String[]
        end
        return interaction_cache[oid] = ids
    end

    pathway_cache = Dict{String,String}()
    function pathway_ids(oid)
        haskey(pathway_cache, oid) && return [pathway_cache[oid]]
        o = objs[oid]
        p = PathwayModule(id = stable_id(oid, o), name = (n = names(o, nothing); isempty(n) ? nothing : n[1]),
                          organism = organism(o), publications = pubs(o))
        pathway_cache[oid] = p.id
        upsert!(kg, p)
        for c in refids(o, "pathwayComponent")
            for x in (objs[c].class == "Pathway" ? pathway_ids(c) : interaction_ids(c))
                x in kg[p.id].components || push!(kg[p.id].components, x)
            end
        end
        return [p.id]
    end

    # Pathways first (so every interaction finds its pathway), then anything not in one.
    for (oid, o) in sort!(collect(objs); by = first)
        o.class == "Pathway" && pathway_ids(oid)
    end
    for (oid, o) in sort!(collect(objs); by = first)
        (o.class in BP_CONVERSIONS || o.class in BP_CONTROLS || o.class in ("TemplateReaction", "MolecularInteraction",
                                                                         "GeneticInteraction")) && interaction_ids(oid)
    end
    report(skipped, source)
    stamp!(kg, source_provenance!(kg, source; file = src isa IO ? nothing : src))
    return kg
end

# --- Reactome, BioCyc, Pathway Commons ---------------------------------------------------------

"""
    fetch_reactome_biopax(id; cache = cache_dir(), refresh = false) -> file

BioPAX Level 3 export of a Reactome pathway (stable id `"R-HSA-69541"` or database id `69541`),
including its sub-pathways.
"""
function fetch_reactome_biopax(id; kwargs...)
    dbid = string(id)
    if !occursin(r"^\d+$", dbid)
        buf = IOBuffer()
        download_retry("https://reactome.org/ContentService/data/query/$(escape_uri(dbid))/dbId", buf)
        dbid = strip(String(take!(buf)))
        occursin(r"^\d+$", dbid) || throw(ArgumentError("Reactome has no pathway \"$id\""))
    end
    return fetch_url("https://reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/$dbid";
                     name = "reactome_$(replace(string(id), r"[^A-Za-z0-9-]" => "_")).owl", kwargs...)
end

"""    read_reactome(src; kwargs...) -> KnowledgeGraph

`read_biopax(src; source = "Reactome", kwargs...)` for Reactome BioPAX exports."""
read_reactome(src; kwargs...) = read_biopax(src; source = "Reactome", kwargs...)

"""    read_biocyc(src; kwargs...) -> KnowledgeGraph

`read_biopax(src; source = "BioCyc", kwargs...)` for BioPAX exported from BioCyc (EcoCyc,
YeastCyc, HumanCyc, ...). BioCyc requires an account and blocks scripted downloads:
export a pathway or a whole database as BioPAX Level 3 from the BioCyc website, then read the file."""
read_biocyc(src; kwargs...) = read_biopax(src; source = "BioCyc", kwargs...)

"""
    fetch_pathwaycommons(uri; cache = cache_dir(), refresh = false) -> file

BioPAX for a Pathway Commons object (pathway URI, as found with Pathway Commons search),
covering Reactome, KEGG, PANTHER, HumanCyc and other databases it integrates.
"""
fetch_pathwaycommons(uri::AbstractString; kwargs...) =
    fetch_url("https://www.pathwaycommons.org/pc2/get?uri=$(escape_uri(uri))&format=BIOPAX";
              name = "pc_" * string(hash(uri); base = 16) * ".owl", kwargs...)
