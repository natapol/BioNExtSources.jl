# KEGG: pathway maps (KGML) and identifier conversion through the KEGG REST API.
# KEGG data are free for academic use only; see https://www.kegg.jp/kegg/legal.html.

# KGML relation subtypes => relationship type.
const KEGG_SUBTYPE = Dict(
    "activation" => :Activation, "inhibition" => :Inhibition, "expression" => :RegulateExpression,
    "repression" => :Repression, "indirect effect" => :Indirect, "state change" => :StateTransition,
    "binding/association" => :Binding, "dissociation" => :Dissociation,
    "missing interaction" => :MissingInteraction, "phosphorylation" => :Phosphorylation,
    "dephosphorylation" => :Dephosphorylation, "glycosylation" => :Glycosylation,
    "ubiquitination" => :Ubiquitination, "methylation" => :Methylation,
)

"""
    read_kegg_conv(src) -> Dict{String,Vector{String}}

Read KEGG REST `conv` output (`hsa:7157<TAB>up:P04637`) into KEGG gene ⇒ target ids
(`"hsa:7157" => ["P04637", ...]`), e.g. for `read_kgml(...; uniprot)`.
"""
function read_kegg_conv(src)
    out = Dict{String,Vector{String}}()
    with_source(src) do io
        for line in eachline(io)
            f = split(strip(line), '\t')
            length(f) == 2 || continue
            push!(get!(out, String(f[1]), String[]), replace(String(f[2]), r"^[a-z]+:" => ""))
        end
    end
    return out
end

"""
    fetch_kegg(path; cache = cache_dir(), refresh = false) -> file

Any KEGG REST operation, e.g. `fetch_kegg("get/hsa04115/kgml")`, `fetch_kegg("conv/uniprot/hsa")`,
`fetch_kegg("list/pathway/hsa")`.
"""
fetch_kegg(op::AbstractString; kwargs...) =
    fetch_url("https://rest.kegg.jp/" * op; name = "kegg_" * replace(op, r"[^A-Za-z0-9._-]+" => "_"), kwargs...)

"""    fetch_kgml(pathway; kwargs...) -> file

KGML of a KEGG pathway (`"hsa04115"`, `"sce04111"`, `"map00010"`)."""
fetch_kgml(pathway::AbstractString; kwargs...) = fetch_kegg("get/$(replace(pathway, "path:" => ""))/kgml"; kwargs...)

"""    fetch_kegg_uniprot(organism; kwargs...) -> Dict{String,Vector{String}}

KEGG gene ⇒ UniProt accessions for a whole KEGG organism (`"hsa"`, `"sce"`, `"eco"`)."""
fetch_kegg_uniprot(org::AbstractString; kwargs...) = read_kegg_conv(fetch_kegg("conv/uniprot/$org"; kwargs...))

"""
    read_kgml(src; uniprot = Dict(), genes = false, taxon = nothing, maplinks = false,
              strict = true) -> KnowledgeGraph

Read a KEGG pathway map (KGML). The map becomes a `PathwayModule` (`kegg.pathway:hsa04115`)
listing everything below as components.

- **Gene boxes** are gene products, so they become `Protein`s, one per KEGG gene in the box
  (KEGG puts paralogs that can fill the same role in one box). A protein is named by its
  UniProt accession when `uniprot` (KEGG gene ⇒ accessions, e.g. from
  `fetch_kegg_uniprot("hsa")`) knows it, and `gene_product/kegg.genes:hsa:7157` otherwise.
  With `genes = true`, each also gets its `Gene` (`kegg.genes:hsa:7157`) and an
  `Expression` edge; this lets "expression"/"repression" relations target genes as the
  ontology expects (without genes they become `:Increase` / `:Decrease`).
- **Compounds** and glycans become `SmallMolecule`s (`kegg.compound:C00002`), **ortholog
  boxes** `OrthologGroup`s (`kegg.orthology:K…`), **groups** `MolecularComplex`es whose
  subunits are the member products, **linked maps** `PathwayModule`s.
- **Relations** become one edge per subtype (`KEGG_SUBTYPE`) for every pair of products in
  the two boxes. ECrel relations become `:SuccessiveReaction`; relations to linked maps
  are skipped unless `maplinks = true`.
- **Reactions** become `Reaction`s (`kegg.reaction:R00710`, `:Reversible` or
  `:Irreversible`) from substrates to products, catalysed (`:Catalysis`) by the products
  of the enzyme box.
- `taxon` (an NCBI taxonomy id) sets the `organism` of the proteins and genes.
"""
function read_kgml(src; uniprot = Dict{String,Vector{String}}(), genes::Bool = false, taxon = nothing,
                   maplinks::Bool = false, strict::Bool = true)
    kg = KnowledgeGraph()
    skipped = Skipped()
    doc = with_source(io -> EzXML.parsexml(read(io, String)), src)
    root = EzXML.root(doc)
    org = taxon === nothing ? nothing : taxon!(kg, taxon)
    pw = PathwayModule(id = curie("kegg pathway", replace(root["name"], "path:" => "")),
                       name = haskey(root, "title") ? root["title"] : nothing)
    members = Dict{String,Vector{String}}()      # KGML entry id => node ids
    entries = EzXML.findall("entry", root)

    function product!(keggid, label)
        accs = get(uniprot, keggid, String[])
        gene = curie("kegg.genes", keggid)
        id = isempty(accs) ? gene_product_id(gene) : protein_id(accs[1])
        upsert!(kg, Protein(; id, name = label, organism = org, db_xrefs = [gene],
                            external_ids = [protein_id(a) for a in accs[2:end]]))
        if genes
            upsert!(kg, Gene(id = gene, name = label, organism = org))
            edge!(kg, :Expression, [gene], [id]; source = "kegg", strict)
        end
        return id
    end
    label_of(e) = (g = EzXML.findfirst("graphics", e); g === nothing || !haskey(g, "name") ? nothing :
                   String(strip(first(split(replace(g["name"], "..." => "", r"^TITLE:" => ""), ',')))))

    for e in entries
        eid, typ, names = e["id"], e["type"], split(e["name"])
        typ == "map" && names == [root["name"]] && continue      # the map's own title box
        ids = if typ == "gene"
            [product!(n, length(names) == 1 ? label_of(e) : nothing) for n in names]
        elseif typ in ("compound", "glycan")
            [(upsert!(kg, SmallMolecule(id = curie(n[1:findfirst(':', n)-1], n), name = label_of(e))); curie(n[1:findfirst(':', n)-1], n))
             for n in names]
        elseif typ == "ortholog"
            [(upsert!(kg, OrthologGroup(id = curie("ko", n), name = label_of(e))); curie("ko", n)) for n in names]
        elseif typ == "map"
            [(upsert!(kg, PathwayModule(id = curie("path", replace(n, "path:" => "")), name = label_of(e)));
              curie("path", replace(n, "path:" => ""))) for n in names]
        else
            String[]
        end
        members[eid] = ids
    end
    for e in entries                                    # groups refer to other entries
        e["type"] == "group" || continue
        subunits = unique(vcat((get(members, c["id"], String[]) for c in EzXML.findall("component", e))...))
        isempty(subunits) && continue
        id = "kegg_complex/" * join(sort(subunits), "+")
        upsert!(kg, MolecularComplex(; id, subunits))
        members[e["id"]] = [id]
    end

    gene_target(id) = genes && kg[id] isa Protein && !isempty(kg[id].db_xrefs) ? kg[id].db_xrefs[1] : id
    for r in EzXML.findall("relation", root)
        a, b = get(members, r["entry1"], String[]), get(members, r["entry2"], String[])
        (isempty(a) || isempty(b)) && (skip!(skipped, "relation with an unsupported entry"); continue)
        if !maplinks && r["type"] == "maplink"
            skip!(skipped, "map link (pass maplinks = true)"); continue
        end
        subtypes = [s["name"] for s in EzXML.findall("subtype", r)]
        types = r["type"] == "ECrel" ? [:SuccessiveReaction] :
                [t for s in subtypes for t in (get(KEGG_SUBTYPE, s, s == "compound" ? :CompoundMediatedStep : nothing),) if t !== nothing]
        isempty(types) && (types = [:Association])
        for t in unique(types), x in a, y in b
            target = t in (:RegulateExpression, :Repression) ? gene_target(y) : y
            push!(pw.components, edge!(kg, t, [x], [target]; source = "kegg", strict).id)
        end
    end
    for rx in EzXML.findall("reaction", root)
        rids = [curie("rn", n) for n in split(rx["name"])]
        sub = [x for s in EzXML.findall("substrate", rx) for x in get(members, s["id"], String[])]
        prod = [x for p in EzXML.findall("product", rx) for x in get(members, p["id"], String[])]
        r = upsert!(kg, Reaction(id = rids[1], type = rx["type"] == "reversible" ? :Reversible : :Irreversible,
                                 left = sub, right = prod, external_ids = rids[2:end]))
        push!(pw.components, r.id)
        for enzyme in get(members, rx["id"], String[])
            kg[enzyme] isa OrthologGroup && (skip!(skipped, "catalysis by an ortholog group"); continue)
            push!(pw.components, edge!(kg, :Catalysis, [enzyme], [r.id]; source = "kegg", strict).id)
        end
    end
    for ids in values(members), id in ids
        kg[id] isa PathwayModule || id in pw.components || push!(pw.components, id)
    end
    unique!(pw.components)
    upsert!(kg, pw)
    report(skipped, "KEGG")
    stamp!(kg, source_provenance!(kg, "KEGG"; url = "https://www.kegg.jp", file = src isa IO ? nothing : src))
    return kg
end
