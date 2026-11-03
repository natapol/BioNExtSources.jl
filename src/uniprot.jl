# UniProtKB: proteins, their genes, locations, modified residues, cofactors.

"""Fields requested from the UniProt REST API by `fetch_uniprot`."""
const UNIPROT_FIELDS = "accession,id,gene_names,protein_name,organism_id,organism_name,sequence," *
    "cc_subcellular_location,go_c,go_f,xref_geneid,xref_hgnc,xref_string,xref_ensembl,ft_mod_res,cc_cofactor"

escape_uri(s) = join(c in ('-', '_', '.', '~') || isletter(c) && isascii(c) || isdigit(c) ? string(c) :
                     join("%" * uppercase(string(b; base = 16, pad = 2)) for b in codeunits(string(c))) for c in s)

"""
    fetch_uniprot(query; fields = UNIPROT_FIELDS, cache = cache_dir(), refresh = false) -> path

Download all UniProtKB entries matching a UniProt query (e.g. `"organism_id:9606 AND reviewed:true"`,
`"accession:P04637 OR accession:Q00987"`) as TSV. Read it with `read_uniprot`.
"""
function fetch_uniprot(query::AbstractString; fields::AbstractString = UNIPROT_FIELDS, kwargs...)
    # Uncompressed: the stream endpoint answers `compressed=true` with HTTP 500 for some field sets.
    url = "https://rest.uniprot.org/uniprotkb/stream?query=$(escape_uri(query))&format=tsv&fields=$fields"
    return fetch_url(url; name = "uniprot_" * string(hash((query, fields)); base = 16) * ".tsv", kwargs...)
end

"""
    read_uniprot(src; genes = true) -> KnowledgeGraph

Read UniProtKB TSV (as from `fetch_uniprot` or the UniProt website's TSV download).
Only the `Entry` column is required; these are used when present:

| Column | Becomes |
|---|---|
| Entry, Entry Name | `Protein` `UniProt:<acc>`; entry name as a synonym |
| Protein names | `name` (recommended name); alternative names as synonyms |
| Gene Names | synonyms; with `genes = true`, the gene's name |
| Organism (ID), Organism | `organism` (`Taxon NCBITaxon:<id>`) |
| Sequence | `sequence` |
| Subcellular location [CC], Gene Ontology (cellular component) | `cellular_location` |
| Modified residue | `modification_sites` (`"Ser15"`) |
| Cofactor | `cofactors` (`SmallMolecule CHEBI:<id>`) |
| STRING | `external_ids` (`STRING:9606.ENSP…`), so STRING-named data merges |
| Ensembl | `db_xrefs` (transcripts) |
| GeneID, HGNC, Gene Ontology (molecular function) | with `genes = true`: a `Gene` (`NCBIGene:<id>`, else `HGNC:<id>`) with the HGNC id as external id, GO functions in `functions`, and an `Expression` edge gene → protein |
"""
function read_uniprot(src; genes::Bool = true)
    kg = KnowledgeGraph()
    with_source(src) do io
        lines = eachline(io)
        header = String.(split(first(lines), '\t'))
        col = Dict(h => i for (i, h) in enumerate(header))
        haskey(col, "Entry") || throw(ArgumentError("UniProt TSV needs an `Entry` column; got: $(join(header, ", "))"))
        for line in lines
            isempty(strip(line)) && continue
            f = split(line, '\t')
            get_(name) = (i = get(col, name, nothing); i === nothing || i > length(f) ? nothing : cell(f[i]))
            uniprot_row!(kg, get_, genes)
        end
    end
    stamp!(kg, source_provenance!(kg, "UniProt"; url = "https://www.uniprot.org", file = src isa IO ? nothing : src))
    return kg
end

function uniprot_row!(kg, get_, genes)
    acc = get_("Entry")
    p = Protein(id = protein_id(acc))
    if (names = get_("Protein names")) !== nothing
        p.name = String(strip(first(split(names, " ("))))
        for m in eachmatch(r"\(([^()]+)\)", names)
            startswith(m[1], "EC ") || push!(p.synonyms, m[1])
        end
    end
    gene_names = split_list(get_("Gene Names"), ' ')
    append!(p.synonyms, gene_names)
    (e = get_("Entry Name")) === nothing || push!(p.synonyms, e)
    tax = get_("Organism (ID)")
    p.organism = tax === nothing ? nothing : taxon!(kg, tax; name = (o = get_("Organism")) === nothing ? nothing : strip(first(split(o, " ("))))
    p.sequence = get_("Sequence")
    for c in ("Subcellular location [CC]", "Gene Ontology (cellular component)")
        (x = get_(c)) === nothing && continue
        texts = c == "Subcellular location [CC]" ? [x] : [m.match for m in eachmatch(r"GO:\d{7}", x)]
        foreach(t -> foreach(l -> l in p.cellular_location || push!(p.cellular_location, l), location_terms(t)), texts)
    end
    if (mods = get_("Modified residue")) !== nothing
        for m in eachmatch(r"MOD_RES (\d+);\s*/note=\"([^\";]+)", mods)
            pos = parse(Int, m[1])
            site = p.sequence !== nothing && pos <= length(p.sequence) ? residue_site(p.sequence[pos], pos) : residue_site(m[2], pos)
            site in p.modification_sites || push!(p.modification_sites, site)
        end
    end
    if (cof = get_("Cofactor")) !== nothing
        for m in eachmatch(r"Name=([^;]+);\s*Xref=ChEBI:(CHEBI:\d+)", cof)
            upsert!(kg, SmallMolecule(id = m[2], name = strip(m[1])))
            m[2] in p.cofactors || push!(p.cofactors, m[2])
        end
    end
    for s in split_list(get_("STRING"), ';')
        push!(p.external_ids, "STRING:" * s)
    end
    for t in split_list(get_("Ensembl"), ';')
        (m = match(r"^(ENS\w+?)(\.\d+)?\b", t)) === nothing || push!(p.db_xrefs, "Ensembl:" * m[1])
    end
    unique!(p.synonyms); unique!(p.db_xrefs)
    upsert!(kg, p)

    genes || return
    geneids, hgnc = split_list(get_("GeneID"), ';'), split_list(get_("HGNC"), ';')
    gid = !isempty(geneids) ? "NCBIGene:" * geneids[1] : !isempty(hgnc) ? hgnc[1] : nothing
    gid === nothing && return
    g = Gene(id = gid, name = isempty(gene_names) ? nothing : gene_names[1], organism = p.organism,
             synonyms = gene_names[2:end], external_ids = filter(!=(gid), String.(hgnc)))
    (mf = get_("Gene Ontology (molecular function)")) === nothing ||
        append!(g.functions, unique(m.match for m in eachmatch(r"GO:\d{7}", mf)))
    upsert!(kg, g)
    edge!(kg, :Expression, [gid], [p.id]; source = "uniprot")
end
