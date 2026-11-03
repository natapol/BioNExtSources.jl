# BioGRID: TAB3 interaction files and chemical-association (chemtab) files.
# (BioGRID's PSI-MITAB release can also be read with `read_mitab(src; source = "BioGRID")`.)

# BioGRID experimental systems => relationship type. Systems that show direct contact
# are :Binding; those that only show co-membership in a complex are :Association.
const BIOGRID_SYSTEM = Dict(
    "two-hybrid" => :Binding, "reconstituted complex" => :Binding, "co-crystal structure" => :Binding,
    "far western" => :Binding, "fret" => :Binding, "pca" => :Binding, "protein-peptide" => :Binding,
    "protein-rna" => :Binding, "affinity capture-western" => :Association, "affinity capture-ms" => :Association,
    "affinity capture-luminescence" => :Association, "affinity capture-rna" => :Association,
    "co-fractionation" => :Association, "co-purification" => :Association, "proximity label-ms" => :Association,
    "cross-linking-ms (xl-ms)" => :Association, "co-localization" => :CoOccurence,
)
# Modifications reported with "Biochemical Activity" (A modifies B).
const BIOGRID_MODIFICATION = Dict(
    "phosphorylation" => :Phosphorylation, "dephosphorylation" => :Dephosphorylation,
    "ubiquitination" => :Ubiquitination, "methylation" => :Methylation, "glycosylation" => :Glycosylation,
)

# The protein for a BioGRID interactor: UniProt (Swiss-Prot, else TrEMBL) when given,
# otherwise the product of the Entrez gene.
function biogrid_protein!(kg, row, side, isoforms)
    s = side == 'A' ? "Interactor A" : "Interactor B"
    entrez = cell(row("Entrez Gene $s"))
    sp = split_list(row("SWISS-PROT Accessions $s"))
    tr = split_list(row("TREMBL Accessions $s"))
    gene = entrez === nothing ? nothing : "NCBIGene:" * entrez
    acc = !isempty(sp) ? sp[1] : !isempty(tr) ? tr[1] : nothing
    id = acc !== nothing ? protein_id(acc; isoforms) : gene !== nothing ? gene_product_id(gene) : nothing
    id === nothing && return nothing
    p = Protein(; id, name = cell(row("Official Symbol $s")), synonyms = split_list(row("Synonyms $s")),
                organism = taxon!(kg, row("Organism ID $s"); name = cell(row("Organism Name $s"))))
    # Other accessions are related, not equivalent: BioGRID lists the same TrEMBL entry
    # (e.g. a read-through product) under several genes, so they must not merge proteins.
    for a in [sp; tr]
        x = protein_id(a; isoforms)
        x == id || x in p.db_xrefs || push!(p.db_xrefs, x)
    end
    gene === nothing || push!(p.db_xrefs, gene)
    (b = cell(row("BioGRID ID $s"))) === nothing || push!(p.db_xrefs, "BioGRID:" * b)
    upsert!(kg, p)
    return id
end

function biogrid_gene!(kg, row, side)
    s = side == 'A' ? "Interactor A" : "Interactor B"
    entrez = cell(row("Entrez Gene $s"))
    entrez === nothing && return nothing
    g = Gene(id = "NCBIGene:" * entrez, name = cell(row("Official Symbol $s")),
             synonyms = split_list(row("Synonyms $s")),
             organism = taxon!(kg, row("Organism ID $s"); name = cell(row("Organism Name $s"))))
    (b = cell(row("BioGRID ID $s"))) === nothing || push!(g.db_xrefs, "BioGRID:" * b)
    upsert!(kg, g)
    return g.id
end

# Iterate over the rows of a tab-separated file with a header line; `f(row)` gets a
# function from column name to cell text.
function each_row(f, io::IO; header_prefix = "#")
    lines = eachline(io)
    header = String.(split(replace(first(lines), r"^#" => ""), '\t'))
    col = Dict(h => i for (i, h) in enumerate(header))
    for line in lines
        isempty(strip(line)) && continue
        fields = split(line, '\t')
        f(name -> (i = get(col, name, nothing); i === nothing || i > length(fields) ? nothing : fields[i]))
    end
end

"""
    read_biogrid(src; physical = true, genetic = false, organisms = nothing, isoforms = :merge,
                 strict = true) -> KnowledgeGraph

Read a BioGRID TAB3 file (`*.tab3.txt`, or a BioGRID `.zip` download; the first entry
with "INTERACTIONS" or "tab3" in its name is used).

- **Physical** interactions connect `Protein`s, named by UniProt accession (Swiss-Prot
  first) or, without one, `gene_product/NCBIGene:<id>`; the Entrez and BioGRID gene ids
  go to `db_xrefs`. The experimental system sets the type (`BIOGRID_SYSTEM`):
  direct methods (two-hybrid, co-crystal, ...) give `:Binding`, affinity capture and
  co-fractionation `:Association`, co-localization `:CoOccurence`. "Biochemical Activity"
  rows with a modification (phosphorylation, ubiquitination, ...) become directed
  `Influence`s from A to B.
- **Genetic** interactions (with `genetic = true`) connect `Gene`s (`NCBIGene:<id>`) as
  `:Association`.
- The score, when present, becomes `weight` (BioGRID scores are method-specific, not
  probabilities); the PubMed id becomes `publications`. Other UniProt accessions of the
  gene go to `db_xrefs`.
"""
function read_biogrid(src; physical::Bool = true, genetic::Bool = false, organisms = nothing,
                      isoforms::Symbol = :merge, strict::Bool = true)
    kg = KnowledgeGraph()
    skipped = Skipped()
    keep_tax = organisms === nothing ? nothing : Set(string.(organisms))
    member = src isa AbstractString && endswith(lowercase(src), ".zip") ?
             something(zip_member(src, r"INTERACTIONS|tab3"i), nothing) : nothing
    with_source(src; member) do io
        each_row(io) do row
            kind = lowercase(something(cell(row("Experimental System Type")), ""))
            if keep_tax !== nothing && !(strip(row("Organism ID Interactor A")) in keep_tax &&
                                         strip(row("Organism ID Interactor B")) in keep_tax)
                return skip!(skipped, "interaction outside `organisms`")
            end
            system = lowercase(something(cell(row("Experimental System")), ""))
            id = "biogrid_interaction/" * strip(row("BioGRID Interaction ID"))
            pubs = String[p for p in (pmid(something(row("Publication Source"), "")),) if p !== nothing]
            # BioGRID scores are method-specific (MS scores, SAINT, ...), not probabilities:
            # they go to `weight`, not `confidence`.
            score = tryparse(Float64, something(cell(row("Score")), ""))
            if kind == "physical" && physical
                a, b = biogrid_protein!(kg, row, 'A', isoforms), biogrid_protein!(kg, row, 'B', isoforms)
                (a === nothing || b === nothing) && return skip!(skipped, "physical interaction without ids")
                mod = lowercase(something(cell(row("Modification")), ""))
                type = system == "biochemical activity" ? get(BIOGRID_MODIFICATION, mod, :StateTransition) :
                       get(BIOGRID_SYSTEM, system, :Association)
                edge!(kg, type, [a], [b]; id, source = "biogrid", strict, weight = score, publications = pubs)
            elseif kind == "genetic" && genetic
                a, b = biogrid_gene!(kg, row, 'A'), biogrid_gene!(kg, row, 'B')
                (a === nothing || b === nothing) && return skip!(skipped, "genetic interaction without Entrez ids")
                edge!(kg, :Association, [a], [b]; id, source = "biogrid", strict, weight = score, publications = pubs)
            else
                skip!(skipped, "$kind interaction (not requested)")
            end
        end
    end
    report(skipped, "BioGRID")
    stamp!(kg, source_provenance!(kg, "BioGRID"; url = "https://thebiogrid.org", file = src isa IO ? nothing : src))
    return kg
end

function zip_member(path, pattern)
    z = ZipFile.Reader(path)
    try
        i = findfirst(e -> occursin(pattern, e.name), z.files)
        return i === nothing ? nothing : z.files[i].name
    finally
        close(z)
    end
end

# BioGRID chemical actions => relationship type (chemical → protein).
const BIOGRID_ACTION = Dict(
    "inhibitor" => :Inhibition, "antagonist" => :Inhibition, "inverse agonist" => :Inhibition,
    "blocker" => :Inhibition, "negative modulator" => :Inhibition, "suppressor" => :Inhibition,
    "agonist" => :Activation, "activator" => :Activation, "positive modulator" => :Activation,
    "inducer" => :Activation, "potentiator" => :Activation, "partial agonist" => :Activation,
    "binder" => :Binding, "binding" => :Binding,
)

"""
    read_biogrid_chemicals(src; organisms = nothing, isoforms = :merge) -> KnowledgeGraph

Read a BioGRID chemical-association file (`*.chemtab.txt`, or the zip containing it).
Each chemical becomes a `SmallMolecule` named by its source database id (DrugBank,
...), with name, synonyms, formula, CAS number and InChIKey; each association an edge
from the chemical to the target protein, typed by its action (`BIOGRID_ACTION`:
inhibitors and antagonists `:Inhibition`, agonists and activators `:Activation`,
otherwise `:Association`).
"""
function read_biogrid_chemicals(src; organisms = nothing, isoforms::Symbol = :merge, strict::Bool = true)
    kg = KnowledgeGraph()
    keep_tax = organisms === nothing ? nothing : Set(string.(organisms))
    member = src isa AbstractString && endswith(lowercase(src), ".zip") ? zip_member(src, r"chemtab"i) : nothing
    with_source(src; member) do io
        each_row(io) do row
            keep_tax === nothing || strip(row("Organism ID")) in keep_tax || return
            entrez = cell(row("Entrez Gene ID"))
            entrez === nothing && return
            gene = "NCBIGene:" * entrez
            target = gene_product_id(gene)
            upsert!(kg, Protein(id = target, name = cell(row("Official Symbol")), synonyms = split_list(row("Synonyms")),
                                organism = taxon!(kg, row("Organism ID"); name = cell(row("Organism"))), db_xrefs = [gene]))
            srcdb, srcid = cell(row("Chemical Source")), cell(row("Chemical Source ID"))
            chem = srcdb !== nothing && srcid !== nothing ? curie(srcdb, srcid) : "BioGRID.chemical:" * strip(row("BioGRID Chemical ID"))
            upsert!(kg, SmallMolecule(id = chem, name = cell(row("Chemical Name")), synonyms = split_list(row("Chemical Synonyms")),
                                      formula = cell(row("Molecular Formula")), cas_number = cell(row("CAS Number")),
                                      inchikey = cell(row("InChIKey"))))
            action = lowercase(something(cell(row("Action")), ""))
            pubs = String[p for p in (pmid(something(row("Pubmed ID"), "")),) if p !== nothing]
            edge!(kg, get(BIOGRID_ACTION, action, :Association), [chem], [target];
                  id = "biogrid_chemical/" * strip(row("BioGRID Chemical Interaction ID")), source = "biogrid",
                  strict, publications = pubs)
        end
    end
    stamp!(kg, source_provenance!(kg, "BioGRID"; url = "https://thebiogrid.org", file = src isa IO ? nothing : src))
    return kg
end

"""
    fetch_biogrid(genes; accesskey, organism = 9606, cache = cache_dir(), refresh = false,
                  limit = 10_000) -> path

Interactions of the given gene symbols from the BioGRID REST service, as TAB3. Needs a
free access key from https://webservice.thebiogrid.org. For whole organisms, download a
release file (`BIOGRID-ORGANISM-…tab3.zip`) from https://downloads.thebiogrid.org instead
and pass it to `read_biogrid`.
"""
function fetch_biogrid(genes::AbstractVector{<:AbstractString}; accesskey::AbstractString, organism::Integer = 9606,
                       limit::Integer = 10_000, kwargs...)
    url = "https://webservice.thebiogrid.org/interactions/?searchNames=true&includeHeader=true&format=tab3" *
          "&geneList=$(escape_uri(join(genes, '|')))&taxId=$organism&max=$limit&accesskey=$(escape_uri(accesskey))"
    return fetch_url(url; name = "biogrid_" * string(hash((genes, organism, limit)); base = 16) * ".tab3.txt", kwargs...)
end
