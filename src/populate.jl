# Combining sources into one graph.

"""
Default identifier preference when sources are combined: the id of a merged entity is the
one whose prefix comes first here (UniProt for proteins, ChEBI for small molecules, ...).
"""
const DEFAULT_PREFER = ["UniProt", "CHEBI", "NCBIGene", "HGNC", "Ensembl", "ComplexPortal", "Reactome",
                        "kegg.pathway", "kegg.reaction", "kegg.compound", "IntAct", "SIGNOR", "pubchem.compound",
                        "DrugBank", "kegg.genes"]

"""
    keep_best(field, kept, other)

Conflict policy used by `populate!`: the higher `confidence` / `weight`, the longer
`sequence`, otherwise the value already in the graph.
"""
keep_best(field, kept, other) =
    field in (:confidence, :weight) && kept isa Number && other isa Number ? max(kept, other) :
    field === :sequence && kept isa AbstractString && other isa AbstractString ?
        (length(other) > length(kept) ? other : kept) : kept

"""
    populate!(kg, graphs...; prefer = DEFAULT_PREFER, on_conflict = keep_best, kwargs...) -> Resolution

Integrate graphs produced by the readers (`read_uniprot`, `read_kgml`, `read_mitab`, ...) into
`kg` one after another with `integrate!`, so the same molecule from different databases
becomes one entity and the same relationship one edge (with the publications, provenance
and best confidence of all its sources). Returns the combined `Resolution`; other keywords
go to `integrate!`.
"""
function populate!(kg::KnowledgeGraph, graphs::KnowledgeGraph...; prefer = DEFAULT_PREFER,
                   on_conflict = keep_best, kwargs...)
    mapping, groups, conflicts = Dict{String,String}(), Vector{String}[], Conflict[]
    for g in graphs
        res = integrate!(kg, g; prefer, on_conflict, kwargs...)
        for (k, v) in mapping
            mapping[k] = get(res.mapping, v, v)
        end
        merge!(mapping, res.mapping)
        append!(groups, res.groups)
        append!(conflicts, res.conflicts)
    end
    return Resolution(mapping, groups, conflicts)
end
