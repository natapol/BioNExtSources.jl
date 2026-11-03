"""
    BioNExtSources

Populate BioNExt knowledge graphs from biological databases: UniProt, KEGG, IntAct
(and any PSICQUIC service), BioGRID, SIGNOR, Reactome, BioCyc and Pathway Commons.
Every reader returns an ordinary `KnowledgeGraph` that follows the BioNExt ontology,
ready to be combined with `populate!` / `integrate!`.
"""
module BioNExtSources

using BioNExt
using BioNExt: Resource, Entity, add!
using Dates, Downloads
using CodecZlib: GzipDecompressorStream
import EzXML
import ZipFile

include("common.jl")
include("uniprot.jl")
include("mitab.jl")
include("biogrid.jl")
include("signor.jl")
include("kegg.jl")
include("biopax.jl")
include("populate.jl")

export
    # shared
    curie, protein_id, gene_product_id, location_terms, residue_site, fit_type, edge!, upsert!,
    taxon!, source_provenance!, fetch_url, cache_dir, with_source,
    # UniProt
    UNIPROT_FIELDS, fetch_uniprot, read_uniprot,
    # PSI-MITAB, IntAct, PSICQUIC
    mi_values, mi_relationship, read_mitab, read_intact, fetch_psicquic, fetch_intact, psicquic_services,
    # BioGRID
    read_biogrid, read_biogrid_chemicals, fetch_biogrid,
    # SIGNOR
    read_signor, read_signor_complexes, fetch_signor,
    # KEGG
    read_kgml, read_kegg_conv, fetch_kegg, fetch_kgml, fetch_kegg_uniprot,
    # BioPAX: Reactome, BioCyc, Pathway Commons
    read_biopax, read_reactome, read_biocyc, fetch_reactome_biopax, fetch_pathwaycommons,
    # combining
    populate!, DEFAULT_PREFER, keep_best

end
