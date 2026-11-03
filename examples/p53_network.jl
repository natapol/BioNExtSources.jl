# Build a p53 signalling network from seven databases' worth of data.
#   julia --project=. examples/p53_network.jl
#
# Downloads go to ~/.cache/bionext (or $BIONEXT_CACHE) and are reused on later runs.
# SIGNOR and BioGRID are read from local files here: SIGNOR's server is often slow to
# reach, and BioGRID's REST service needs a (free) access key; see the README.

using BioNExt, BioNExtSources

seeds = ["P04637", "Q00987", "Q13315", "P38936"]           # p53, MDM2, ATM, p21
fixtures(x) = joinpath(@__DIR__, "..", "test", "fixtures", x)

uniprot  = read_uniprot(fetch_uniprot(join(("accession:" * a for a in seeds), " OR ")))
intact   = read_intact(fetch_intact("P04637 AND species:9606"; limit = 500); organisms = [9606])
kegg     = read_kgml(fetch_kgml("hsa04115"); uniprot = fetch_kegg_uniprot("hsa"), taxon = 9606)
reactome = read_reactome(fetch_reactome_biopax("R-HSA-69541"))
signor   = read_signor(fixtures("signor_handmade.tsv"))                # or read_signor(fetch_signor())
biogrid  = read_biogrid(fixtures("biogrid_sample.tab3.txt"))           # or a BioGRID release file

kg = KnowledgeGraph()
res = populate!(kg, uniprot, intact, kegg, reactome, signor, biogrid)
index!(kg)

show(stdout, MIME"text/plain"(), kg); println()
println(res)
println("issues: ", length(validate(kg)))

p53 = kg["UniProt:P04637"]
println("\np53 is described by: ", join(p53.provenance, ", "))
println("p53 takes part in ", degree_of(kg, p53), " edges; its biggest neighbourhood component has ",
        length(first(network_components(subgraph(kg, p53; hops = 1)))), " members")
println("ATM → MDM2: ", shortest_path(kg, "UniProt:Q13315", "UniProt:Q00987"; direction = :both))

mkpath(joinpath(@__DIR__, "output"))
save(joinpath(@__DIR__, "output", "p53_network.ttl"), kg)
save(joinpath(@__DIR__, "output", "p53_network.graphml"), kg; reify = true)
println("\nwrote examples/output/p53_network.{ttl,graphml}")
