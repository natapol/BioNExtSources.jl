using BioNExt, BioNExtSources, Test

const F(x) = joinpath(@__DIR__, "fixtures", x)
typecounts(kg) = (d = Dict{Symbol,Int}(); foreach(e -> d[e.type] = get(d, e.type, 0) + 1, edges(kg)); d)

@testset "BioNExtSources" begin
    @testset "shared helpers" begin
        @test curie("uniprotkb", "P04637") == "UniProt:P04637"
        @test curie("ChEBI", "CHEBI:15422") == "CHEBI:15422"
        @test curie("entrez gene/locuslink", "7157") == "NCBIGene:7157"
        @test curie("SomeDB", "x1") == "SomeDB:x1"
        @test protein_id("P04637-2") == "UniProt:P04637"
        @test protein_id("uniprotkb:P04637-2"; isoforms = :keep) == "UniProt:P04637-2"
        @test BioNExtSources.isoform_parent("UniProt:P04637-2") == "UniProt:P04637"
        @test gene_product_id("NCBIGene:7157") == "gene_product/NCBIGene:7157"
        @test expand(KnowledgeGraph(), gene_product_id("NCBIGene:7157")) == BioNExt.BNR * "gene_product/NCBIGene:7157"
        @test location_terms("GO:0005654") == [:Nucleus]
        @test location_terms("GO:0005634") == [:Nucleus]
        @test location_terms("nucleoplasm") == [:Nucleus]
        @test location_terms("SUBCELLULAR LOCATION: Nucleus, nucleoplasm {ECO:1}. Cytoplasm {ECO:2}. Note=x") ==
              [:Nucleus, :Cytoplasm]
        @test location_terms("endoplasmic reticulum") == Symbol[]       # ambiguous: rough or smooth
        @test residue_site('S', 15) == "Ser15"
        @test residue_site("O-phospho-L-serine", 15) == "Ser15"
        @test residue_site("N6-acetyllysine", 120) == "Lys120"
        @test residue_site("Phosphoisoleucine?", 3) == "Ile3"        # longest name wins over "leucine"

        kg = KnowledgeGraph()
        add!(kg, Protein(id = "P"), Protein(id = "Q"), Gene(id = "G"), MolecularComplex(id = "C"))
        @test fit_type(kg, :Ubiquitination, ["P"], ["Q"]) === :Ubiquitination
        @test fit_type(kg, :Ubiquitination, ["C"], ["Q"]) === :StateTransition
        @test fit_type(kg, :RegulateExpression, ["P"], ["G"]) === :RegulateExpression
        @test fit_type(kg, :RegulateExpression, ["P"], ["Q"]) === :Increase
        @test fit_type(kg, :RegulateExpression, ["P"], ["Q"]; strict = false) === :RegulateExpression
        e = edge!(kg, :Binding, ["P"], ["Q"]; source = "test", confidence = 0.5)
        @test e isa Interaction && e.id == "test/Binding/P/Q" && e.confidence == 0.5
        @test edge!(kg, :Activation, ["P"], ["Q"]; source = "test") isa Influence

        vs = mi_values("psi-mi:\"MI:0915\"(physical association)|uniprotkb:TP53(gene name)")
        @test vs[1].db == "psi-mi" && vs[1].id == "MI:0915" && vs[1].text == "physical association"
        @test vs[2].id == "TP53" && vs[2].text == "gene name"
        @test mi_relationship(mi_values("psi-mi:\"MI:0407\"(direct interaction)")) == (:Binding, false)
        @test mi_relationship(mi_values("psi-mi:\"MI:0217\"(phosphorylation reaction)")) == (:Phosphorylation, true)
        @test mi_relationship(mi_values("psi-mi:\"MI:0192\"(acetylation reaction)")) == (:StateTransition, true)
        @test mi_relationship(BioNExtSources.MIValue[]; causal = mi_values("psi-mi:\"MI:2241\"(down-regulates activity)")) == (:Inhibition, true)
    end

    @testset "UniProt" begin
        kg = read_uniprot(F("uniprot_p53_mdm2.tsv"))
        p = kg["UniProt:P04637"]
        @test p.name == "Cellular tumor antigen p53"
        @test :Nucleus in p.cellular_location && :Cytoplasm in p.cellular_location
        @test "Ser15" in p.modification_sites && "Ser20" in p.modification_sites
        @test p.cofactors == ["CHEBI:29105"] && kg["CHEBI:29105"].name == "Zn(2+)"
        @test "STRING:9606.ENSP00000269305" in p.external_ids
        @test p.organism == "NCBITaxon:9606" && kg["NCBITaxon:9606"].name == "Homo sapiens"
        g = kg["NCBIGene:7157"]
        @test g.name == "TP53" && g.external_ids == ["HGNC:11998"] && "GO:0003700" in g.functions
        @test [(e.type, e.left, e.right) for e in edges_of(kg, p)] == [(:Expression, ["NCBIGene:7157"], ["UniProt:P04637"])]
        @test p.provenance == ["prov/uniprot"]
        @test isvalid(kg; strict = true)
        @test !any(r -> r isa Gene, read_uniprot(F("uniprot_p53_mdm2.tsv"); genes = false))
    end

    @testset "PSI-MITAB / IntAct" begin
        kg = read_intact(F("intact_p53.mitab27.txt"))
        @test length(edges(kg)) == 20
        # IntAct records keep their id; spoke-expanded pairs get a minted "intact/EBI-…/A+B" id
        @test all(e -> startswith(e.id, "IntAct:EBI-") || startswith(e.id, "intact/EBI-"), edges(kg))
        p = kg["UniProt:P04637"]
        @test p.name == "TP53" && "IntAct:EBI-366083" in p.external_ids && p.organism == "NCBITaxon:9606"
        e = kg["IntAct:EBI-11297132"]
        @test e isa Interaction && e.type == :Association && e.confidence == 0.98 && e.publications == ["PMID:25402006"]
        @test isvalid(kg; strict = true)
        human = read_intact(F("intact_p53.mitab27.txt"); organisms = [9606])
        @test all(n -> n.organism == "NCBITaxon:9606", nodes(human))
        @test length(edges(human)) <= 20
    end

    @testset "BioGRID" begin
        kg = read_biogrid(F("biogrid_sample.tab3.txt"))
        e = kg["biogrid_interaction/418"]
        @test e.type == :Binding && e.left == ["UniProt:P15927"] && e.right == ["UniProt:P40763"]
        @test e.publications == ["PMID:10875894"]
        p = kg["UniProt:P15927"]
        @test p.name == "RPA2" && isempty(p.external_ids) && "UniProt:B4DUL2" in p.db_xrefs && "NCBIGene:6118" in p.db_xrefs
        @test "RPA32" in p.synonyms
        @test Set(keys(typecounts(kg))) ⊆ Set([:Binding, :Association, :CoOccurence])
        @test isvalid(kg; strict = true)

        c = read_biogrid_chemicals(F("biogrid_sample.chemtab.txt"))
        x = c["biogrid_chemical/1912"]
        @test x.left == ["DrugBank:DB00126"] && x.right == ["gene_product/NCBIGene:112399"]
        m = c["DrugBank:DB00126"]
        @test m.name == "Vitamin C" && m.formula == "C6H8O6" && m.cas_number == "50-81-7"
        @test m.inchikey == "CIWBSHSKHKDKBQ-JLAZNSOCSA-N"
        @test isvalid(c; strict = true)
    end

    @testset "SIGNOR" begin
        cx = read_signor_complexes(F("signor_handmade_complexes.csv"))
        @test cx["SIGNOR:SIGNOR-C1"] == ("MRN complex", ["UniProt:P49959", "UniProt:Q92878", "UniProt:O60934"])
        kg = read_signor(F("signor_handmade.tsv"); complexes = cx)
        @test [(e.id, e.type) for e in edges(kg)] == [
            ("SIGNOR:SIGNOR-TEST-1", :Phosphorylation), ("SIGNOR:SIGNOR-TEST-2", :Ubiquitination),
            ("SIGNOR:SIGNOR-TEST-3", :Increase), ("SIGNOR:SIGNOR-TEST-4", :Inhibition),
            ("SIGNOR:SIGNOR-TEST-5", :Activation)]                 # the stimulus row is skipped
        @test kg["SIGNOR:SIGNOR-TEST-1"].confidence == 0.9 && kg["SIGNOR:SIGNOR-TEST-1"].publications == ["PMID:9733515"]
        @test kg["UniProt:P04637"].modification_sites == ["Ser15"]
        @test kg["SIGNOR:SIGNOR-C1"].subunits == ["UniProt:P49959", "UniProt:Q92878", "UniProt:O60934"]
        @test kg["CHEBI:46742"] isa SmallMolecule
        @test isvalid(kg; strict = true)
        noheader = read_signor(F("signor_handmade_noheader.tsv"); complexes = cx)
        @test [e.id for e in edges(noheader)] == [e.id for e in edges(kg)]
        loose = read_signor(F("signor_handmade.tsv"); strict = false)
        @test loose["SIGNOR:SIGNOR-TEST-3"].type == :RegulateExpression
    end

    @testset "KEGG" begin
        conv = read_kegg_conv(F("kegg_conv_uniprot.txt"))
        @test conv["hsa:7157"] == ["P04637", "K7PPA8", "Q53GA5"]

        kg = read_kgml(F("kegg_hsa04115.kgml"); taxon = 9606)
        pw = kg["kegg.pathway:hsa04115"]
        @test pw.name == "p53 signaling pathway"
        @test haskey(kg, "gene_product/kegg.genes:hsa:7157")
        @test typecounts(kg)[:Increase] > 0 && !haskey(typecounts(kg), :RegulateExpression)   # no genes
        @test all(id -> haskey(kg, id), pw.components)
        @test length(resources(kg, MolecularComplex)) == 3
        @test isvalid(kg; strict = true)

        g = read_kgml(F("kegg_hsa04115.kgml"); uniprot = conv, genes = true, taxon = 9606)
        p53 = g["UniProt:P04637"]
        @test p53.name == "TP53" && p53.external_ids == ["UniProt:K7PPA8", "UniProt:Q53GA5"]
        @test "kegg.genes:hsa:7157" in p53.db_xrefs
        reg = edges_of(g, "kegg.genes:hsa:4193"; side = :right, type = :RegulateExpression)
        @test any(e -> e.left == ["UniProt:P04637"], reg)                 # p53 drives MDM2 expression
        @test isvalid(g; strict = true)

        m = read_kgml(F("kegg_hsa00010.kgml"))
        r = m["kegg.reaction:R00710"]
        @test r.type == :Reversible && r.left == ["kegg.compound:C00084"] && r.right == ["kegg.compound:C00033"]
        @test length(edges_of(m, r; side = :right, type = :Catalysis)) == 5
        @test typecounts(m)[:SuccessiveReaction] > 0
        @test isvalid(m; strict = true)
    end

    @testset "BioPAX / Reactome" begin
        kg = read_reactome(F("reactome_R-HSA-69541.owl.gz"))
        @test [p.id for p in pathways(kg)] == ["Reactome:R-HSA-69541", "Reactome:R-HSA-349425"]
        @test kg["Reactome:R-HSA-69541"].name == "Stabilization of p53"
        rx = kg["Reactome:R-HSA-5693609"]
        @test rx.name == "ATM phosphorylates TP53 at S15" && rx.type == :Irreversible
        @test rx.left == ["Reactome:R-HSA-3209194", "CHEBI:30616"] && rx.right == ["Reactome:R-HSA-349474", "CHEBI:456216"]
        @test "EC:2.7.11.1" in rx.ec_numbers && "PMID:9733514" in rx.publications
        cat = only(edges_of(kg, rx; side = :right))
        @test cat.type == :Catalysis && cat.left == ["UniProt:Q13315"]
        @test kg["UniProt:Q13315"].name == "ATM"
        p = kg["UniProt:P04637"]
        @test p.name == "TP53" && :Nucleus in p.cellular_location && "Ser15" in p.modification_sites
        @test kg["Reactome:R-HSA-3209194"].subunits == ["UniProt:P04637"]
        @test isvalid(kg; strict = true)

        keep = read_reactome(F("reactome_R-HSA-69541.owl.gz"); states = :keep)
        forms = filter(n -> n isa Protein && "UniProt:P04637" in n.db_xrefs, nodes(keep))
        @test length(forms) > 1 && !haskey(keep, "UniProt:P04637")
        @test isvalid(keep; strict = true)
    end

    @testset "populate!" begin
        kg = KnowledgeGraph()
        res = populate!(kg,
            read_uniprot(F("uniprot_p53_mdm2.tsv")),
            read_intact(F("intact_p53.mitab27.txt")),
            read_biogrid(F("biogrid_sample.tab3.txt")),
            read_kgml(F("kegg_hsa04115.kgml"); uniprot = read_kegg_conv(F("kegg_conv_uniprot.txt")), taxon = 9606),
            read_reactome(F("reactome_R-HSA-69541.owl.gz")),
            read_signor(F("signor_handmade.tsv"); complexes = read_signor_complexes(F("signor_handmade_complexes.csv"))))
        p = kg["UniProt:P04637"]
        @test p.provenance == ["prov/uniprot", "prov/intact", "prov/kegg", "prov/reactome", "prov/signor"]
        @test p.name == "Cellular tumor antigen p53"                    # UniProt came first
        @test :Nucleus in p.cellular_location
        @test count(n -> n isa Protein && occursin("P04637", n.id), nodes(kg)) == 1
        @test isempty(duplicates(kg))
        @test isvalid(kg; strict = true)
        @test length(resources(kg, Provenance)) == 6
        @test keep_best(:confidence, 0.5, 0.9) == 0.9 && keep_best(:name, "a", "b") == "a"
    end

    @testset "imports pass the data-quality audit" begin
        graphs = [
            "UniProt" => read_uniprot(F("uniprot_p53_mdm2.tsv")),
            "IntAct" => read_intact(F("intact_p53.mitab27.txt")),
            "BioGRID" => read_biogrid(F("biogrid_sample.tab3.txt")),
            "BioGRID chemicals" => read_biogrid_chemicals(F("biogrid_sample.chemtab.txt")),
            "SIGNOR" => read_signor(F("signor_handmade.tsv"); complexes = read_signor_complexes(F("signor_handmade_complexes.csv"))),
            "KEGG" => read_kgml(F("kegg_hsa04115.kgml"); uniprot = read_kegg_conv(F("kegg_conv_uniprot.txt")), genes = true),
            "KEGG metabolism" => read_kgml(F("kegg_hsa00010.kgml")),
            "Reactome" => read_reactome(F("reactome_R-HSA-69541.owl.gz")),
        ]
        for (name, g) in graphs
            r = audit(g)
            count_findings(r; severity = :error) == 0 || @info "$name audit" r
            @test passed(r)
            # identifiers the readers mint must look right, except the hand-written SIGNOR ids
            name == "SIGNOR" || @test count_findings(r; check = :identifiers) == 0
        end
        signor = audit(last(graphs[5]); checks = [:identifiers])
        @test count_findings(signor) == 5                        # SIGNOR-TEST-n are not real SIGNOR ids
    end

    if get(ENV, "BIONEXT_ONLINE_TESTS", "false") == "true"
        @testset "online" begin
            cache = mktempdir()
            u = read_uniprot(fetch_uniprot("accession:P04637"; cache))
            @test u["UniProt:P04637"].name == "Cellular tumor antigen p53"
            @test length(edges(read_kgml(fetch_kgml("hsa04115"; cache)))) > 50
            @test fetch_kegg_uniprot("hsa"; cache)["hsa:7157"][1] == "P04637"
            @test haskey(psicquic_services(; cache), "IntAct")
            i = read_intact(fetch_intact("P04637"; page = 50, limit = 120, cache))
            @test length(edges(i)) == 120
            @test pathways(read_reactome(fetch_reactome_biopax("R-HSA-69541"; cache)))[1].name == "Stabilization of p53"
        end
    else
        @info "Skipping online tests (set BIONEXT_ONLINE_TESTS=true to run them)"
    end
end
