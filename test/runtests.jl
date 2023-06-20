using CUBScout
using Test

@testset "CUBScout.jl" begin
    @test isapprox(b(EXAMPLE_DATA_PATH).self[1], 0.2091269922)
    @test isapprox(enc(EXAMPLE_DATA_PATH).ENC[1], 56.7872822025) 
    @test isapprox(enc_p(EXAMPLE_DATA_PATH).self[1], 61.0)
    @test isapprox(mcb(EXAMPLE_DATA_PATH).self[1], 0.0872112376)
    @test isapprox(milc(EXAMPLE_DATA_PATH).self[1], 0.494825732)
    ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal")
    @test isapprox(melp(EXAMPLE_DATA_PATH, ribosomal_genes).MELP[1], 0.9294138732)
    @test isapprox(cai(EXAMPLE_DATA_PATH, ribosomal_genes).CAI[1], 0.8449667854)
    @test isapprox(fop(EXAMPLE_DATA_PATH, ribosomal_genes).FOP[1], 0.567816092)
    @test isapprox(gcb(EXAMPLE_DATA_PATH).GCB[1], -0.0587654329)
    codon_counts = count_codons(EXAMPLE_DATA_PATH);
    @test isapprox(codon_frequency(codon_counts[1], "net_genomic")[1],  0.04941242971710299)
    @test isapprox(codon_frequency(codon_counts[1], "net_gene")[1],  0.07158836689038031)
    @test isapprox(codon_frequency(codon_counts[1], "byAA_gene")[1],  0.8421052631578947)
    @test isapprox(codon_frequency(codon_counts[1], "byAA_genomic")[1],  0.7016640025759265)

end
