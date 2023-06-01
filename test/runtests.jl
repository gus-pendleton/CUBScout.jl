using CUBScout
using Test

@testset "CUBScout.jl" begin
    @test round(b(example_data_path).self[1], digits = 10) == 0.2091269922
    @test round(enc(example_data_path).ENC[1], digits = 10) == 56.7872822025
    @test round(enc_p(example_data_path).self[1], digits = 10) == 61.0
    @test round(mcb(example_data_path).self[1], digits = 10) == 0.0872112376
    @test round(milc(example_data_path).self[1], digits = 10) == 0.494825732
    @test round(scuo(example_data_path).SCUO[1], digits = 10) == 0.1431209294
    ribosomal_genes = find_seqs(example_data_path, r"ribosomal")
    @test round(e(example_data_path, ribosomal_genes).E[1], digits = 10) == 0.7623168751
    @test round(melp(example_data_path, ribosomal_genes).MELP[1], digits = 10) == 0.9294138732
    @test round(cai(example_data_path, ribosomal_genes).CAI[1], digits = 10) == 0.8449667854
    @test round(fop(example_data_path, ribosomal_genes).FOP[1], digits = 10) == 0.567816092
    @test round(gcb(example_data_path).GCB[1], digits = 10) == -0.0587654329
end

