struct codon_dict
    codons::Vector{LongDNA{2}} # 64 codons (in alphabetical order)
    AA::Vector{String} # Amino acids corresponding to each codon
    AA_nostops::Vector{String} # Amino acids withou stop codons
    uniqueAA::Vector{String} # Unique amino acids (21, generally)
    uniqueAA_nostops::Vector{String} # Amino acids, stops removed
    uniqueI::Vector{Vector{Int32}} # Indices of codons for each amino acid
    uniqueI_nostops::Vector{Vector{Int32}} # Indices after stops are removed
    deg::Vector{Int32} # Degeneracy of each amino acid
    deg_nostops::Vector{Int32} # Degeneracy with each amino acid removed
    stop_mask::Vector{Bool} # Boolean vector length 64 which removes stop codons
end


function make_codon_dict(filepath::AbstractString, delimiter::AbstractChar = '\t')
    cod_mat = readdlm(filepath, delimiter)
    alph_cod_mat = sortslices(cod_mat,dims=1,by=x->x[1],rev=false) # Alphabetize codons
    uniqueAA = unique(alph_cod_mat[:,2])
    stop = map(x->occursin(r"Stop|stop|STOP|\*", x), alph_cod_mat[:,2]) # Search for stop codons
    AA_nostops = alph_cod_mat[.!stop,2]
    uniqueAA_nostops = unique(AA_nostops)
    indices = [findall(alph_cod_mat[:,2] .== aa) for aa in uniqueAA]
    indices_nostops = [findall(AA_nostops .== aa) for aa in uniqueAA_nostops]
    codon_dict(
        LongDNA{2}.(alph_cod_mat[:,1]), 
        alph_cod_mat[:,2], 
        AA_nostops,
        uniqueAA, 
        uniqueAA_nostops,
        indices, 
        indices_nostops,
        length.(indices),
        length.(indices_nostops),
        .!stop
        )
end

default_dict = make_codon_dict(pwd()*"/Data_Artifacts/codon_dict.txt")
alt_dict = make_codon_dict(pwd()*"/codon_dict_altstart.txt")
