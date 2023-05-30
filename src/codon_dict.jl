"""
    codon_dict
The `codon_dict` type defines how codons are translated, and is essential for calculating codon usage bias as it identifies stop codons and each amino acid's degeneracy. A default codon dictionary is provided (`default_codon_dict`), or a user can make their own using the `make_codon_dict` function.

### Fields

- `codons`: the 64 codons, in alphabetical order
- `AA`: corresponding amino acid for each codon (64 entries long)
- `AA_nostops`: same as AA, but with stop codons removed
- `uniqueAA`: unique amino acid names including stop codons. Under a standard translation table, this is 21 amino acids long
- `uniqueAA`: same as uniqueAA, but with stop codons removed
- `uniqueI`: a vector of the same length as uniqueAA, containing vectors of the indices of each codon for that amino acid. For instance, the first entry corresponds to Lysine, and contains the vector `[1, 3]`, corresponding to the positions of codons AAA and AAG in the codons field
- `uniqueI_nostops`: same as uniqueI, but with stop codons removed
- `deg`: a vector of the same length as uniqueAA, containing the degeneracy for each amino acid.
- `deg_nostops`: same as deg, but with stop codons removed
- `stop_mask`: a Boolean vector of length 64 which is false for stop codons. This is used to remove stop codons when calculating codon usage bias.

### Notes

Generally, CUBScout users shouldn't need to interact with the `codon_dict` type, as the standard genetic code is applied by default. Details for constructing a custom `codon_dict` are documented under the `make_codon_dict` function.
"""
struct codon_dict
    codons::Vector{LongDNA{2}} # 64 codons (in alphabetical order)
    AA::Vector{String} # Amino acids corresponding to each codon
    AA_nostops::Vector{String} # Amino acids without stop codons
    uniqueAA::Vector{String} # Unique amino acids (21, generally)
    uniqueAA_nostops::Vector{String} # Amino acids, stops removed
    uniqueI::Vector{Vector{Int32}} # Indices of codons for each amino acid
    uniqueI_nostops::Vector{Vector{Int32}} # Indices after stops are removed
    deg::Vector{Int32} # Degeneracy of each amino acid
    deg_nostops::Vector{Int32} # Degeneracy with each amino acid removed
    stop_mask::Vector{Bool} # Boolean vector length 64 which removes stop codons
end

"""
    make_codon_dict(filepath::AbstractString, delimiter::AbstractChar = '\t')
Make a custom codon dictionary for organisms with non-standard genetic code. `filepath` points to a delimited file with two columns and no header. the first column should be codons, and the second column their corresponding amino acid. Avoid spaces and special characters (e.g., write GlutamicAcid instead of Glutamic Acid). Stop codons can be coded as Stop, stop, STOP, or *. If delimited using any character outside of tab, supply the delimiter as the second argument as Char, not a string (e.g. ',' not ","). `make_codon_dict` uses `readdlm` from `DelimitedFiles`; it's a good idea to check whether `readdlm` parses your file correctly before passing to `make_codon_dict`

# Examples
```jldoctest
julia> my_codon_dict = make_codon_dict(CUBScout.codon_dict_path)
codon_dict(BioSequences.LongSequence{BioSequences.DNAAlphabet{2}}[AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC  …  TCG, TCT, TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT], ["Lysine", "Asparagine", "Lysine", "Asparagine", "Threonine", "Threonine", "Threonine", "Threonine", "Arginine", "Serine"  …  "Serine", "Serine", "Stop", "Cysteine", "Tryptophan", "Cysteine", "Leucine", "Phenylalanine", "Leucine", "Phenylalanine"], ["Lysine", "Asparagine", "Lysine", "Asparagine", "Threonine", "Threonine", "Threonine", "Threonine", "Arginine", "Serine"  …  "Serine", "Serine", "Serine", "Cysteine", "Tryptophan", "Cysteine", "Leucine", "Phenylalanine", "Leucine", "Phenylalanine"], ["Lysine", "Asparagine", "Threonine", "Arginine", "Serine", "Isoleucine", "Methionine", "Glutamine", "Histidine", "Proline"  …  "Glutamicacid", "Asparticacid", "Alanine", "Glycine", "Valine", "Stop", "Tyrosine", "Cysteine", "Tryptophan", "Phenylalanine"], ["Lysine", "Asparagine", "Threonine", "Arginine", "Serine", "Isoleucine", "Methionine", "Glutamine", "Histidine", "Proline", "Leucine", "Glutamicacid", "Asparticacid", "Alanine", "Glycine", "Valine", "Tyrosine", "Cysteine", "Tryptophan", "Phenylalanine"], Vector{Int32}[[1, 3], [2, 4], [5, 6, 7, 8], [9, 11, 25, 26, 27, 28], [10, 12, 53, 54, 55, 56], [13, 14, 16], [15], [17, 19], [18, 20], [21, 22, 23, 24]  …  [33, 35], [34, 36], [37, 38, 39, 40], [41, 42, 43, 44], [45, 46, 47, 48], [49, 51, 57], [50, 52], [58, 60], [59], [62, 64]], Vector{Int32}[[1, 3], [2, 4], [5, 6, 7, 8], [9, 11, 25, 26, 27, 28], [10, 12, 51, 52, 53, 54], [13, 14, 16], [15], [17, 19], [18, 20], [21, 22, 23, 24], [29, 30, 31, 32, 58, 60], [33, 35], [34, 36], [37, 38, 39, 40], [41, 42, 43, 44], [45, 46, 47, 48], [49, 50], [55, 57], [56], [59, 61]], Int32[2, 2, 4, 6, 6, 3, 1, 2, 2, 4  …  2, 2, 4, 4, 4, 3, 2, 2, 1, 2], Int32[2, 2, 4, 6, 6, 3, 1, 2, 2, 4, 6, 2, 2, 4, 4, 4, 2, 2, 1, 2], Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  1, 1, 0, 1, 1, 1, 1, 1, 1, 1])

julia> typeof(my_codon_dict)
codon_dict

julia> fieldnames(codon_dict)
(:codons, :AA, :AA_nostops, :uniqueAA, :uniqueAA_nostops, :uniqueI, :uniqueI_nostops, :deg, :deg_nostops, :stop_mask)
```
"""
function make_codon_dict(filepath::AbstractString, delimiter::AbstractChar = '\t')
    cod_mat = readdlm(filepath, delimiter)
    alph_cod_mat = sortslices(cod_mat,dims=1,by=x->x[1],rev=false) # Alphabetize codons
    uniqueAA = unique(alph_cod_mat[:,2]) # Find unique amino acid names
    stop = map(x->occursin(r"Stop|stop|STOP|\*", x), alph_cod_mat[:,2]) # Search for stop codons
    AA_nostops = alph_cod_mat[.!stop,2] # Filter out stops
    uniqueAA_nostops = unique(AA_nostops) # Filter out stops
    indices = [findall(alph_cod_mat[:,2] .== aa) for aa in uniqueAA] # Find codon indices for each 
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

const codon_dict_path = joinpath(artifact"codon_dict", "codon_dict.txt")

const default_codon_dict = make_codon_dict(codon_dict_path)

const altstart_codon_dict = make_codon_dict(joinpath(artifact"codon_dict_altstart", "codon_dict_altstart.txt"))

"""
    example_data_path
The path to an example dataset, stored as an artifact within the package. This is an .fna file containing coding sequences from Bacillus subtilis subsp. subtilis str. 168, NCBI Accession # NC_000964.3.
"""
const example_data_path = joinpath(artifact"example_genome", "B_subtilis.fna")
