module CUBScout

using FASTX: FASTAReader, sequence, identifier, description
using DelimitedFiles: readdlm
using BioSequences: BioSequences, LongDNA, NucSeq, encoded_data, @dna_str, @rna_str
using Artifacts
using Statistics: mean

export CodonDict,
    make_CodonDict,
    CodonDict_PATH,
    DEFAULT_CodonDict,
    ALTSTART_CodonDict,
    EXAMPLE_DATA_PATH,
    find_seqs,
    seq_names,
    seq_descriptions,
    count_codons,
    codon_frequency
export b, enc, enc_p, mcb, milc, scuo, all_cub
export e, melp, cai, fop, gcb

include("codon_dict.jl")
include("accessory_functions.jl")

include("CUB_measures.jl")
include("EXP_measures.jl")


end
