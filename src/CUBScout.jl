module CUBScout

    using Dataframes
    using FASTX
    using DelimitedFiles
    using BioSequences

    export codon_dict, make_codon_dict, find_seqs, seq_names, seq_descriptions
    export b, enc, enc_p, mcb, milc, scuo
    export e, melp, cai, fop, gcb

    include("accessory_functions.jl")
    include("codon_dict.jl")
    include("CUB_measures.jl")
    include("EXP_measures.jl")


end
