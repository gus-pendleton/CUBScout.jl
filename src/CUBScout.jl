module CUBScout

    using FASTX
    using DelimitedFiles
    using BioSequences
    using Artifacts
    using Statistics

    export CodonDict, make_CodonDict, CodonDict_PATH, DEFAULT_CodonDict, ALTSTART_CodonDict, EXAMPLE_DATA_PATH, find_seqs, seq_names, seq_descriptions
    export b, enc, enc_p, mcb, milc, scuo, all_cub
    export e, melp, cai, fop, gcb

    include("codon_dict.jl")
    include("accessory_functions.jl")
    
    include("CUB_measures.jl")
    include("EXP_measures.jl")


end
