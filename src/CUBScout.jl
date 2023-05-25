module CUBScout

    using DataFrames
    using FASTX
    using DelimitedFiles
    using BioSequences
    using Artifacts

    export codon_dict, make_codon_dict, codon_dict_path, default_codon_dict, altstart_codon_dict, example_data_path, find_seqs, seq_names, seq_descriptions
    export b, enc, enc_p, mcb, milc, scuo
    export e, melp, cai, fop, gcb

    include("codon_dict.jl")
    include("accessory_functions.jl")
    
    include("CUB_measures.jl")
    include("EXP_measures.jl")


end
