# B
"""
    b(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate B from Karlin and Mrazek, 1996. 

# Arguments
-`filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
-`dict`: codon dictionary of type codon_dict. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
-`ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, ref_seqs takes a named tuple in the form `("subset_name" = Bool[],)``, where `Bool[]`` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset.
-`rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` and `false`.
-`rm_stop`: whether to remove stop codons from calculations of codon usage bias.
-`threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold`` to 0.
-`dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
julia> result = b(example_data_path) # Run B on example dataset
(self = [0.20912699220973896, 0.3289759448740455, 0.22365336363593893, 0.5391135258658497, 0.24919594143501034, 0.2880358413249049, 0.31200964304415874, 0.34858035204347476, 0.2455189361074733, 0.4690734561271221  …  0.3629137353834403, 0.3621330537227321, 0.4535285720373026, 0.3357858047622507, 0.28183191395624935, 0.2668809561422238, 0.22381338105820905, 0.4034837015709619, 0.3594626865160133, 0.3724863965444541],)

julia> result_300 = b(example_data_path, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> b(example_data_path, altstart_codon_dict) # Use alternative start codons
(self = [0.20897234061622738, 0.33515000264964157, 0.23659038285006437, 0.5444798345895256, 0.2510726777670733, 0.2931440715811394, 0.32097661134289895, 0.35705614480228676, 0.25452296343822073, 0.513313870450466  …  0.34414371567411556, 0.38229515825882665, 0.4592524704597901, 0.3399477982926337, 0.29297757306048133, 0.2680028918895221, 0.2168486105068708, 0.414543030746344, 0.3829702745346273, 0.39870546723886807],)

julia> b(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ self      Identifier                         File                              
      │ Float64   String                             String                            
──────┼────────────────────────────────────────────────────────────────────────────────
    1 │ 0.209127  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.328976  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.223653  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.539114  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮                      ⋮                                  ⋮
 3799 │ 0.403484  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.359463  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.372486  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                      3794 rows omitted

julia> all_genes = find_seqs(example_data_path, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> b(example_data_path, ref_seqs = (ribosomal = ribosomal_genes,), dataframe = true) # Calculate B using ribosomal genes as a reference subset
3801×3 DataFrame
  Row │ ribosomal  Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │  0.274331  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │  0.32069   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │  0.255325  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │  0.546493  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │  0.406673  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │  0.375857  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │  0.437981  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

julia> b(example_data_path, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,), dataframe = true) # Calculate B using all genes and ribosomal genes as a reference subset
3801×4 DataFrame
  Row │ self      ribosomal  Identifier                         File                              
      │ Float64   Float64    String                             String                            
──────┼───────────────────────────────────────────────────────────────────────────────────────────
    1 │ 0.209127   0.274331  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.328976   0.32069   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.223653   0.255325  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.539114   0.546493  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮          ⋮                      ⋮                                  ⋮
 3799 │ 0.403484   0.406673  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.359463   0.375857  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.372486   0.437981  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                                 3794 rows omitted
"""
function b(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end

    return b(filepath, ref_seqs, uniqueI, stop_mask, rm_start, threshold, dataframe)
end

function b(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = b(filepaths[i], ref_seqs, uniqueI, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = b(filepaths[i], ref_seqs[i], uniqueI, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function b(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    lengths =  @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),))
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Same as in MILC
    pa = map(x-> x ./ lengths, eachrow(countAA))
    pa = transpose(reduce(hcat, pa))
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Same as in MILC
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs)
    dts = map(x->abs.((normfreq .- x)), normsetfreqs) # This is good
    dts = map(dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end
    bas = map(dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end
    bs = map(bas) do ba
         vec(sum(ba .* pa, dims = 1))
    end
    dataframe && return (df = DataFrame(bs); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df)
    return bs
end

# ENC
function enc(filepaths::String, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return enc(filepaths, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function enc(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    Threads.@threads for i in 1:len
       @inbounds results[i] = enc(filepaths[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
    end
    dataframe && return reduce(vcat, results)
    return results
end

function enc(fasta_seq::String, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Same as in MILC
    pi_vec =  enc_pi(count_matrix, countAA, seqs, dict_uniqueI)
    fa =  @. (countAA * pi_vec - 1) / (countAA - 1)
    fa[isnan.(fa)] .= 0.0
    res = vec(eFFNc(fa, dict_deg))
    dataframe && return DataFrame(ENC = res, Identifier = names, File = fasta_seq)
    return res
end

# ENC Prime
function enc_p(filepaths::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return enc_p(filepaths, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function enc_p(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = enc_p(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = enc_p(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end


function enc_p(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),))
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Same as in MILC
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Same as in MILC
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    @inbounds dts = map(x->(@. (normfreq - x) ^ 2 / x), normsetfreqs) # Calculate Ma for each codon
    dts = map(dts) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    chisums = map(dts) do dt
        chisum = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            chisum[i,:] = row
        end
        chisum
    end
    chisqs = map(x-> x .* countAA, chisums) # This also looks good
    fas = map(chisqs) do chisq
            fa = @. (chisq + countAA - dict_deg) / ((countAA - 1) * dict_deg)
            fa[countAA .< 5] .= 0.0
            fa
        end
        res = map(x->vec(eFFNc(x, dict_deg)), fas)
    dataframe && return (df = DataFrame(res); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df)
    return res
end

# MCB 
function mcb(filepaths::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return mcb(filepaths, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function mcb(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = mcb(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = mcb(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function mcb(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),))
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Same as in MILC
    A = countAA[dict_deg .> 1,:] .> 0
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Same as in MILC
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    @inbounds dts = map(x->(@. (normfreq - x) ^ 2 / x), normsetfreqs)
    dts = map(dts) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    no_counts = count_matrix .<= 0 
    dts = map(dts) do dt
        dt[no_counts] .= 0
        dt
    end

    bas = map(dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end

    mcbs = map(bas) do ba
        mat1 = @. ba[dict_deg > 1,:] * log10(countAA[dict_deg > 1, :])
        mat1 = map((x) -> isnan(x) ? 0.0 : x, mat1)
        vec(sum(mat1, dims = 1) ./ sum(A, dims = 1))
    end
    dataframe && return (df = DataFrame(mcbs); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df)
    return mcbs
end

# MILC
function milc(filepaths::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return milc(filepaths, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function milc(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
    Threads.@threads for i in 1:len
       @inbounds results[i] = milc(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
    end
    else
        Threads.@threads for i in 1:len
            @inbounds results[i] = milc(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
         end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function milc(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    lengths =  @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),))
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Count the codons per amino acid
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find the frequency of each codon for each gene
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    @inbounds cor = correction_term(countAA, lengths, dict_deg) # Calculate correction term
    @inbounds per_codon_mas = map(x->(@. log(normfreq / x) * count_matrix), normsetfreqs) # Calculate Ma for each codon
    fixed_per_codon_mas = map(per_codon_mas) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    mas = map(fixed_per_codon_mas) do pcmas
        ma = Array{Float64}(undef, size(countAA, 1), size(countAA, 2)) # Pre-allocate matrix for Ma across each amino acid
        for (i,aa) in enumerate(dict_uniqueI)
            row = 2 .* sum(pcmas[aa,:],dims = 1) # Calculate ma for each amino acid
            ma[i,:] = row
        end
        ma
    end
    milcs = map(mas) do ma
        @views vec(([sum(ma, dims = 1)...] ./ lengths) .- cor) # Calculate MILC for each gene
    end
    dataframe && return (df = DataFrame(milcs); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df)
    return milcs
end

# SCUO
function scuo(filepaths::String, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return scuo(filepaths, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function scuo(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Vector{Float64}}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    Threads.@threads for i in 1:len
       @inbounds results[i] = scuo(filepaths[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
    end
    dataframe && return reduce(vcat, results)
    return results
end


function scuo(fasta_seq::String, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Same as in MILC
    Ha = scuo_freq(count_matrix, countAA, seqs, dict_uniqueI)
    
    Hmax = log10.(dict_deg)

    Oa = map(x -> (Hmax .- x) ./ Hmax, eachcol(Ha))
    Oa = reduce(hcat, Oa)

    Fa = countAA ./ sum(countAA[dict_deg .> 1, :], dims = 1)

    mult = Oa .* Fa
    mult = map((x) -> isnan(x) ? 0.0 : x, mult)

    res =  vec(sum(mult, dims = 1))
    dataframe && return DataFrame(SCUO = res, Identifier = names, File = fasta_seq)
    return res
end



function all_cub(filepaths::String, dict::codon_dict = default_codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return all_cub(filepaths, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function all_cub(filepaths::Vector{String}, dict::codon_dict = default_codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = all_cub(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = all_cub(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end


function all_cub(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    seqnames = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    lengths =  @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),))
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # This is the same for all measures
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) 
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs)
    # Up to this point, all of the measures should have been the same

    # Calculating B
    b_pa = map(x-> x ./ lengths, eachrow(countAA))
    b_pa = transpose(reduce(hcat, b_pa))
    
    b_dts = map(x->abs.((normfreq .- x)), normsetfreqs) # This is good
    b_dts = map(b_dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end
    b_bas = map(b_dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end
    B_result = map(b_bas) do ba
         vec(sum(ba .* b_pa, dims = 1))
    end
    # End calculating B (B is the result)

    # Calculate ENC
    pi_vec =  enc_pi(count_matrix, countAA, seqs, dict_uniqueI)
    enc_fa =  @. (countAA * pi_vec - 1) / (countAA - 1)
    enc_fa[isnan.(enc_fa)] .= 0.0
    ENC_result = (ENC = vec(eFFNc(enc_fa, dict_deg)),)
    # End calculating ENC (ENC is the result)

    # Calculate ENC'
    encp_dts = map(x->(@. (normfreq - x) ^ 2 / x), normsetfreqs) # Calculate Ma for each codon
    encp_dts = map(encp_dts) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    encp_chisums = map(encp_dts) do dt
        chisum = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            chisum[i,:] = row
        end
        chisum
    end
    encp_chisqs = map(x-> x .* countAA, encp_chisums) # This also looks good
    encp_fas = map(encp_chisqs) do chisq
            fa = @. (chisq + countAA - dict_deg) / ((countAA - 1) * dict_deg)
            fa[countAA .< 5] .= 0.0
            fa
        end
    ENCP_result = map(x->vec(eFFNc(x, dict_deg)), encp_fas)
    #End calculating ENC' (ENCP is the result)

    #Start calculating MCB (can use ENCP dts)
    no_counts = count_matrix .<= 0 
    mcb_dts = map(encp_dts) do dt
        dt[no_counts] .= 0
        dt
    end

    mcb_bas = map(mcb_dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end
    A = countAA[dict_deg .> 1,:] .> 0
    MCB_result = map(mcb_bas) do ba
        mat1 = @. ba[dict_deg > 1,:] * log10(countAA[dict_deg > 1, :])
        mat1 = map((x) -> isnan(x) ? 0.0 : x, mat1)
        vec(sum(mat1, dims = 1) ./ sum(A, dims = 1))
    end
        
    # End calculating MCB (MCB is the result)

    # Start calculating MILC
    cor = correction_term(countAA, lengths, dict_deg) # Calculate correction term
    per_codon_mas = map(x->(@. log(normfreq / x) * count_matrix), normsetfreqs) # Calculate Ma for each codon
    fixed_per_codon_mas = map(per_codon_mas) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    mas = map(fixed_per_codon_mas) do pcmas
        ma = Array{Float64}(undef, size(countAA, 1), size(countAA, 2)) # Pre-allocate matrix for Ma across each amino acid
        for (i,aa) in enumerate(dict_uniqueI)
            row = 2 .* sum(pcmas[aa,:],dims = 1) # Calculate ma for each amino acid
            ma[i,:] = row
        end
        ma
    end
    MILC_result = map(mas) do ma
        @views vec(([sum(ma, dims = 1)...] ./ lengths) .- cor) # Calculate MILC for each gene
    end
    # End calculating MILC (MILC is the result)

    # Start calculating SCUO
    Ha = scuo_freq(count_matrix, countAA, seqs, dict_uniqueI)
    
    Hmax = log10.(dict_deg)

    Oa = map(x -> (Hmax .- x) ./ Hmax, eachcol(Ha))
    Oa = reduce(hcat, Oa)

    Fa = countAA ./ sum(countAA[dict_deg .> 1, :], dims = 1)

    mult = Oa .* Fa
    mult = map((x) -> isnan(x) ? 0.0 : x, mult)

    SCUO_result = (SCUO = vec(sum(mult, dims = 1)),)

    if dataframe
        milc_df = DataFrame(MILC_result); rename!(milc_df, names(milc_df) .* "_MILC")
        mcb_df = DataFrame(MCB_result); rename!(mcb_df, names(mcb_df) .* "_MCB")
        encp_df = DataFrame(ENCP_result); rename!(encp_df, names(encp_df) .* "_ENCP")
        b_df = DataFrame(B_result); rename!(b_df, names(b_df) .* "_B")
        enc_df = DataFrame(ENC_result)
        scuo_df = DataFrame(SCUO_result) ; scuo_df[!,:Identifier] = seqnames; scuo_df[!,:File] .= fasta_seq
        return hcat(b_df, enc_df, encp_df, mcb_df, milc_df, scuo_df)
    end
    return (B = B_result, ENC = ENC_result, ENCP = ENCP_result, MCB = MCB_result, MILC = MILC_result, SCUO = SCUO_result)

end


