# B

using DataFrames: DataFrame

function b(filepath::String, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end

    return b(filepath, ref_seqs, uniqueI, stop_mask, rm_start, threshold, dataframe)
end

function b(filepaths::Vector{String}, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
function enc(filepaths::String, dict::codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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

function enc(filepaths::Vector{String}, dict::codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
function enc_p(filepaths::String, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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

function enc_p(filepaths::Vector{String}, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
function mcb(filepaths::String, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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

function mcb(filepaths::Vector{String}, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe)
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
function milc(filepaths::String, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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

function milc(filepaths::Vector{String}, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
function scuo(filepaths::String, dict::codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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

function scuo(filepaths::Vector{String}, dict::codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
scuo(path, default_dict, dataframe = true)


function all_cub(filepaths::String, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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

function all_cub(filepaths::Vector{String}, dict::codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe)
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


