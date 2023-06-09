
# Calculating Codon Usage Bias

## Under default parameters

Codon usage bias can be calculated using the functions `b()`, `enc()`, `enc_p()`, `mcb()`, `milc()` and `scuo()`. These functions accept fasta-formatted files or vectors of BioSequences. If analyzing fasta files, these functions can accept filepaths as a string, an open IO, or an open FASTAReader. If providing a vector of BioSequences, the sequence can use either DNA or RNA alphabets. The output is a vector giving the codon usage bias of each coding sequence in the fasta file. 

!!! warning

    `CUBScout` does not identify ORFs, pause at stop codons, or parse non-nucleotide characters. It is assumed the coding sequences you provide are in-frame and don't contain 5' or 3' untranslated regions. Codons which have non-specific nucleotides, like "W", are skipped. Sequences with characters outside of those recognized by BioSequences will throw an error.

`CUBScout` is loaded with an example dataset, which can accessed at `CUBScout.EXAMPLE_DATA_PATH`. This string points to a .fna of coding sequences from *B. subtilis*. Let's calculate ENC for the genes in this file. 

```julia-repl
julia> EXAMPLE_DATA_PATH
"/your/path/to/file/B_subtilis.fna"

julia> enc_result = enc(EXAMPLE_DATA_PATH);

julia> enc_result.ENC
3801-element Vector{Float64}:
 56.787282202547104
 52.725946690067296
 59.287948966886226
 52.29668642771212
 55.26298060679466
 53.44161579771853
  ⋮
 50.30390962534221
 56.29539618087172
 55.229391962859935
 52.58401385627267
 60.19275631834157

julia> enc_result.Identifier
3801-element Vector{String}:
 "lcl|NC_000964.3_cds_NP_387882.1_1"
 "lcl|NC_000964.3_cds_NP_387883.1_2"
 "lcl|NC_000964.3_cds_NP_387885.1_4"
 "lcl|NC_000964.3_cds_NP_387886.2_5"
 "lcl|NC_000964.3_cds_NP_387887.1_6"
 "lcl|NC_000964.3_cds_NP_387888.1_7"
 ⋮
 "lcl|NC_000964.3_cds_NP_391981.1_4232"
 "lcl|NC_000964.3_cds_NP_391982.1_4233"
 "lcl|NC_000964.3_cds_NP_391983.1_4234"
 "lcl|NC_000964.3_cds_NP_391984.1_4235"
 "lcl|NC_000964.3_cds_NP_391985.1_4236"
```

ENC and SCUO calculate codon usage bias against a theoretical, unbiased distribution, and so simply return a named tuple containing ENC/SCUO and then the gene identifiers. B, ENC', MCB, and MILC calculate an expected codon frequency using a reference set of the genome, and then calculate codon usage bias for each gene against that reference set. As such, these functions return a named tuple which describes which reference set was used, alongside gene identifiers. By default, the codon usage bias is calculated against the codon usage bias of the genome as a whole, which we typically refer to as "self".

```julia-repl
julia> b_result = b(EXAMPLE_DATA_PATH)
(self = [0.20912699220973896, 0.3289759448740455, 0.22365336363593893, 0.5391135258658497, 0.24919594143501034, 0.2880358413249049, 0.31200964304415874, 0.34858035204347476, 0.2455189361074733, 0.4690734561271221  …  0.3629137353834403, 0.3621330537227321, 0.4535285720373026, 0.3357858047622507, 0.28183191395624935, 0.2668809561422238, 0.22381338105820905, 0.4034837015709619, 0.3594626865160133, 0.3724863965444541],)

julia> b_result.self
3801-element Vector{Float64}:
 0.20912699220973896
 0.3289759448740455
 0.22365336363593893
 0.5391135258658497
 0.24919594143501034
 0.2880358413249049
 ⋮
 0.2668809561422238
 0.22381338105820905
 0.4034837015709619
 0.3594626865160133
 0.3724863965444541
```

Many of these measures rely on the same initial calculations. If you want to calculate all six measures at the same time, use the function `all_cub()`. This only runs these initial calculations once before calculating individual codon usage measures, and as such is more efficient than running all the functions separately. By default, all_cub returns a named tuple, each key of which corresponds to a different codon usage bias measure.

```julia-repl
julia> all_cub_result = all_cub(EXAMPLE_DATA_PATH);

julia> all_cub_result.B.self
3801-element Vector{Float64}:
 0.20912699220973896
 0.3289759448740455
 0.22365336363593893
 0.5391135258658497
 0.24919594143501034
 0.2880358413249049
 ⋮
 0.2668809561422238
 0.22381338105820905
 0.4034837015709619
 0.3594626865160133
 0.3724863965444541

julia> all_cub_result.ENC.ENC
3801-element Vector{Float64}:
 56.787282202547104
 52.725946690067296
 59.287948966886226
 52.29668642771212
 55.26298060679466
 53.44161579771853
  ⋮
 50.30390962534221
 56.29539618087172
 55.229391962859935
 52.58401385627267
 60.19275631834157
```

## Codon Dictionaries

If you are working with genomes that use the standard genetic code, than feel free to skip this section - you should not need to worry about it. By default, `CUBScout` translates sequences using the standard code, as loaded in `CUBScout.DEFAULT_CodonDict`. However, if your sequences are translated differently, you will need to provide a custom codon dictionary to `CUBScout`. 

Codon dictionaries are of a custom type `CodonDict`. You can use `?CodonDict` to see the information this struct holds, which our codon usage bias functions need to correctly translate codons and calculate codon frequency. However, I recommend you **do not** construct a `CodonDict` manually, but instead make one using the `make_CodonDict()` function. 

`make_CodonDict` reads a plain text delimited file which lists the 64 codons and their corresponding amino acid. The file should look something like this:

```
AAA    Lysine
AAC    Asparagine
AAG    Lysine 
...    ...
```     

Please follow these formatting guidelines to make sure the table is parsed correctly:

- The first column should be codons, and the second column their corresponding amino acid. 
- Do not include headers and avoid trailing whitespace. 
- Codons do not need to be alphabetized. 
- Avoid spaces and special characters (e.g., write GlutamicAcid instead of Glutamic Acid). 
- Stop codons can be coded as Stop, stop, STOP, or *. 
- If delimited using any character outside of tab, supply the delimiter as the second argument as a `Char`, not a `String` (e.g. `','` not `","`). 

`make_CodonDict` uses `readdlm` from `DelimitedFiles`; it's a good idea to check whether `readdlm` parses your file correctly before passing to `make_CodonDict`.

For demonstration purposes, `CUBScout` includes the delimited file used to construct the `DEFAULT_CodonDict`. 

```julia-repl
julia> CodonDict_PATH
"your/path/to/codon_dict.txt"

julia> our_CodonDict = make_CodonDict(CodonDict_PATH);

julia> our_CodonDict.codons
64-element Vector{BioSequences.LongSequence{BioSequences.DNAAlphabet{2}}}:
 AAA
 AAC
 AAG
 AAT
[...]

julia> our_CodonDict.AA
64-element Vector{String}:
 "Lysine"
 "Asparagine"
 "Lysine"
 "Asparagine"
 [...]
```

You can supply your custom codon dictionary to any of the codon usage bias functions as the second argument.

```julia-repl
julia> milc(EXAMPLE_DATA_PATH, our_CodonDict)
(self = [0.49482573202153163, 0.5839439121281993, 0.49947166558087047, 0.6354929447434434, 0.5439352548027006, 0.6104721251245075, 0.6256398806438782, 0.6228376952086359, 0.5355298113407091, 0.7832276821181443  …  0.5968814155010973, 0.5964500002803941, 0.5930680822246766, 0.5412999510428169, 0.49866919389111675, 0.5830959504630727, 0.5139438478694085, 0.6164434557282711, 0.6018041071661588, 0.48775477465069617],)
```

## Alternative Start Codons

`CUBScout` provides three options to handle start codons outside of ATG. By default, alternative codons are translated as they would be anywhere else in the sequence. As such, TTG would be counted as leucine, even when in the first position.

If you would like to disregard start codons entirely, set the argument `rm_start = true`. This will decrease the length of each gene sequence by one, but is my preferred method to dealing with alternative start codons.

Other packages to calculate codon usage bias, such as [coRdon](https://www.bioconductor.org/packages/release/bioc/html/coRdon.html), handle alternative start codons differently. They encode all TTG and CTG codons as methionine, regardless of their location in the gene. While I disagree with this approach from a biological perspective, you can implement it using the pre-loaded `ALTSTART_CodonDict`. 

```julia-repl
julia> scuo_result = scuo(EXAMPLE_DATA_PATH, ALTSTART_CodonDict);

julia> scuo_result.SCUO
3801-element Vector{Float64}:
 0.14286111587263958
 0.19315278493814017
 0.0966128845976179
 0.3473543659821751
 0.10792236840320082
 0.12039525638448735
 ⋮
 0.152064610300728
 0.11200912387676948
 0.18952246579743504
 0.16473723774598686
 0.24160824180945173
```
## Custom Reference Sets with `ref_seqs`

B, ENC', MCB, and MILC all calculate an expected codon frequency using a reference set of the genome, and then calculate codon usage bias for each gene against that reference set. By default, this is the entire genome ("self"). However, you can provide your own reference subset(s) to these functions. 

First, you'll need a Boolean vector, whose length matches the number of sequences in your fasta file. Genes which you want included in your subset should be `true`; the rest of the vector should be `false`. One way to make this vector is with the `find_seqs` function to look for genes with specific functions.

```julia-repl
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal")
4237-element Vector{Bool}:
 0
 0
 0
 0
 0
 0
 ⋮
 0
 0
 0
 0
 1
```

!!! tip
    `CUBScout` is designed not to hold the data from your fasta file as an object in your Julia environment. If you want to get sequence identifiers or descriptions outside of codon usage bias functions, there are the convenience functions `seq_names` and `seq_descriptions`:
    
    ```julia-repl
    julia> seq_names(EXAMPLE_DATA_PATH)[1:5]
    5-element Vector{String}:
     "lcl|NC_000964.3_cds_NP_387882.1_1"
     "lcl|NC_000964.3_cds_NP_387883.1_2"
     "lcl|NC_000964.3_cds_NP_387884.1_3"
     "lcl|NC_000964.3_cds_NP_387885.1_4"
     "lcl|NC_000964.3_cds_NP_387886.2_5"

    julia> seq_descriptions(EXAMPLE_DATA_PATH)[1]
    "lcl|NC_000964.3_cds_NP_387882.1_1 [gene=dnaA] [locus_tag=BSU_00010] [db_xref=EnsemblGenomes-Gn:BSU00010,EnsemblGenomes-Tr:CAB11777,GOA:P05648,InterPro:IPR001957,InterPro:IPR003593,InterPro:IPR010921,InterPro:IPR013159,InterPro:IPR013317,InterPro:IPR018312,InterPro:IPR020591,InterPro:IPR024633,InterPro:IPR027417,PDB:4TPS,SubtiList:BG10065,UniProtKB/Swiss-Prot:P05648] [protein=chromosomal replication initiator informational ATPase] [protein_id=NP_387882.1] [location=410..1750] [gbkey=CDS]"
    ```
Once you have your reference vector, you can supply an argument to `ref_seqs` as a named tuple. If you have multiple reference sets you want to use, those can be included as additional entries in the `ref_seqs` tuple.

```julia-repl
julia> b_ribo_result = b(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes,));

julia> b_ribo_result.ribosomal
3801-element Vector{Float64}:
 0.27433079214149625
 0.3206897249908304
 0.25532544766240484
 0.5464925047248634
 0.22424329272203575
 0.22684609299155567
 ⋮
 0.2561376033448253
 0.2217345501228918
 0.40667338789742696
 0.3758568749612823
 0.4379807676614555

julia> dna_genes = find_seqs(EXAMPLE_DATA_PATH, r"dna|DNA|Dna")
4237-element Vector{Bool}:
 1
 1
 0
 1
 0
 1
 ⋮
 0
 1
 0
 0
 0

julia> b_multi_result = b(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes, DNA = dna_genes));

julia> b_multi_result.ribosomal
3801-element Vector{Float64}:
 0.27433079214149625
 0.3206897249908304
 0.25532544766240484
 0.5464925047248634
 0.22424329272203575
 0.22684609299155567
 ⋮
 0.2561376033448253
 0.2217345501228918
 0.40667338789742696
 0.3758568749612823
 0.4379807676614555

julia> b_multi_result.DNA
3801-element Vector{Float64}:
 0.2148821062833632
 0.3182032724315858
 0.23577274334969703
 0.5371269155669846
 0.2684310325581909
 0.2860168153422687
 ⋮
 0.273137416897346
 0.21136319951043495
 0.3866134722044515
 0.3510891124098759
 0.3668966776242405
```

## Other Arguments

### `rm_stop`

Whether to remove stop codons from the calculation of codon usage bias. Default is `false`

### `threshold`

The minimum length of a gene in codons to be used when calculating codon usage bias. The default is 80; all genes under that length are discarded. If you want to discard no genes, set `threshold = 0`. You do **not** need to adjust your reference sequence vector when adjusting threshold values.

```julia-repl
julia> b_result_0 = b(EXAMPLE_DATA_PATH, threshold = 0);

julia> b_result_300 = b(EXAMPLE_DATA_PATH, threshold = 300);

julia> length(b_result_0.self)
4237

julia> length(b_result_300.self)
1650
```

### `names`

If providing a vector of BioSequences, `CUBScout` won't be able to provide identifiers for codon usage bias results. As such, you can optionally provide a vector of identifiers as an argument, so you can link results to the original input sequences. 

## Analyzing Multiple Files

Often, you might have a directory containing multiple .fna files, each of which you want to analyze. You can provide a vector of filepaths (or FASTAReaders, or IOStreams) to any `CUBScout` function, which will return a vector of results. If using BioSequences, each vector of sequences is considered a genome; if you provide a `Vector{<:Vector{<:NucSeq}}`, this will function the same as providing multiple filepaths. If supplying `ref_seqs`, provide a vector of named tuples corresponding to each file. The same goes for providing `names` - provide a `Vector{Vector{String}}` where each vector of names corresponds to each vector of Biosequences. `CUBScout` is multi-threaded, and if Julia is started with multiple threads, will assign individual threads to process individual files. This means you *should not broadcast* `CUBScout` codon usage bias functions as it will reduce efficiency. Also each file is only ever processed by a single thread, so using more threads than you have files is unnecessary. 

```julia-repl
julia> enc_p([EXAMPLE_DATA_PATH,EXAMPLE_DATA_PATH])
2-element Vector{Any}:
 self = [61.0, 59.36979815371983, 60.7494622549966, 61.0, ...],
 Identifier = ["lcl|NC_000964.3_cds_NP_387882.1_1", "lcl|NC_000964.3_cds_NP_387883.1_2", ...]),
 self = [61.0, 59.36979815371983, 60.7494622549966, 61.0, ...],
 Identifier = ["lcl|NC_000964.3_cds_NP_387882.1_1", "lcl|NC_000964.3_cds_NP_387883.1_2", ...])

julia> enc_p([EXAMPLE_DATA_PATH,EXAMPLE_DATA_PATH], ref_seqs = [(ribosomal = ribosomal_genes,), (ribosomal = ribosomal_genes,)])
2-element Vector{Any}:
 self = [61.0, 58.88817312982425, 56.41038374603565, 61.0, ...],
 Identifier = ["lcl|NC_000964.3_cds_NP_387882.1_1", "lcl|NC_000964.3_cds_NP_387883.1_2", ...]),
 self = [61.0, 58.88817312982425, 56.41038374603565, 61.0, ...],
 Identifier = ["lcl|NC_000964.3_cds_NP_387882.1_1", "lcl|NC_000964.3_cds_NP_387883.1_2", ...])
```