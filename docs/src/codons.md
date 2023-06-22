# Counting Codons and Calculating Codon Frequency

In the process of calculating codon usage bias, `CUBScout` developed some handy functions for counting codons, calculating codon frequency, and sorting through a fasta file. These accessory functions underlie the codon usage bias functions, but are available to users as well. If you just want to calculate codon usage bias, you won't need to interact with these functions and can skip this section.

## Counting Codons

The first step in calculating codon usage bias is counting the occurence of each codon within each coding sequence. The `count_codons` performs this function quickly and can accept numerous datatypes.

#### Inputs and Outputs

`count_codons` can accept a filepath to a fasta file of sequences, an IOStream, an open FASTAReader, or a BioSequences nucleotide sequence, or a vector of nucleotide sequences. If providing a single sequence, the result will be a 64x1 Matrix, which corresponds to the 64 codons in alphabetical order. If you want a list of the codons in alphabetical order, this is stored in `CUBScout.DEFAULT_CodonDict.codons`. If analyzing a fasta file or a vector of sequences, the result will be a tuple. The first element of the tuple is a 64xn matrix, where n is the number of sequences above the length threshold. The second element is a vector of names corresponding to each column. The third element is a Boolean vector where `true` corresponds to sequences which did pass the threshold, and `false` is sequences which did not pass the threshold and so are not included in the results matrix. 

`CUBScout` is loaded with an example dataset, which can accessed at `CUBScout.EXAMPLE_DATA_PATH`. This string points to a .fna of coding sequences from *B. subtilis*. Let's go ahead and run `count_codons` on our example data set:

```julia-repl
julia> result = count_codons(EXAMPLE_DATA_PATH);

julia> size(result[1])
(64, 4237)

julia> result[1]
64×4237 Matrix{Int32}:
 32  21  1  17  7  33  38  …  33  12  13  14  11  6
  7   8  2   7  2  20  19     10   7   5   3   1  2
  6   3  2   9  1  10   9      7   8  10   8   8  0
 14   7  2   7  1  13  16     12   7   4   4   5  1
 11   9  2   8  2  20  23     25  11   5   2   2  1
  ⋮                 ⋮      ⋱                   ⋮  
  0   0  0   0  0   1   0      1   2   0   1   0  0
 12  11  3   6  3  10  12  …   4  11   3  10   3  2
  4   6  1   3  0   6  14      3   2   1   8   1  2
  6   6  0   7  1   5   8     12   6   7   9   2  0
 17   8  3   6  6   9   2      5   5   1   4   4  0

julia> first(result[2], 5)
5-element Vector{String}:
 "lcl|NC_000964.3_cds_NP_387882.1_1"
 "lcl|NC_000964.3_cds_NP_387883.1_2"
 "lcl|NC_000964.3_cds_NP_387884.1_3"
 "lcl|NC_000964.3_cds_NP_387885.1_4"
 "lcl|NC_000964.3_cds_NP_387886.2_5" 
```
Each column of the matrix corresponds to a gene in the genome, and each row is a codon. The value of each entry corresponds to the count of that codon. So there are 32 AAA codons in the first gene of the example genome. This gene has the identifier "lcl|NC_000964.3_cds_NP_387882.1_1".

There are a few more arguments you can use to tune how count_codons works.

#### rm_start

Whether to remove the first codon from any counts. A more thorough discussion of `CUBScout`'s treatment of start codons can be found in the Codon Usage Bias section "Alternative Start Codons". 

#### threshold

This is the minimum length of a sequence, in codons, to be included in the results. If you would like all sequences analyzed, this can be set to 0. If a sequence does not pass the threshold, it won't be included in the results matrix, and its identifier won't be included in the names vector. As we increase the threshold, our matrix becomes smaller as sequences are filtered out. Note that when supplying a filepath as input, `count_codons` dispatches on argument position, so in order to edit the threshold I also need to include the rm_start argument (`false`):

```julia-repl
julia> result_300 = count_codons(EXAMPLE_DATA_PATH, false, 300);

julia> size(result_300[1])
(64, 1650)
```

#### names

When using `count_codons` on a fasta file, identifiers are automatically pulled from the sequence headers. However, if providing a vector of BioSequences, there aren't names linked to each sequence, and instead identifiers might be stored in a separate vector of the same length. Because `count_codons` removes sequences below the threshold, this poses a quandary: how do you know which sequences were kept and which were removed? By providing a names vector to `count_codons`, as sequences are discarded, their corresponding identifiers will also be discarded from the results. This is an optional argument, and is only relevant if providing `Vector{<:NucSeq}` as an argument. For this reason, however, `count_codons` dispatches on kwargs if provided a `Vector{<:NucSeq}`, instead of position:

```julia-repl
julia> using BioSequences: @dna_str

julia> example_dna1 = dna"ATGAAAATGAACTTTTGA"
18nt DNA Sequence:
ATGAAAATGAACTTTTGA

julia> count_codons(example_dna)
64×1 Matrix{Int64}:
 1
 1
 0
 0
 0
 ⋮
 0
 0
 0
 0
 1

julia> example_dna2 = dna"ATGAAAATGAACTTTTGAATGAAAATGAACTTTTGAATGAAAATGAACTTTTGA"
54nt DNA Sequence:
ATGAAAATGAACTTTTGAATGAAAATGAACTTTTGAATGAAAATGAACTTTTGA

julia> count_codons([example_dna1, example_dna2])
(Int32[1 3; 1 3; … ; 0 0; 1 3], nothing, Bool[1, 1])

julia> count_codons([example_dna1, example_dna2], threshold = 8)
(Int32[3; 3; … ; 0; 3;;], nothing, Bool[0, 1])

julia> count_codons([example_dna1, example_dna2], threshold = 8, names = ["Example1","Example2"])
(Int32[3; 3; … ; 0; 3;;], ["Example2"], Bool[0, 1])
```

## Codon Frequency

Once you've counted your codons, you may want to calculate the *frequency* at which each codon occurs. `codon_frequency` accepts a count matrix from `count_codons` and calculates codon frequency. The input must be a 64xn matrix of integers, where each row corresponds to the 64 codons in order, and each column corresponds do a gene sequence. 

There are four ways to calculate codon frequency. Let's use an example to illustrate their differences by making a mock genome with two genes in it.

```julia-repl
julia> gene1 = dna"AAAAAGAAATTT"
12nt DNA Sequence:
AAAAAGAAATTT

julia> gene2 = dna"AAGAAGTTTTTCTTC"
15nt DNA Sequence:
AAGAAGTTTTTCTTC

julia> genome = [gene1, gene2]
2-element Vector{LongSequence{DNAAlphabet{4}}}:
 AAAAAGAAATTT
 AAGAAGTTTTTCTTC

julia> count_result = count_codons(genome)
(Int32[2 0; 0 0; … ; 0 0; 1 1], nothing, Bool[1, 1])

julia> count_matrix = count_result[1]
64×2 Matrix{Int32}:
 2  0
 0  0
 1  2
 0  0
 0  0
 ⋮  
 0  0
 0  0
 0  2
 0  0
 1  1
```

We can see that gene1 (column 1) had two counts of AAA (row 1), one of AAG (row 3), and one of TTT (row 64). For gene2 (column 2), we have 2 counts of AAG, 2 counts of TTC (row 62), and 1 count of TTT. AAA and AAG code for lysine, while TTT and TTC code for phenylalanine.

#### "net_genomic"

This will calculate the cumulative codon frequency of each codon across the entire matrix (genome), as a percentage of all codon counts in the matrix. So for example:

```julia-repl
julia> codon_frequency(count_matrix, "net_genomic")
64×1 Matrix{Float64}:
 0.2222222222222222
 0.0
 0.3333333333333333
 0.0
 0.0
 ⋮
 0.0
 0.0
 0.2222222222222222
 0.0
 0.2222222222222222
```
Across our entire genome, we counted 9 codons. AAA, TTC, and TTT all had two counts, and so at rows 1, 62, and 64 we have 0.222 = 2/9. AAG occurred three times, so we have 0.33 at row 3. 

#### "net_gene"

This will calculate the codon frequency within each gene (column) in the matrix, as a percentage of all codon counts in that gene.

```julia-repl
julia> codon_frequency(count_matrix, "net_gene")
64×2 Matrix{Float64}:
 0.5   0.0
 0.0   0.0
 0.25  0.4
 0.0   0.0
 0.0   0.0
 ⋮     
 0.0   0.0
 0.0   0.0
 0.0   0.4
 0.0   0.0
 0.25  0.2
```
First, we now have a matrix instead of a vector, because we are getting results for each gene instead of summarized across the genome. There were 4 codons in gene 1. There were two counts of AAA, so its codon frequency was 0.5 = 2/4 (row 1). AAG and TTT both had a by-gene frequency of 0.25.  There were 5 codons in gene 2. Both AAG (row 3) and TTC (row 62) had 2 counts, so had a by-gene frequency of 0.4 = 2/5.

#### "byAA_genomic"

This will calculate the codon frequency of each codon within each amino acid across the entire matrix (genome). Let's see what happens in this scenario in our mock genome:

```julia-repl
julia> codon_frequency(count_matrix, "byAA_genomic")
64-element Vector{Real}:
 0.4
 0
 0.6
 0
 0
 ⋮
 0
 0
 0.5
 0
 0.5
```
Because we've summarized across the genome we're back to a vector as output. We had 5 codons which coded for lysine across the genome. If we look at our lysine codons (AAA and AAG), we see AAA had a byAA genomic frequency of 0.4 = 2 / 5, while AAG had a byAA genomic frequency of 0.6 = 3/5. For our phenylalanine codons (TTC and TTT), they both occurred twice across the genome, and so had a byAA genomic frequency of 0.5.

#### "byAA_gene"

Finally, we can calculate the codon frequency of each codon within each amino acid within each gene. For example:

```julia-repl
julia> codon_frequency(count_matrix, "byAA_gene")
64×2 Matrix{Real}:
 0.666667  0.0
 0         0
 0.333333  1.0
 0         0
 0         0
 ⋮         
 0         0
 0         0
 0.0       0.666667
 0         0
 1.0       0.333333
```
If we look at gene1 (column 1), it exclusively used TTT to code for phenylalanine, but used AAA to code for lysine 66% of the time. In contrast, gene2 used AAG to code for lysine 100% of the time. 

#### Codon Dictionarys

If calculating byAA codon frequencies, codons need to be translated into amino acids, which is done by supplying a codon dictionary. A more complete description of codon dictionaries, including using alternative genetic codes, can be found in the Codon Usage Bias section "Codon Dictionaries".