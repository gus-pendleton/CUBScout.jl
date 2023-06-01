# CUBScout

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gus-pendleton.github.io/CUBScout.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gus-pendleton.github.io/CUBScout.jl/dev/)
[![Build Status](https://github.com/gus-pendleton/CUBScout.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gus-pendleton/CUBScout.jl/actions/workflows/CI.yml?query=branch%3Amain)

**C**odon **U**sage **B**ias (CUB) in Julia

`CUBScout` calculates **C**odon **U**sage **B**ias (CUB) and related expressivity predictions. Currently, `CUBScout` calculates:

- Six measures of codon usage bias:
    - B, from Karlin and Mrazek, 1996
    - ENC, from Wright 1990
    - ENC', from Novembre, 2002
    - MCB, from Urrutia and Hurst, 2001
    - MILC, from Supek and Vlahovicek, 2005
    - SCUO, from Wan et al., 2004
- Five expressivity measures based on codon usage bias:
    - CAI, from Sharp and Li, 1987
    - E, from Karlin and Mrazek, 1996
    - FOP, from Ikemura, 1981
    - GCB, from Merkl, 2003
    - MELP, from Supek and Vlahovicek, 2005

`CUBScout` is based off of the fabulous [coRdon](https://www.bioconductor.org/packages/release/bioc/html/coRdon.html) package in R by Anamaria Elek, Maja Kuzman, and Kristian Vlahovicek. I am grateful for their clear code and would encourage you to cite coRdon as well when using `CUBScout`.

You can install `CUBScout` by:

```julia
] add CUBScout#main
```
or
```julia
using Pkg
Pkg.add(url = "https://github.com/gus-pendleton/CUBScout.jl")
```


CUBScout is under active development, and I welcome contributions or suggestions! Additional features I'm working on/would like to incorporate:

- Performance improvements
- Plotting support (e.g. BPlots)
- Additional CUB measures, including S, RCDI, CDC, RCA, RCSU, and RCBS
- Growth predictions derived from CUB, such as those in growthpred and gRodon