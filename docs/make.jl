using CUBScout
using Documenter

DocMeta.setdocmeta!(CUBScout, :DocTestSetup, :(using CUBScout); recursive=true)

makedocs(;
    modules=[CUBScout],
    authors="Augustus Pendleton",
    repo="https://github.com/gus-pendleton/CUBScout.jl/blob/{commit}{path}#{line}",
    sitename="CUBScout.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gus-pendleton.github.io/CUBScout.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gus-pendleton/CUBScout.jl",
    devbranch="main",
)
