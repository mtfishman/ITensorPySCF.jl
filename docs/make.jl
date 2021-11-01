using ITensorPySCF
using Documenter

DocMeta.setdocmeta!(ITensorPySCF, :DocTestSetup, :(using ITensorPySCF); recursive=true)

makedocs(;
    modules=[ITensorPySCF],
    authors="Matthew Fishman <mfishman@flatironinstitute.org> and contributors",
    repo="https://github.com/mtfishman/ITensorPySCF.jl/blob/{commit}{path}#{line}",
    sitename="ITensorPySCF.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mtfishman.github.io/ITensorPySCF.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mtfishman/ITensorPySCF.jl",
    devbranch="main",
)
