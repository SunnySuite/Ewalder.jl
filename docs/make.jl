using Ewalder
using Documenter

DocMeta.setdocmeta!(Ewalder, :DocTestSetup, :(using Ewalder); recursive=true)

makedocs(;
    modules=[Ewalder],
    authors="Kipton Barros <kbarros@gmail.com> and contributors",
    repo="https://github.com/kbarros/Ewalder.jl/blob/{commit}{path}#{line}",
    sitename="Ewalder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kbarros.github.io/Ewalder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kbarros/Ewalder.jl",
    devbranch="main",
)
