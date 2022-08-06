using VectorSpinVlasovMaxwell1D1V
using Documenter

DocMeta.setdocmeta!(VectorSpinVlasovMaxwell1D1V, :DocTestSetup, :(using VectorSpinVlasovMaxwell1D1V); recursive=true)

makedocs(;
    modules=[VectorSpinVlasovMaxwell1D1V],
    authors="Julia Vlasov",
    repo="https://github.com/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl/blob/{commit}{path}#{line}",
    sitename="VectorSpinVlasovMaxwell1D1V.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaVlasov.github.io/VectorSpinVlasovMaxwell1D1V.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl",
    devbranch="main",
)
