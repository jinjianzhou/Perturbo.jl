using Documenter, Perturbo

makedocs(;
    modules=[Perturbo],
    format=Documenter.HTML(prettyurls = false),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jinjianzhou/Perturbo.jl/blob/{commit}{path}#L{line}",
    sitename="Perturbo.jl",
    authors="Jin-Jian Zhou <jjchou.comphy@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/jinjianzhou/Perturbo.jl",
)
