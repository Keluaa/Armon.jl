using Armon
using Documenter

CI = get(ENV, "CI", "false") == "true"

DocMeta.setdocmeta!(Armon, :DocTestSetup, :(using Armon); recursive=true)

makedocs(;
    modules=[Armon],
    authors="Luc Briand <luc.briand35@gmail.com> and contributors",
    repo=Remotes.GitHub("Keluaa", "Armon.jl"),
    sitename="Armon.jl",
    format=Documenter.HTML(;
        prettyurls=CI,
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "About" => "index.md",
        "Solver" => "solver.md",
        "Blocking" => "blocking.md",
        "Kernels" => "kernels.md",
        "Utils" => "utils.md"
    ],
)

if CI
    deploydocs(
        repo = "github.com/Keluaa/Armon.jl.git",
        push_preview = true
    )
end
