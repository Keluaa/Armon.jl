using Armon
using Documenter

include("render_plantuml.jl")

render_plantuml_files([
    "src/assets/structure.puml",
    "src/assets/global_time_step.puml",
    "src/assets/halo_exchange.puml",
    "src/assets/block_grid.puml"
])

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
        "Home" => "index.md",
        "Parameters & Device" => "parameters_and_device.md",
        "Grid" => "grid.md",
        "Blocks" => "blocks.md",
        "Block state" => "block_states.md",
        "Kernels" => "kernels.md",
        "Utils" => "utils.md",
        "UML Structure" => "structure.md"
    ],
)

if CI
    deploydocs(
        repo = "github.com/Keluaa/Armon.jl.git",
        push_preview = true
    )
end
