using Documenter
import Mosek

makedocs(
    format = :html,
    sitename = "Mosek",
    pages = [
        "Index"         => "index.md",
        "API Reference" => "Mosek-Functions.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaOpt/Mosek.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
