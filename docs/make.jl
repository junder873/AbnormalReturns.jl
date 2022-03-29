using AbnormalReturns
using Documenter

Documenter.makedocs(
    modules = [AbnormalReturns],
    sitename = "AbnormalReturns.jl",
    pages = [
        "Introduction" => "index.md",
        "Example" => "example.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/junder873/AbnormalReturns.jl.git",
    target = "build",
)