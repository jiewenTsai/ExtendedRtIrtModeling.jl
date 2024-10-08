using ExtendedRtIrtModeling
using Documenter

DocMeta.setdocmeta!(ExtendedRtIrtModeling, :DocTestSetup, :(using ExtendedRtIrtModeling); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
  file for
  file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [ExtendedRtIrtModeling],
    authors = "Jie-Wen Tsai <tsai.jiewen@gmail.com>",
    repo = "https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/blob/{commit}{path}#{line}",
    sitename = "ExtendedRtIrtModeling.jl",
    format = Documenter.HTML(; canonical = "https://jiewenTsai.github.io/ExtendedRtIrtModeling.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/jiewenTsai/ExtendedRtIrtModeling.jl")
