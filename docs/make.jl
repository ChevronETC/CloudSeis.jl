using CloudSeis, Documenter, TeaSeis

makedocs(sitename="CloudSeis", modules=[CloudSeis])
 
deploydocs(
    repo = "git@github.com:ChevronETC/CloudSeis.jl.git"
)