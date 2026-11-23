import Pkg
Pkg.activate(;temp=true)
Pkg.develop(path=joinpath(@__DIR__,"NTL.jl"))
for (name,version) in (("Nemo","0.57.2"),("Polynomials","4.1.6"),("BitIntegers","0.3.7"))
    Pkg.add(Pkg.PackageSpec(;name,version))
end
include("optional.jl")
