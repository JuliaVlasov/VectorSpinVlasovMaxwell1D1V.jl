module VectorSpinVlasovMaxwell1D1V

using DocStringExtensions
using FFTW

include("initialfunction.jl")
include("numeint.jl")
include("diagnostics.jl")
include("recover.jl")
include("translation.jl")
include("h2fh.jl")
include("he.jl")
include("haa.jl")
include("h3fh.jl")
include("h1f.jl")

end
