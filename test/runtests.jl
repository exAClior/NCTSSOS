using NCTSSOS
using Test

@testset "NCTSSOS.jl" begin
	# include("test.jl")
	include("pop.jl")
	include("solver_utils.jl")
	include("moment_solver.jl")
	include("sos_solver.jl")
end

