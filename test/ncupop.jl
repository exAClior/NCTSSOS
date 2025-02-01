using Test, NCTSSOS, DynamicPolynomials
using NCTSSOS: poly_info

# @testset "Unconstrained Non-commuting Term Sparse Sum of Square" begin
    n = 3
    @ncpolyvar x[1:n]
    f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

@testset "poly_info" begin
    n, supp, coe = poly_info(f, x)
    @test n == 3

    @test map(x -> Int.(x), supp) == [[1, 2, 2, 1], [3, 2, 2, 3], [1, 2, 1], [2, 2, 3], [3, 2, 2], [3, 2, 3], [1, 1], [1, 2], [2, 2], [2, 1], [2, 3], [3, 3], [3, 2]]

    @test coe == [2, 142, -2, 9, 9, -54, 1, -1, 3, -1, -1, 6, -1]
end


    opt, data = nctssos_first(f, x, newton=true, reducebasis=true, TS="MD", obj="eigen", QUIET=true)
	@test opt ≈  -0.0035512

	opt,data = nctssos_first(f, x, newton=true, TS="MD", obj="trace", QUIET=true)
    @test opt ≈ -0.0035512


# end