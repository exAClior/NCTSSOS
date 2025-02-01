using Test, NCTSSOS, DynamicPolynomials

using NCTSSOS: poly_info
using NCTSSOS: bfind
using NCTSSOS: _cyclic_canon, cyclic_canon, _sym_canon, sym_canon

using NCTSSOS: star

#TODO Why? isn't our variable supposed to be noncommuting?
# Why can you take and combine cyclic terms?
@testset "cyclic_canon" begin

	@ncpolyvar x[1:4]
    f = x[1]*x[2]^2*x[1] + x[2]*x[1]^2*x[2] + x[3]*x[2]*x[1]
	n, supp, coe = poly_info(f, x)

	@test _cyclic_canon(UInt16[3,4,1,2])== [1,2,3,4]
	@test _cyclic_canon(UInt16[]) == UInt16[]

	#  This behaves different from last
	# due to [1,2,3,4] < [4,1,2,3]
	# but x[1]*x[2]*x[3]*x[4] > x[4]*x[1]*x[2]*x[3]
    @test _cyclic_canon(x[3] * x[4] * x[1] * x[2]) == x[4] * x[1] * x[2] * x[3]

	empty_monomial = x[3]^0 * x[4]^0
    @test _cyclic_canon(empty_monomial) == empty_monomial


    nsupp, ncoe = cyclic_canon(supp, coe)

	@test nsupp == [UInt16[1,1,2,2],UInt16[1,2,3]]
	@test ncoe == [2.0,1.0]

end

@testset "Sym Canon" begin
    @ncpolyvar x[1:3]
    f = x[1] * x[2] + x[2] * x[1] + x[3] * x[2]
	n, supp, coe = poly_info(f, x)

    @test _sym_canon(UInt16[3, 1, 2]) == UInt16[2, 1, 3]
    @test _sym_canon(x[3] * x[1] * x[2]) == x[3] * x[1] * x[2]

	#TODO the use of star causes two methods _sym_canon to have different coefs
	# in case star(w) is used and coefficient is imaginary
    # @test _sym_canon(UInt16[3, 1, 2]) == UInt16[2, 1, 3]
    # @test _sym_canon(1im * x[2] * x[1] * x[3]) == -im * x[3] * x[1] * x[2]

    @test sym_canon(supp, coe, type=ComplexF64) == ([UInt16[1, 2], UInt16[2, 3]], [2, 1])

end

@testset "Hermitian Conjugate" begin
	@ncpolyvar x[1:3]
	@test star(x[1]*x[2]) == x[2] * x[1]
    @test star(1im * x[1] * x[2] + 2 * x[2] * x[3]) == -1im * x[2] * x[1] + 2 * x[3] * x[2]
end


@testset "bfind" begin
	len_vars = 10
	@ncpolyvar x[1:len_vars]
	@ncpolyvar y

	@test bfind(x, len_vars, x[7],rev=true) == 7
	@test isnothing(bfind(x,len_vars,y,rev=true))
end
