using Test, NCTSSOS, DynamicPolynomials

using NCTSSOS: _cyclic_canon, cyclic_canon, _sym_canon, sym_canon

using NCTSSOS: star

using NCTSSOS: newton_cyclic, get_basis, cc
using NCTSSOS: is_sym
using NCTSSOS: get_ncbasis
using NCTSSOS: reduce_cons!, reduce_cons, reduce!, reduce
using NCTSSOS: _comm
using NCTSSOS: bfind
using NCTSSOS: permutation, _permutation
using NCTSSOS: poly_info, polys_info
using NCTSSOS: isless_td
using NCTSSOS: sym_cyclic, _cyclic_basis, sym, iscomm, gind, res_com!, issym
using NCTSSOS: add_SOHS!, add_poly!, arrange

@testset "utils.jl Tests" begin

    @ncpolyvar x[1:4]
    # f = x[1] * x[2]^2 * x[1] + x[2] * x[1]^2 * x[2] + x[3] * x[2] * x[1]
    order = 2
	n, supp, coe = 4, [UInt16[1,2,2,1],UInt16[2,1,1,2],UInt16[3,2,1]], [1,1,1]

    @testset "get_basis" begin

        basis = get_basis(n, order)

        @test size(basis) == (n, binomial(n + order, order))
        unique_basis = unique(eachcol(basis))
        @test length(unique_basis) == size(basis, 2)
        @test all(sum(basis, dims=1) .<= order)
    end

    @testset "cyclic_canon" begin

        @test _cyclic_canon(UInt16[3, 4, 1, 2]) == [1, 2, 3, 4]
        @test _cyclic_canon(UInt16[]) == UInt16[]

        #  This behaves different from last
        # due to [1,2,3,4] < [4,1,2,3]
        # but x[1]*x[2]*x[3]*x[4] > x[4]*x[1]*x[2]*x[3]
        @test _cyclic_canon(x[3] * x[4] * x[1] * x[2]) == x[4] * x[1] * x[2] * x[3]

        empty_monomial = x[3]^0 * x[4]^0
        @test _cyclic_canon(empty_monomial) == empty_monomial


        nsupp, ncoe = cyclic_canon(supp, coe)

        @test nsupp == [UInt16[1, 1, 2, 2], UInt16[1, 2, 3]]
        @test ncoe == [2.0, 1.0]

    end

    @testset "sym_canon" begin

        # @ncpolyvar x[1:3]
        # f = x[1] * x[2] + x[2] * x[1] + x[3] * x[2]
        supp = [UInt16[1, 2], UInt16[2, 1], UInt16[3, 2]]
        coe = [1, 1, 1]

        @test _sym_canon(UInt16[3, 1, 2]) == UInt16[2, 1, 3]
        @test _sym_canon(x[3] * x[1] * x[2]) == x[3] * x[1] * x[2]

        #TODO the use of star causes two methods _sym_canon to have different coefs
        # in case star(w) is used and coefficient is imaginary
        # @test _sym_canon(UInt16[3, 1, 2]) == UInt16[2, 1, 3]
        # @test _sym_canon(1im * x[2] * x[1] * x[3]) == -im * x[3] * x[1] * x[2]

        @test sym_canon(supp, coe, type=ComplexF64) == ([UInt16[1, 2], UInt16[2, 3]], [2, 1])
    end

	@testset "is_sym" begin
		@test is_sym(UInt16[1,2,1]) 
		@test !is_sym(UInt16[1,2,3])
	end

    @testset "get_ncbasis" begin
        # Test basic case
        basis = get_ncbasis(2, 2)
        @test length(basis) == 7  # 1 (degree 0) + 2 (degree 1) + 4 (degree 2)
		@test sort(basis) == sort([UInt16[], UInt16[1], UInt16[2], UInt16[1,1], UInt16[1,2], UInt16[2,1], UInt16[2,2]])

        # Test with custom indices
        basis = get_ncbasis(2, 1, ind=[10, 20])
        @test length(basis) == 3
		@test sort(basis) == sort([UInt16[], UInt16[10], UInt16[20]])

        # Test binary case
        basis = get_ncbasis(2, 3, binary=true)
        @test length(basis) == 7  
        @test sort(basis) == sort([UInt16[], UInt16[1], UInt16[2], UInt16[2, 1], UInt16[1, 2], UInt16[1, 2, 1], UInt16[2, 1, 2]])

        # Test edge cases
        @test get_ncbasis(0, 0) == [UInt16[]]
        @test get_ncbasis(2, 0) == [UInt16[]]

    end


    @testset "reduce_cons!" begin
		# do you mean unipotent or involutory?
        word = UInt16[1,2,1,2,2,1]
        reduce_cons!(word)
        @test word == UInt16[1, 2, 1, 1]
        @test_broken word == UInt16[1,2] #is this intended?
		word = UInt16[1,2,1,2,2,1]
		reduce_cons!(word,constraint=nothing)
		@test word == UInt16[1,2,1,2,1]
    end

	@testset "reduce_cons" begin
		@ncpolyvar a[1:2]
		@test reduce_cons(a[1]*a[2]^2*a[1],constraint=nothing) == a[1]*a[2]*a[1]
		@test reduce_cons(a[1]*a[2]^2*a[1]) == a[1]*a[1]
		@test_broken reduce_cons(a[1]*a[2]^2*a[1]) == 1 # is this intended?
	end

    @testset "_comm" begin
        # Test _comm with Vector{UInt16}
        @test _comm(UInt16[1, 2, 3, 4], 2) == UInt16[1, 2, 3, 4]
        @test _comm(UInt16[3, 1, 4, 2], 2) == UInt16[1, 2, 3, 4]
        @test _comm(UInt16[5, 3, 1], 1) == UInt16[1, 5, 3]
        @test _comm(UInt16[2, 2, 1], 1) == UInt16[1, 2, 2]

        # Test _comm with Monomial
		@ncpolyvar a[1:4]
        w1 = a[1] * a[3] * a[2]
        w2 = a[4] * a[1] * a[2]
        w3 = a[3]^2 * a[1]

        @test _comm(w1, a, 2) == a[1] * a[2] * a[3]
        @test _comm(w2, a, 2) == a[1] * a[2] * a[4]
        @test _comm(w3, a, 1) == a[1] * a[3]^2
    end

	#TODO, need to add more test
    @testset "reduce!" begin
        word = UInt16[2,3,1]
        @test reduce!(word, obj="trace") == [1,2,3]

    end

    @testset "bfind" begin
        len_vars = 10
        @ncpolyvar a[1:len_vars]
        @ncpolyvar y

        @test bfind(a, len_vars, a[7], rev=true) == 7
        @test isnothing(bfind(a, len_vars, y, rev=true))
    end

    @testset "permutation" begin
        a = [1,1,2]
        perms = permutation(a)
        @test length(perms) == 3
        @test all(x -> isa(x, Vector{UInt16}), perms)
    end

    @testset "poly_info" begin
        @ncpolyvar a[1:2]
        f = a[1]*a[2] + 2*a[2]*a[1]
        n, supp, coe = poly_info(f, a)
        @test n == 2
        @test length(supp) == 2
        @test coe ≈ [1.0, 2.0]
    end

    @testset "isless_td" begin
        @test isless_td([1], [1,2]) == true
        @test isless_td([1,2], [1]) == false
        @test isless_td([1,2], [1,3]) == true
    end

    @testset "star" begin
        @ncpolyvar a[1:3]
        @test star(a[1] * a[2]) == a[2] * a[1]
        @test star(1im * a[1] * a[2] + 2 * a[2] * a[3]) == -1im * a[2] * a[1] + 2 * a[3] * a[2]
    end

    @testset "add_SOHS!" begin
        # This will need a mock model setup
        # Test basic functionality and constraints
    end

    @testset "add_poly!" begin
        # This will need a mock model setup
        # Test polynomial generation with different parameters
    end

    @testset "arrange" begin
        x = @ncpolyvar x[1:2]
        p = x[1]*x[2] + 2*x[2]*x[1]
        nmons, ncoe = arrange(p, x)
        @test length(nmons) == 1
        @test ncoe[1] ≈ 3.0
    end

end