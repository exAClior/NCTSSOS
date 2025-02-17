using Test, NCTSSOS, DynamicPolynomials

using NCTSSOS: 
    _cyclic_canon, cyclic_canon, _sym_canon, sym_canon,
    star,
    newton_cyclic, get_basis, cc,
    is_sym,
    get_ncbasis,
    reduce_cons!, reduce_cons, reduce!, reduce,
    _comm,
    bfind,
    permutation, _permutation,
    poly_info, polys_info,
    isless_td,
    sym_cyclic, _cyclic_basis, sym, iscomm, gind, res_com!, issym,
    add_SOHS!, add_poly!, arrange, res_comm!

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
        a = [1, 1, 1] # this support is x1, x2,x3
        perms = permutation(a)
        @test length(perms) == factorial(length(a))
		@test sort(perms) == sort(Vector{UInt16}[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]])
    end

    @testset "poly_info" begin
        @ncpolyvar a[1:2]
        f = a[1]*a[2] + 2*a[2]*a[1]
        n, supp, coe = poly_info(f, a)
        @test n == 2
        @test length(supp) == 2
        @test coe ≈ [1.0, 2.0]
    end

	@testset "polys_info" begin
		@ncpolyvar a[1:2]
		f = a[1]*a[2] + 2*a[2]*a[1]
		g1 = a[1]^2
		g2 = a[2]^2
		n, supp, coe = polys_info([f,g1,g2],a)
		@test n == 2
		@test supp == Vector{Vector{UInt16}}[[[1,2],[2,1]],[[1,1]],[[2,2]]]
		@test coe == [[1.0,2.0],[1.0],[1.0]]
	end

    @testset "isless_td" begin
        @test isless_td([1], [1,2]) == true
        @test isless_td([1,2], [1]) == false
        @test isless_td([1,2], [1,3]) == true
    end

	@testset "sym_cyclic" begin
		@test sym_cyclic(UInt16[3,2,1]) == UInt16[1,2,3]
		@test sym_cyclic(UInt16[3,4,5,1]) == UInt16[1,3,4,5]
	end

	@testset "_cyclic_basis" begin
        #x1 * x2^2 * x3
        @test sort(_cyclic_basis(UInt16[1, 2, 3])) == sort(Vector{UInt16}[[1, 2, 3], [1, 3, 2]])
	end

	@testset "sym" begin
		# TODO
	end

	#TODO what is the purpose? RIP?
    @testset "iscomm" begin
        # Test case 1: Commuting variables within same group
        vargroup1 = [2, 2]
        a1 = [1, 2, 3, 4]
        @test iscomm(a1, vargroup1) == true

        # Test case 2: Non-commuting variables across different groups
        vargroup2 = [2, 2]
        a2 = [3, 1, 4, 2]
        @test iscomm(a2, vargroup2) == false

        # Test case 3: Empty array
        vargroup3 = [2, 2]
        a3 = []
        @test iscomm(a3, vargroup3) == true

        # Test case 4: Single element array
        vargroup4 = [2, 2]
        a4 = [1]
        @test iscomm(a4, vargroup4) == true

        # Test case 5: Variables in correct order across groups
        vargroup5 = [3, 2]
        a5 = [1, 2, 3, 4, 5]
        @test iscomm(a5, vargroup5) == true

        # Test case 6: Variables in wrong order within same group
        vargroup6 = [3, 2]
        a6 = [1, 3, 2, 4, 5]
        @test iscomm(a6, vargroup6) == true
    end

	@testset "res_com!" begin
		a = [1,4,3,5]
		vargroup = [3,2]

		@test !iscomm(a,vargroup)
		res_comm!(a,vargroup)
		@test iscomm(a,vargroup)
	end

	@testset "gind" begin 
		@test isnothing(gind(10,[1,2,3]))
		@test gind(4,[3,4,100]) == 2
		@test gind(7,[3,4,100]) == 2
		@test gind(100,[3,4,100]) == 3
	end
		

    @testset "issym" begin
        # Test case 1: Symmetric word with single group
        word1 = [1, 2, 2, 1]
        vargroup1 = [4]
        @test issym(word1, vargroup1) == true

        # Test case 2: Non-symmetric word with single group
        word2 = [1, 2, 3, 1]
        vargroup2 = [4]
        @test issym(word2, vargroup2) == false

        # Test case 3: Symmetric word with multiple groups
        word3 = [1, 2, 2, 1, 3, 4, 4, 3]
        vargroup3 = [2, 2, 4]
        @test issym(word3, vargroup3) == true

        # Test case 4: Non-symmetric word with multiple groups
        word4 = [1, 2, 3, 1, 3, 4, 4, 3]
        vargroup4 = [2, 2, 4]
        @test issym(word4, vargroup4) == false

        # Test case 5: Empty word
        word5 = []
        vargroup5 = [2, 2]
        @test issym(word5, vargroup5) == true

        # Test case 6: Single element word
        word6 = [1]
        vargroup6 = [1]
        @test issym(word6, vargroup6) == true
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
        @ncpolyvar x[1:2]
        p = x[1]*x[2] + 2*x[2]*x[1]
        nmons, ncoe = arrange(p, x)
        @test length(nmons) == 1
        @test ncoe[1] ≈ 3.0
    end

end