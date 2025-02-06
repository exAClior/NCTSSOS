using Test, NCTSSOS, DynamicPolynomials
using NCTSSOS: remove

@testset "ncupop.jl Tests" begin

    @testset "Unconstrained Non-commuting" begin
        n = 3
        @ncpolyvar x[1:n]
        f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] + 2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
            6x[3]^2 + 9x[2]^2 * x[3] + 9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]


        opt, data = nctssos_first(f, x, newton=true, reducebasis=true, TS="MD", obj="eigen", QUIET=true)
        @test opt ≈ -0.0035512 atol=1e-7

        opt, data = nctssos_first(f, x, newton=true, TS="MD", obj="trace", QUIET=true)
        @test opt ≈ -0.0035512 atol=1e-7
    end


    @testset "nctssos_higher!" begin
        # Setup initial data
        x = @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        _, data = nctssos_first(f, x)
        
        # Test higher step
        opt, data = nctssos_higher!(data)
        @test opt isa Union{Float64, Nothing}
        # Add more specific tests
    end

    @testset "cc function" begin
        a = UInt16[1, 2, 2, 3, 3, 3]
        n = 4
        result = cc(a, n)
        @test result == UInt16[1, 2, 3, 0]
    end

    @testset "remove" begin
        #TODO what exactly is the purpose?
        n = 3
        order = 2
        pbasis = get_basis(n,order)

        supp = Vector{UInt16}[[1,2],[3,1],[1]]

        csupp = cc.(supp, n)

        csupp
        csupp[1] .- 2 * pbasis[:, 1]

        pbasis

        @test !remove(csupp, pbasis[:,1],n)

    end

    @testset "newton_cyclic" begin
        supp = Vector{UInt16}[[1, 2], [2, 1], [1], [2, 2]]
        n = 2
        d = 2
        basis = newton_cyclic(supp, n, d)
        basis
        @test basis == Vector{UInt16}[[],[2]]
        # Add more specific tests
    end

    @testset "newton_ncbasis" begin
        supp = [[1, 2, 2, 1], [1, 1]]
        basis = newton_ncbasis(supp)
        @test basis == Vector{UInt16}[[], [1], [2, 1]]
    end

    @testset "get_graph" begin
        ksupp = [[1, 2], [2, 1]]
        basis = [[1], [2]]
        G = get_graph(ksupp, basis)
        @test nv(G) == length(basis)
        # Add more specific tests
    end

    @testset "get_blocks" begin
        ksupp = [[1, 2], [2, 1]]
        basis = [[1], [2]]
        blocks, cl, blocksize = get_blocks(ksupp, basis)
        @test blocks isa Vector{Vector{Int}}
        @test cl isa Int
        @test blocksize isa Vector{Int}
    end

    @testset "solvesdp" begin
        # Setup test data
        supp = [[1, 2], [2, 1]]
        coe = [1.0, 1.0]
        basis = [[1], [2]]
        blocks = [[1], [2]]
        cl = 2
        blocksize = [1, 1]
        
        # Test basic functionality
        objv, ksupp, moment, GramMat = solvesdp(supp, coe, basis, blocks, cl, blocksize)
        @test objv isa Union{Float64, Nothing}
        @test ksupp isa Vector{Vector{UInt16}}
        # Add more specific tests
    end


end