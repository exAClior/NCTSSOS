using Test
using NCTSSOS  # Assuming this is your module name

@testset "ncupop.jl Tests" begin

    @testset "ncupop_type struct" begin
        # Test initialization and field types
        pop = ncupop_type([], [], 0, nothing, [], "eigen", [], [], [])
        @test pop isa ncupop_type
        # Add more specific field tests
    end

    @testset "cosmo_para struct" begin
        # Test default constructor and field values
        params = cosmo_para()
        @test params.eps_abs == 1e-5
        @test params.eps_rel == 1e-5
        @test params.max_iter == 1e4
    end

    @testset "nctssos_first" begin
        # Setup test data
        x = @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        
        # Test basic functionality
        opt, data = nctssos_first(f, x)
        @test opt isa Float64
        @test data isa ncupop_type
        
        # Test different options
        opt2, data2 = nctssos_first(f, x, obj="trace", partition=1)
        @test opt2 isa Float64
        # Add more specific tests
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
        a = [1, 2, 2, 3, 3, 3]
        n = 4
        result = cc(a, n)
        @test result == [0, 1, 2, 3]
    end

    @testset "newton_cyclic" begin
        supp = [[1, 2], [2, 1]]
        n = 2
        d = 2
        basis = newton_cyclic(supp, n, d)
        @test basis isa Vector{Vector{UInt16}}
        # Add more specific tests
    end

    @testset "newton_ncbasis" begin
        supp = [[1, 2, 2, 1], [1, 1]]
        basis = newton_ncbasis(supp)
        @test basis isa Vector{Vector{UInt16}}
        # Add more specific tests
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