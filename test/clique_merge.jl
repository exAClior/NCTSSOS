using Test
using NCTSSOS  # Assuming this is your module name
using LightGraphs
using MetaGraphs

@testset "clique_merge.jl Tests" begin

    @testset "fweight function" begin
        a = [1, 2, 3]
        b = [2, 3, 4]
        
        # Test basic weight calculation
        @test fweight(a, b) == 27 + 27 - 64  # 3^3 + 3^3 - 4^3
        
        # Test with different dimension parameter
        @test fweight(a, b, d=2) == 9 + 9 - 16  # 3^2 + 3^2 - 4^2
        
        # Test with empty sets
        @test fweight([], []) == 0
    end

    @testset "mergeclique! function" begin
        cliques = [[1, 2], [2, 3], [4, 5]]
        cliqueG = MetaGraph(3)
        add_edge!(cliqueG, 1, 2)
        add_edge!(cliqueG, 2, 3)
        
        # Test merging cliques
        new_cliques, new_graph = mergeclique!(copy(cliques), copy(cliqueG), 1, 2)
        
        @test length(new_cliques) == 2
        @test new_cliques[1] == [1, 2, 3]
        @test nv(new_graph) == 2
        @test ne(new_graph) == 1
        
        # Test weight updates
        weight = get_prop(new_graph, 1, 2, :weight)
        expected_weight = fweight(new_cliques[1], new_cliques[2])
        @test weight == expected_weight
    end

    @testset "clique_merge! function" begin
        cliques = [[1, 2], [2, 3], [3, 4], [5, 6]]
        
        # Test basic merging
        merged_cliques, cql, cliquesize = clique_merge!(copy(cliques))
        
        @test length(merged_cliques) <= length(cliques)
        @test cql == length(merged_cliques)
        @test cliquesize == length.(merged_cliques)
        
        # Test with QUIET=false
        merged_cliques, cql, cliquesize = clique_merge!(copy(cliques), QUIET=false)
        @test length(merged_cliques) <= length(cliques)
        
        # Test with different dimension parameter
        merged_cliques_d2, _, _ = clique_merge!(copy(cliques), d=2)
        @test length(merged_cliques_d2) != length(merged_cliques)
        
        # Test with non-overlapping cliques
        non_overlapping = [[1, 2], [3, 4], [5, 6]]
        merged_non_overlapping, _, _ = clique_merge!(copy(non_overlapping))
        @test merged_non_overlapping == non_overlapping
    end

end