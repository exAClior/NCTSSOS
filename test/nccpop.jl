using Test, NCTSSOS, DynamicPolynomials

@testset "Constrained Non-commuting Term Sparse Sum of Square" begin
    n = 2
    @ncpolyvar x[1:2]
    f = 2 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4 - x[1]^2 - x[2]^2
    h = x[1] * x[2] + x[2] * x[1] - 2
    pop = [f, g, h]

    opt, data = nctssos_first(pop, x, 2, numeq=1, TS="MD", obj="eigen", QUIET=true)
    # optimum = -0.9999999

    opt, data = nctssos_higher!(data, TS="MD", QUIET=true)
    # optimum = -0.9999999

    opt, data = nctssos_first(pop, x, 2, numeq=1, TS="MD", obj="trace", QUIET=true)
    # optimum = -0.9999999

    opt, data = nctssos_higher!(data, TS="MD", QUIET=true)
    # optimum = -0.9999999
end


@testset "Constrained Non-commuting Correlative Sparse Sum of Square" begin

n = 10
@ncpolyvar x[1:n]
f = 0.0
for i = 1:n
    jset = max(1,i-5):min(n,i+1)
    jset = setdiff(jset,i)
    f += (2x[i]+5*x[i]^3+1)^2
    f -= sum([4x[i]*x[j]+10x[i]^3*x[j]+2x[j]+4x[i]*x[j]^2+10x[i]^3*x[j]^2+2x[j]^2 for j in jset])
    f += sum([x[j]*x[k]+2x[j]^2*x[k]+x[j]^2*x[k]^2 for j in jset for k in jset])
end

opt,data = cs_nctssos_first([f], x, 3, TS="MD", obj="trace", QUIET=true)
# optimum = 6.5e-8

opt,data = cs_nctssos_higher!(data, TS="MD", QUIET=true)
# optimum = 6.5e-7

pop = [f]
for i = 1:n
    push!(pop, 1-x[i]^2)
    push!(pop, x[i]-1/3)
end

opt,data = cs_nctssos_first(pop, x, 3, TS="MD", obj="eigen", QUIET=true)
# optimum = 3.011288

opt,data = cs_nctssos_first(pop, x, 3, TS="MD", obj="trace", QUIET=true)
# optimum = 3.011288

end

using Test
using DynamicPolynomials
using MosekTools
using COSMO

@testset "NCTSSOS Tests" begin

    @testset "Basic Polynomial Optimization" begin
        # Define non-commuting variables
        @ncpolyvar x[1:2]
        
        # Define a simple polynomial
        p = x[1]^2 + x[2]^2 - x[1]*x[2] + x[2]*x[1]
        
        # First step of NCTSSOS
        opt, data = nctssos_first([p], x, 2, QUIET=true)
        
        @test opt !== nothing
        @test typeof(data) == nccpop_data
        @test data.n == 2
        @test data.m == 0
        @test data.order == 2
    end

    @testset "Constrained Optimization" begin
        @ncpolyvar x[1:2]
        
        # Objective and constraint
        p = x[1]^2 + x[2]^2
        g = x[1] + x[2] - 1
        
        # First step with constraint
        opt, data = nctssos_first([p, g], x, 2, QUIET=true)
        
        @test opt !== nothing
        @test data.m == 1
    end

    @testset "Higher Order Relaxation" begin
        @ncpolyvar x[1:2]
        p = x[1]^2 + x[2]^2 - x[1]*x[2] + x[2]*x[1]
        
        # First step
        opt, data = nctssos_first([p], x, 2, QUIET=true)
        
        # Higher order step
        opt_higher, data = nctssos_higher!(data, QUIET=true)
        
        @test opt_higher !== nothing
        @test opt_higher <= opt  # Higher order should be tighter or equal
    end

    @testset "Different Solvers" begin
        @ncpolyvar x[1:2]
        p = x[1]^2 + x[2]^2
        
        # Test Mosek solver
        opt_mosek, _ = nctssos_first([p], x, 2, solver="Mosek", QUIET=true)
        
        # Test COSMO solver
        opt_cosmo, _ = nctssos_first([p], x, 2, solver="COSMO", QUIET=true)
        
        @test opt_mosek !== nothing
        @test opt_cosmo !== nothing
        @test abs(opt_mosek - opt_cosmo) < 1e-4  # Results should be similar
    end

    @testset "Trace vs Eigen Objective" begin
        @ncpolyvar x[1:2]
        p = x[1]^2 + x[2]^2 - x[1]*x[2] + x[2]*x[1]
        
        # Eigen objective
        opt_eigen, _ = nctssos_first([p], x, 2, obj="eigen", QUIET=true)
        
        # Trace objective
        opt_trace, _ = nctssos_first([p], x, 2, obj="trace", QUIET=true)
        
        @test opt_eigen !== nothing
        @test opt_trace !== nothing
        @test opt_eigen <= opt_trace  # Eigen should be tighter or equal
    end

    @testset "Partitioned Variables" begin
        @ncpolyvar x[1:3]
        p = x[1]^2 + x[2]^2 + x[3]^2 - x[1]*x[2] + x[2]*x[1]
        
        # First step with partition
        opt, data = nctssos_first([p], x, 2, partition=2, QUIET=true)
        
        @test opt !== nothing
        @test data.partition == 2
    end
end