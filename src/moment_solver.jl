# C: is the varaibles commuting
# T: type of the coefficients
# order: order of the moment problem
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
struct MomentProblem{C,T} <: OptimizationProblem 
    order::Int
    model::GenericModel{T}
    monomap::Dict{Monomial{C},GenericVariableRef{T}}
end

function replace_variable_with_jump_variables(
    poly::Polynomial{C,T}, monomap::Dict{Monomial{C},GenericVariableRef{T2}}
) where {C,T,T2}
    return mapreduce(
        x -> coefficient(x) * monomap[remove_zero_degree(monomial(x))], +, terms(poly)
    )
end

function moment_relax(pop::PolynomialOptimizationProblem{C,T}, order::Int) where {C,T}
    # construct model
    objective = (symmetric_canonicalize ∘ (x -> getfield(x, :objective)))(pop)

    T2 = (T <: Integer) ? Float64 : T
    model = GenericModel{T2}()

    # total_basis that is the union of support of all constraints and objective and moment matrix and localizing matrices
    total_basis = begin
        moment_mtx_idcs = get_basis(pop.variables, order)
        sort(
            unique([
                remove_zero_degree(row_idx * col_idx) for
                row_idx in star.(moment_mtx_idcs) for col_idx in moment_mtx_idcs
            ]),
        )
    end

    monomap = init_moment_vector!(model, total_basis)

    constraint_matrices = [
        constrain_moment_matrix!(
            model,
            cur_poly,
            get_basis(
                pop.variables, order - ceil(Int, maxdegree(cur_poly) / 2)
            ),
            monomap,
        ) for cur_poly in vcat(pop.constraints..., one(objective))
    ]

    model[:mtx_constraints] = constraint_matrices

    @objective(
        model, Min, replace_variable_with_jump_variables(objective, monomap)
    )

    return MomentProblem(order, model, monomap)
end


function init_moment_vector!(
    model::GenericModel{T}, basis::Vector{Monomial{C}}
) where {C,T}
    moment_vector = @variable(model, y[1:length(basis)])
    @constraint(model, first(y) - 1 in Zeros())
    return Dict(basis .=> moment_vector)
end


function constrain_moment_matrix!(
    model::GenericModel{T},
    poly::Polynomial{C,T2},
    local_basis::Vector{Monomial{C}},
    basis2var_dict::Dict{Monomial{C},GenericVariableRef{T}},
) where {C,T,T2}
    moment_mtx = [
        replace_variable_with_jump_variables(poly * row_idx * col_idx, basis2var_dict) for
        row_idx in star.(local_basis), col_idx in local_basis
    ]
    return @constraint(
        model,
        moment_mtx in PSDCone(),
        base_name = isconstant(poly) ? "moment_matrix" : "localizing_matrix_$(hash(poly))"
    )
end
