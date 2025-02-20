using OMEinsumContractionOrders, TreeWidthSolver
using OMEinsumContractionOrders: adjacency_matrix
using OMEinsumContractionOrders: NestedEinsum
using OMEinsum

function get_cliques(tree::DecompositionTreeNode{T}) where T<:Integer
	isempty(tree.children) && return [collect(tree.bag)]
	return vcat([collect(tree.bag)], [get_cliques(child) for child in tree.children ]...)
end

# copied from https://github.com/ArrogantGao/TensorBranching.jl/blob/782826c657eb9158b82d1086b98691cc19a3dc8a/src/tree_decompose.jl#L32


function graph2eincode(g)
    ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
    return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
end

function tree_decomposition(g, method)
    code_interm = graph2eincode(g)
    if method == :sa
        ordered_eincode = OMEinsumContractionOrders.optimize_sa(
        code_interm,
        uniformsize(code_interm, 2);
        sc_target=10,
        βs=0.1:0.2:20.0,
        ntrials=20,
        initializer=:greedy,
    )
    elseif method == :kahypar
        ordered_eincode = OMEinsumContractionOrders.optimize_kahypar(code_interm,uniformsize(code_interm, 2); max_group_size=50, sc_target=10)
	elseif method == :exact
        ordered_eincode = OMEinsumContractionOrders.optimize_exact_treewidth(
            ExactTreewidth(), code_interm, uniformsize(code_interm, 2)
        )
    end

    elim_order = reverse!(label_elimination_order(ordered_eincode)) 

    return decomposition_tree(g, elim_order)
end