function get_basis(n, d)
    lb = binomial(n+d, d)
    basis = zeros(UInt8, n, lb)
    i = 0
    t = 1
    while i < d+1
        t += 1
        if basis[n,t-1] == i
           if i < d
              basis[1,t] = i+1
           end
           i += 1
        else
            j = findfirst(x->basis[x,t-1]!=0, 1:n)
            basis[:,t] = basis[:,t-1]
            if j == 1
               basis[1,t] -= 1
               basis[2,t] += 1
            else
               basis[1,t] = basis[j,t] - 1
               basis[j,t] = 0
               basis[j+1,t] += 1
            end
        end
    end
    return basis
end

function _cyclic_canon(a::Vector{UInt16})
    if isempty(a)
        return a
    else
        return minimum([[a[i+1:length(a)]; a[1:i]] for i=0:length(a)-1])
    end
end

function cyclic_canon(supp, coe; type=Float64)
    nsupp = [min(_cyclic_canon(word), _cyclic_canon(reverse(word))) for word in supp]
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(type, l)
    for (i,item) in enumerate(supp)
        bi = min(_cyclic_canon(item), _cyclic_canon(reverse(item)))
        Locb = bfind(nsupp, l, bi)
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

function _sym_canon(a::Vector{UInt16})
    i = 1
    while i <= Int(ceil((length(a)-1)/2))
        if a[i] < a[end+1-i]
            return a
        elseif a[i] > a[end+1-i]
            return reverse(a)
        else
            i += 1
        end
    end
    return a
end

function is_sym(a::Vector{UInt16})
    l = Int(ceil((length(a)-1)/2))
    return isequal(a[1:l], a[end:-1:end-l+1])
end

function sym_canon(supp, coe; type=Float64)
    nsupp = copy(supp)
    nsupp = _sym_canon.(nsupp)
    sort!(nsupp)
    unique!(nsupp)
    l = length(nsupp)
    ncoe = zeros(type, l)
    for (i,item) in enumerate(supp)
        Locb = bfind(nsupp, l, _sym_canon(item))
        ncoe[Locb] += coe[i]
    end
    return nsupp,ncoe
end

function get_ncbasis(n, d; ind=Vector{UInt16}(1:n), binary=false)
    basis = [UInt16[]]
    for i = 1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind, binary=binary))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=Vector{UInt16}(1:n), binary=false)
    if d > 0
        basis = Vector{UInt16}[]
        for i = 1:n
            temp = _get_ncbasis_deg(n, d-1, ind=ind, binary=binary)
            if binary == false || d == 1
                push!.(temp, ind[i])
                append!(basis, temp)
            else
                for item in temp
                    if item[end] != ind[i]
                        push!(basis, [item;ind[i]])
                    end
                end
            end
        end
        return basis
    else
        return [UInt16[]]
    end
end

function constraint_reduce!(word::Vector{UInt16}; constraint="unipotent")
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
            if constraint == "unipotent"
                deleteat!(word, i)
                if i > 1
                     i -= 1
                end
            end
        else
            i += 1
        end
    end
    return word
end

function reduce!(word::Vector{UInt16}; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    if obj == "trace"
        word = min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
    else
        if constraint === nothing
            word = min(_comm(word, partition, comm_var), _comm(reverse(word), partition, comm_var))
        else
            word = min(constraint_reduce!(_comm(word, partition, comm_var), constraint = constraint), constraint_reduce!(_comm(reverse(word), partition, comm_var), constraint = constraint))
        end
    end
    return word
end

function reduce(w::Mono, x; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    n = length(x)
    ind = w.z .> 0
    vars = w.vars[ind]
    exp = w.z[ind]
    word = UInt16[]
    for j = 1:length(vars)
        k = bfind(x, n, vars[j], rev=true)
        append!(word, k*ones(UInt16, exp[j]))
    end
    word = reduce!(word, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
    return prod(x[word])
end

function _comm(word::Vector{UInt16}, partition, comm_var)
    if partition > 0
        w1 = copy(word[word .<= partition])
        w2 = word[word .> partition]
    else
        w1 = copy(word)
        w2 = UInt16[]
    end
    if comm_var > 0
        i = 1
        while i < length(w1)
            if w1[i] <= comm_var && w1[i+1] <= comm_var && w1[i] > w1[i+1]
                w1[i],w1[i+1] = w1[i+1],w1[i]
                if i > 1
                    i -= 1
                else
                    i = 2
                end
            else
                i += 1
            end
        end
    end
    return [w1; w2]
end

function bfind(A, l, a; lt=isless, rev=false)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif lt(A[mid], a)
            if rev == false
                low = mid+1
            else
                high = mid-1
            end
        else
            if rev == false
                high = mid-1
            else
                low = mid+1
            end
        end
    end
    return nothing
end

function permutation(a)
    b = sparse(a)
    ua = convert(Vector{UInt16}, b.nzind)
    na = convert(Vector{UInt16}, b.nzval)
    return _permutation(ua, na)
end

function _permutation(ua, na)
    if !isempty(ua)
        perm = Vector{UInt16}[]
        for i = 1:length(ua)
            nua = copy(ua)
            nna = copy(na)
            if na[i] == 1
                deleteat!(nua, i)
                deleteat!(nna, i)
            else
                nna[i] -= 1
            end
            temp = _permutation(nua, nna)
            push!.(temp, ua[i])
            append!(perm, temp)
        end
        return perm
    else
        return [UInt16[]]
    end
end

function polys_info(pop, x)
    n = length(x)
    m = length(pop)-1
    coe = Vector{Vector{Float64}}(undef, m+1)
    supp = Vector{Vector{Vector{UInt16}}}(undef, m+1)
    for k = 1:m+1
        mon = monomials(pop[k])
        coe[k] = coefficients(pop[k])
        supp[k] = [UInt16[] for i=1:length(mon)]
        for i = 1:length(mon)
            ind = mon[i].z .> 0
            vars = mon[i].vars[ind]
            exp = mon[i].z[ind]
            for j = 1:length(vars)
                l = bfind(x, n, vars[j], rev=true)
                append!(supp[k][i], l*ones(UInt16, exp[j]))
            end
        end
    end
    return n,supp,coe
end

function poly_info(f, x)
    n = length(x)
    mon = monomials(f)
    coe = coefficients(f)
    lm = length(mon)
    supp = [UInt16[] for i=1:lm]
    for i = 1:lm
        ind = mon[i].z .> 0
        vars = mon[i].vars[ind]
        exp = mon[i].z[ind]
        for j = 1:length(vars)
            k = bfind(x, n, vars[j], rev=true)
            append!(supp[i], k*ones(UInt16, exp[j]))
        end
    end
    return n,supp,coe
end

function isless_td(a, b)
    if length(a) < length(b)
        return true
    elseif length(a) > length(b)
        return false
    else
        return a < b
    end
end

function sym_cyclic(word)
    return min(_cyclic_canon(word), _cyclic_canon(reverse(word)))
end

function _cyclic_basis(var)
    basis = _permutation(var, ones(length(var)))
    return unique(_cyclic_canon.(basis))
end

function sym(word, vargroup)
    cword = copy(word)
    ind = [gind(word[i], vargroup) for i = 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i = 1:length(uind)]
    k = 0
    for i = 1:length(uind)
        cword[k+1:k+nind[i]] = reverse(cword[k+1:k+nind[i]])
        k += nind[i]
    end
    return min(word, cword)
end

function iscomm(a, vargroup)
    for i = 1:length(a)-1
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            return false
        end
    end
    return true
end

function gind(k, vargroup)
    return findfirst(i -> k <= sum(vargroup[1:i]), 1:length(vargroup))
end

function res_comm!(a, vargroup)
    i = 1
    while i < length(a)
        if a[i] > a[i+1] && gind(a[i], vargroup) != gind(a[i+1], vargroup)
            temp = a[i]
            a[i] = a[i+1]
            a[i+1] = temp
            if i > 1
                i -= 1
            end
        else
            i += 1
        end
    end
    return a
end

function issym(word, vargroup)
    ind = [gind(word[i], vargroup) for i = 1:length(word)]
    uind = unique(ind)
    nind = [count(ind .== uind[i]) for i = 1:length(uind)]
    k = 0
    for i = 1:length(uind)
        temp = word[k+1:k+nind[i]]
        if reverse(temp) != temp
            return false
        end
        k += nind[i]
    end
    return true
end

function star(w::Mono)
    return prod(reverse(w.vars).^reverse(w.z))
end

function star(p::Poly{T}) where T<:Union{Number, AffExpr}
    return coefficients(p)'*star.(monomials(p))
end

# generate an SOHS polynomial with variables vars and degree 2d
function add_SOHS!(model, vars, d; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i,item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(j -> vs[j] < vars[partition] && vs[j+1] >= vars[partition], 1:length(vs)-1) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    sohs = 0
    pos = @variable(model, [1:length(basis), 1:length(basis)], PSD)
    for j = 1:length(basis), k = j:length(basis)
        word = reduce(star(basis[j])*basis[k], vars, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint)
        if j == k
            @inbounds sohs += pos[j,k]*word
        else
            @inbounds sohs += 2*pos[j,k]*word
        end
    end
    return sohs
end

# generate a polynomial with variables vars and degree d
function add_poly!(model, vars, d; partition=0, constraint=nothing)
    basis = vcat([MultivariatePolynomials.monomials(vars, i) for i = 0:d]...)
    if constraint !== nothing
        basis = basis[[all(item.z .< 2) for item in basis]]
    end
    if partition > 0
        ind = Int[]
        for (i,item) in enumerate(basis)
            vs = item.vars[item.z .> 0]
            if findfirst(j -> vs[j] < vars[partition] && vs[j+1] >= vars[partition], 1:length(vs)-1) === nothing
                push!(ind, i)
            end
        end
        basis = basis[ind]
    end
    free = @variable(model, [1:length(basis)])
    poly = free'*basis
    return poly
end

function arrange(p, vars; obj="eigen", partition=0, comm_var=0, constraint=nothing)
    mons = monomials(p)
    coe = coefficients(p)
    mons = [reduce(mon, vars, obj=obj, partition=partition, comm_var=comm_var, constraint=constraint) for mon in mons]
    nmons = unique(sort(mons))
    ncoe = zeros(typeof(coe[1]), length(nmons))
    for (i,item) in enumerate(coe)
        Locb = bfind(nmons, length(nmons), mons[i])
        ncoe[Locb] += item
    end
    return nmons,ncoe
end
