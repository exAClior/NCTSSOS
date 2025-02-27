using DynamicPolynomials, NCTSSOS

n = 10
supp = [UInt16[]]
coe = Float64[n]
for i = 2:n
    push!(supp, UInt16[i-1;i-1;i-1;i-1], UInt16[i-1;i-1;i], UInt16[i], UInt16[i;i])
    push!(coe, 100, -200, -2, 101)
end

@time begin
opt,data = nctssos_first(supp, coe, n, newton=true, reducebasis=false, QUIET=true, solve=true, TS=false, obj="eigen")
end
@time begin
opt,data = cs_nctssos_first([supp], [coe], n, 2, QUIET=true, TS=false, obj="eigen")
end

@time begin
opt,data = nctssos_first([supp], [coe], n, 2, reducebasis=false, TS="MD", obj="eigen")
end

@time begin
opt,data = nctssos_higher!(data, TS="MD")
end