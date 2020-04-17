#fermionic.jl

module Fermionic_module
using SparseArrays

function operators(n)
    rowb, colb, basis, vacio = integer_digits(n)
    l = 2^n
    lb = n*2^(n-1)
    row = []
    col = []
    data = []
    for i in 1:lb
        j = rowb[i]-2^(n-colb[i])
        sign = (-1)^(sum(basis[rowb[i],1:colb[i]])+1)
        append!(row, j)
        append!(col, rowb[i]+(colb[i]-1)*l)
        append!(data, sign)
    end
    append!(row, l)
    append!(col, 1)
    append!(data, 0)
    cm_tot = sparse(row, col, data)
    cd_tot = sparse(cm_tot')
    return cm_tot, cd_tot, l
end

function integer_digits(n)
    rowb = []
    colb = []
    data = []
    for i in 0:(2^n-1)
        binary_base = bitstring(i)[65-n:64]
        bin_vector = splitter(binary_base)
        for j in 1:n
            if bin_vector[j] == 1
                append!(colb,j)
                append!(data,1)
                append!(rowb,i+1)
            end
        end
    end
    basis = sparse(rowb, colb, data)
    vacio = [0 for i in 1:2^n]
    return rowb, colb, basis, vacio
end


function splitter(x)
    vectorized = []
    for i in x
        str_to_int = parse(Int, i)
        append!(vectorized, str_to_int)
    end
    return vectorized
end

end
