#fermionic.jl
using SparseArrays

function operators(n)
    basis, vacio = integer_digits(n)
    l = length(basis)
    row = []
    col = []
    data = []
    for k in 1:n
        for i in 1:l
            cop = deepcopy(basis)
            if cop[i][k] != 0
                for j in 1:l
                    cop[i][k] = 0
                    if basis[j] == cop[i]
                        sign = (-1)^(sum(basis[i][1:k])+1)
                        append!(row, j)
                        append!(col, i+(k-1)*l)
                        append!(data, sign)
                    end
                end
            end
        end
    end
    append!(row, l)
    append!(col, 1)
    append!(data, 0)
    cm_tot = sparse(row, col, data)
    cd_tot = sparse(cm_tot')
    return cm_tot, cd_tot, l
end

function integer_digits(n)
    basis = []
    for i in 0:(2^n-1)
        binary_base = bitstring(i)[65-n:64]
        bin_vector = splitter(binary_base)
        append!(basis, [bin_vector])
    end
    vacio = [0 for i in 1:length(basis)]
    return basis, vacio
end


function splitter(x)
    vectorized = []
    for i in x
        str_to_int = parse(Int, i)
        append!(vectorized, str_to_int)
    end
    return vectorized
end

#=
------------------- Operators -------------------------
=#

@time begin
n = parse(Int,input("Choose the dimension of the system: "))
cm_tot, cd_tot, l = operators(n)
end


function cm(i, cm_tot = cm_tot, l = l)
    return cm_tot[1:l,((i-1)*l+1):i*l]
end

function cdm(i, cd_tot = cd_tot, l = l)
     #name cd is already taken
    return cd_tot[((i-1)*l+1):i*l,1:l]
end

function cdcm(i,j)
    return cdm(i)*cm(j)
end

function cmcd(i,j)
    return cm(i)*cdm(j)
end

function cmcm(i,j)
    return cm(i)*cm(j)
end

function cdcd(i,j)
    return cdm(i)*cdm(j)
end
