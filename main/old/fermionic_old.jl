#fermionic.jl
using SparseArrays


function operators(n)
    basis, vacio = integer_digits(n)
    l = length(basis)
    row = []
    col = []
    data = []
    for i in 1:n
        for j in 1:l
            if eltype(c(i,j,basis)) == Int64
                append!(row, c(i, j, basis)[1])
                append!(col, j+(i-1)*l)
                append!(data, c(i, j, basis)[2])
            end
        end
    end

    #method sparse is not recognizing the dimensions
    #and setting them as the maximum non zero values
    #I create a ficticious 0 entry for fixing This
    append!(row, l)
    append!(col, l*n)
    append!(data, 0)
    #cm_tot = sparse(row, col, data, [l,l*n,combine])
    cm_tot = sparse(row, col, data)
    cd_tot = sparse(cm_tot')
    return cm_tot, cd_tot, l
end


function c(mode, pos, basis)
    cop = deepcopy(basis)
    l = length(basis)
    if cop[pos][mode] != 0
        cop[pos][mode] = 0
        for i in 1:l
            if basis[i] == cop[pos]
                if mode == 0
                    sign = 1
                else
                    sign = (-1)^(sum(cop[pos][1:mode]))
                end
                return i, sign
            end
        end
    end
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
