#fermionic.jl
using SparseArrays
using LinearAlgebra

function operators(n)
    rowb, colb, base, vacio = integer_digits(n)
    l = 2^n
    lb = n*2^(n-1)
    row = spzeros(lb+1)
    col = spzeros(lb+1)
    data = spzeros(lb+1)
    row[lb+1] = l
    col[lb+1] = 1
    data[lb+1] = 0
    for i in 1:lb
        j = floor(Int, rowb[i])-2^(n-floor(Int, colb[i]))
        sign = (-1)^(sum(base[floor(Int, rowb[i]),1:floor(Int, colb[i])])+1)
        row[i] = j
        col[i] = floor(Int, rowb[i])+(floor(Int, colb[i])-1)*l
        data[i] = sign
    end
    cm_tot = sparse(row, col, data)
    cd_tot = sparse(cm_tot')
    return cm_tot, cd_tot, l, base
end

function integer_digits(n)
    #rowb = []
    #colb = []
    #data = []
    rowb = spzeros(n*2^(n-1))
    colb = spzeros(n*2^(n-1))
    data = spzeros(n*2^(n-1))
    counter = 1
    for i in 0:(2^n-1)
        binary_base = bitstring(i)[65-n:64]
        bin_vector = splitter(binary_base)
        for j in 1:n
            if bin_vector[j] == 1
                rowb[counter] = (i+1)::Int
                colb[counter] = j::Int
                data[counter] = 1
                counter = counter + 1
                #append!(colb,j)
                #append!(data,1)
                #append!(rowb,i+1)
            end
        end
    end
    base = sparse(rowb, colb, data)
    vacio = spzeros(2^n)
    vacio[1] = 1
    return rowb, colb, base, vacio
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


#in order to define an operator of dim 4
#you must write op4 = Op(4)

struct Op
    dim::Int
    cmtot::SparseMatrixCSC{Float64,Int64}
    cdtot::SparseMatrixCSC{Float64,Int64}
    le::Int
    basis::SparseMatrixCSC{Float64,Int64}
    Op(dim) = new(dim, operators(dim)...)
end

#after being defined, one can
#retrieve the dimension by typping
# dim(op4)
dim(o::Op) = o.dim
cmtot(o::Op) = o.cmtot
cdtot(o::Op) = o.cdtot
le(o::Op) = o.le
basis(o::Op) = o.basis

#These are the fermionic operators we wanted
#By calling cm(op, 1) we get the sparse matrix
#corresponding to the destruction of the
#first fermionic mode
cm(o::Op, i::Int) = cmtot(o)[1:le(o),((i-1)*le(o)+1):i*le(o)]
cdm(o::Op, i::Int) = cdtot(o)[((i-1)*le(o)+1):i*le(o),1:le(o)]
cdcm(o::Op, i::Int, j::Int) = cdm(o,i)*cm(o,j)
cmcd(o::Op, i::Int, j::Int) = cm(o,i)*cdm(o,j)
cmcm(o::Op, i::Int, j::Int) = cm(o,i)*cm(o,j)
cdcd(o::Op, i::Int, j::Int) = cdm(o,i)*cdm(o,j)

#you can visualize the complete matrix instead of the sparse
#by writing Matrix(cm(op,n))

#=
------------------- States -------------------------
=#

#I am using the multi dispatch structure, so states can be
#defines both as vectors or as sparse vectors

struct State
    st::Array{Int64,1}
    ope::Op
end

struct State_sparse
    st::SparseVector{Float64,Int64}
    ope::Op
end

st(s::State) = s.st
st(s::State_sparse) = s.st

ope(s::State) = s.ope
ope(s::State_sparse) = s.ope

function rhospv(sta::State)
    n = dim(ope(sta))
    rhospv = zeros(n,n)
    for i in 1:n
        for j in 1:n
            rhospv[i,j] = st(sta)'*cdcm(ope(sta), i, j)*st(sta)
        end
    end
    return rhospv
end

function rhosps(sta::State_sparse)
    n = dim(ope(sta))
    rhosps = spzeros(n,n)
    for i in 1:n
        for j in 1:n
            rhosps[i,j] = st(sta)'*cdcm(ope(sta), i, j)*st(sta)
        end
    end
    return rhosps
end

rhosp(s::State) = rhospv(s)
rhosp(s::State_sparse) = rhosps(s)

eigensp(s::State) = eigvals(rhosp(s))
#For some reason, I can not compute eigenvalues
#directly from sparse matrices
eigensp(s::State_sparse) = eigvals(Matrix(rhosp(s)))

function ssp(sta::State)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(eigen[i]) + (1 - eigen[i])*log(1-eigen[i]))
        end
    end
    return s
end

function ssp(sta::State_sparse)
    eigen = eigensp(sta)
    lene = length(eigen)
    s = 0
    for i in 1:lene
        if eigen[i] != 0 && eigen[i] != 1
            s = s - (eigen[i]*log(eigen[i]) + (1 - eigen[i])*log(1-eigen[i]))
        end
    end
    return s
end
