struct FixedWidthMatrix{N, T}
    data::NTuple{N, T}
    function FixedWidthMatrix(xs::T...) where {T}
        @assert all(length.(xs) .== length(xs[1])) "Not all the same length"
        new{length(xs), T}(xs)
    end
    function FixedWidthMatrix(xs::NTuple{N, T}) where {N, T}
        @assert all(length.(xs) .== length(xs[1])) "Not all the same length"
        new{N, T}(xs)
    end
end

function Base.view(data::NTuple{N, DataVector}, dates::ClosedInterval{Date}) where {N}
    FixedWidthMatrix(view.(data, Ref(data)))
end

function Base.view(data::NTuple{N, DataVector}, dates::ClosedInterval{Date}, mssngs::SparseVector{Bool, Int}) where {N}
    FixedWidthMatrix(view.(data, Ref(data), Ref(mssngs)))
end

function Base.getindex(data::FixedWidthMatrix{N}, i::Int) where {N}
    @assert i <= N "Out of bounds error"
    data.data[i]
end

function Base.getindex(data::FixedWidthMatrix{N}, i::Int, j::Int) where {N}
    data.data[j][i]
end

Base.length(data::FixedWidthMatrix{N}) where {N} = N * length(data[1])
Base.size(data::FixedWidthMatrix{N}) where {N} = (length(data[1]), N)
Base.size(data::Adjoint{X, FixedWidthMatrix{N}}) where {X, N} = (N, data.parent[1])

LinearAlgebra.adjoint(x::FixedWidthMatrix) = Adjoint{eltype(x[1]), typeof(x)}(x)

# These are created to minimize the amount of allocations Julia does
# Julia typically allocates a vector for each loop, which when using
# so many loops, can create real garbage collection problems
# As it turns out, doing sum(abs2, resp - mean(resp)) also does
# an allocation, which could mean allocating a huge amount
# caluclating rss was even worse, so these functions are only
# meant to be used internally but do not allocate if passed a view
function calc_tss(resp::AbstractVector)
    out = 0.0
    m = mean(resp)
    @simd for x in resp
        out += (x - m) ^ 2
    end
    out
end
function point_pred(pred::FixedWidthMatrix{N}, coef::SVector{N}, i::Int) where {N}
    out = 0.0
    @simd for j in 1:N
        @inbounds out += pred[i, j] * coef[j]
    end
    out
end
function calc_rss(resp::AbstractVector, pred::FixedWidthMatrix{N}, coef::SVector{N}) where {N}
    @assert length(resp) == size(pred)[1] "Response is not same length as prediction matrix"
    out = 0.0
    @simd for i in 1:length(resp)
        @inbounds out += (resp[i] - point_pred(pred, coef, i)) ^ 2
    end
    out
end


function mult_add(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    @assert length(x) == length(y) "Vectors are not the same length"
    out = zero(T)
    @simd for i in eachindex(x, y)
        @inbounds out += x[i] * y[i]
    end
    out
end


function mult_square(x::FixedWidthMatrix{N, <:AbstractVector{T}}) where {N, T}
    out = MMatrix{N, N, T} |> zeros
    for i in 1:N
        for j in 1:i
            out[i, j] = mult_add(
                x[i],
                x[j]
            )
            if j != i
                out[j, i] = out[i, j]
            end
        end
    end
    SMatrix(out)
end

function Base.:(*)(x::Adjoint{T, <:FixedWidthMatrix{N1, <:AbstractVector{T}}}, y::FixedWidthMatrix{N2, <:AbstractVector{T}}) where {T, N1, N2}
    if x.parent === y
        return mult_square(y)
    end
    out = MMatrix{N1, N2, T} |> zeros
    for i in 1:N1
        for j in 1:N2
            out[i, j] = mult_add(x.parent[i], y[j])
        end
    end
    SMatrix(out)
end

function Base.:(*)(x::Adjoint{T, <:FixedWidthMatrix{N, <:AbstractVector{T}}}, y::FixedWidthMatrix{1, <:AbstractVector{T}}) where {T, N}
    out = MVector{N, T} |> zeros
    for i in 1:N
        out[i] = mult_add(x.parent[i], y[1])
    end
    SVector(out)
end

function Base.:(*)(x::Adjoint{T, <:FixedWidthMatrix{N, <:AbstractVector{T}}}, y::AbstractVector{T}) where {T, N}
    out = MVector{N, T} |> zeros
    for i in 1:N
        out[i] = mult_add(x.parent[i], y)
    end
    SVector(out)
end