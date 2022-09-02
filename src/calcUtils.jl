

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
calc_tss(resp::FixedTable{1}) = calc_tss(resp[:, 1])
function point_pred(pred::FixedTable{N}, coef::SVector{N}, i::Int) where {N}
    out = 0.0
    @simd for j in 1:N
        @inbounds out += pred[i, j] * coef[j]
    end
    out
end
function calc_rss(resp::AbstractVector, pred::FixedTable{N}, coef::SVector{N}) where {N}
    @assert length(resp) == size(pred)[1] "Response is not same length as prediction matrix"
    out = 0.0
    @simd for i in 1:length(resp)
        @inbounds out += (resp[i] - point_pred(pred, coef, i)) ^ 2
    end
    out
end
calc_rss(resp::FixedTable{1}, pred::FixedTable{N}, coef::SVector{N}) where {N} = calc_rss(resp[:, 1], pred, coef)



function mult_add(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    @assert length(x) == length(y) "Vectors are not the same length"
    out = zero(T)
    @simd for i in eachindex(x, y)
        @inbounds out += x[i] * y[i]
    end
    out
end


function mult_square(x::FixedTable{N, T}) where {N, T}
    out = MMatrix{N, N, T} |> zeros
    for i in 1:N
        for j in 1:i
            out[i, j] = mult_add(
                x[:, i],
                x[:, j]
            )
            if j != i
                out[j, i] = out[i, j]
            end
        end
    end
    SMatrix(out)
end

function Base.:(*)(x::Adjoint{T, <:FixedTable{N1, T}}, y::FixedTable{N2, T}) where {T, N1, N2}
    if x.parent === y
        return mult_square(y)
    end
    out = MMatrix{N1, N2, T} |> zeros
    for i in 1:N1
        for j in 1:N2
            out[i, j] = mult_add(x.parent[:, i], y[:, j])
        end
    end
    SMatrix(out)
end

function Base.:(*)(x::Adjoint{T, <:FixedTable{N, T}}, y::FixedTable{1, T}) where {T, N}
    out = MVector{N, T} |> zeros
    for i in 1:N
        out[i] = mult_add(x.parent[:, i], y[:, 1])
    end
    SVector(out)
end

function Base.:(*)(x::Adjoint{T, <:FixedTable{N, T}}, y::AbstractVector{T}) where {T, N}
    out = MVector{N, T} |> zeros
    for i in 1:N
        out[i] = mult_add(x.parent[:, i], y)
    end
    SVector(out)
end

function resp_matrix(data::FixedTable{N, T, AV}) where {N, T, AV}
    @assert N >= 2 "Not enough columns"
    FixedTable(
        NTuple{N-1, AV}(data.data[2:end]),
        SVector{N-1}(names(data)[2:end])
    )
end