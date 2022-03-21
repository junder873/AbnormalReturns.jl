# Iteration Util

# To try to speed up the vectorized regressions, this implements
# a simple iteration method. The idea is only to access the MarketData type
# as little as possible and use the update_dates! method to iterate
# through all the components of that firm id. This would mean the
# TimelineTable is only constructed once per firm id, so firms with multiple
# events will make this method faster relative to one firm per event


struct IterateMarketData{T, MNames, FNames, N1, N2}
    data::MarketData{T, MNames, FNames, N1, N2}
    index_dict::Dict{T, Vector{Tuple{Int, Date, Date}}}
    index_vec::Vector{T}
    function IterateMarketData(data::MarketData{T, MNames, FNames, N1, N2}, index_dict, index_vec) where {T, MNames, FNames, N1, N2}
        @assert Set(index_vec) ⊆ Set(keys(index_dict)) "Some Keys are Missing"
        @assert Set(keys(index_dict)) ⊆ Set(keys(data.firmdata)) "Some Keys are not in the Parent Data"
        new{T, MNames, FNames, N1, N2}(data, index_dict, index_vec)
    end
end
function construct_id_dict(ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    @assert length(ids) == length(date_starts) == length(date_ends) "All Vectors must be of same length"
    out = Dict{T, Vector{Tuple{Int, Date, Date}}}()
    sizehint!(out, length(ids))
    for (i, id) in enumerate(ids)
        if !haskey(out, id)
            out[id] = Tuple{Int, Date, Date}[]
        end
        push!(out[id], (i, date_starts[i], date_ends[i]))
    end
    out
end
function IterateMarketData(data::MarketData{T}, ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    out = construct_id_dict(ids, date_starts, date_ends)
    IterateMarketData(
        data,
        out,
        unique(ids)
    )
end

function Base.iterate(iter::IterateMarketData{T, MNames, FNames}, state=1) where {T, MNames, FNames}
    if state > length(iter.index_vec)
        return nothing
    end
    (iter[state], state+1)
end

function Base.eltype(::IterateMarketData{T, MNames, FNames, N1, N2}) where {T, MNames, FNames, N1, N2}
    Tuple{T, Vector{Tuple{Int, Date, Date}}}
end
function Base.length(iter::IterateMarketData)
    length(iter.index_vec)
end

function validate_iterator(
    data::Dict,
    out_vector::Vector
)
    all_ids = first.(vcat(values(data)...)) |> sort
    all_ids == 1:length(out_vector)
end

function Base.getindex(data::IterateMarketData, i::Int)
    1 <= i <= length(data) || throw(BoundsError(data, i))
    id = data.index_vec[i]
    (id, data.index_dict[id])
end

Base.firstindex(data::IterateMarketData) = 1
Base.lastindex(data::IterateMarketData) = length(data)