
struct IterateTimelineTable{T, MNames, FNames, N1, N2, NCol}
    parent::MarketData{T, MNames, FNames, N1, N2}
    cols::SVector{NCol, TimelineColumn}
    key_vec::Vector{T}
    ranges::Vector{UnitRange{Int}}
    missing_vecs::Dict{T, Union{Nothing, SparseVector{Bool, Int}}}
    function IterateTimelineTable(
        data::MarketData{T, MNames, FNames, N1, N2},
        cols::SVector{NCol, TimelineColumn},
        key_vec,
        ranges,
        missing_vecs
    ) where {T, MNames, FNames, N1, N2, NCol}
        @assert Set(key_vec) ⊆ Set(keys(missing_vecs)) "Some Keys are Missing"
        @assert Set(key_vec) ⊆ Set(keys(data.firmdata)) "Some Keys are not in the Parent Data"
        new{T, MNames, FNames, N1, N2, NCol}(data, cols, key_vec, ranges, missing_vecs)
    end
end

parent(data::IterateTimelineTable) = data.parent
iter_id(data::IterateTimelineTable) = data.key_vec
iter_range(data::IterateTimelineTable) = data.ranges
iter_missings(data::IterateTimelineTable) = data.missing_vecs
iter_cols(data::IterateTimelineTable) = data.cols

function construct_id_dict(
    ids::Vector{T},
    date_starts::Vector{Date},
    date_ends::Vector{Date},
) where {T}
    @assert length(ids) == length(date_starts) == length(date_ends) "All Vectors must be of same length"
    out = Dict{T, Vector{IterateOutput}}()
    sizehint!(out, length(ids))
    for (i, id) in enumerate(ids)
        if !haskey(out, id)
            out[id] = IterateOutput[]
        end
        push!(out[id], IterateOutput(i, date_starts[i] .. date_ends[i]))
    end
    out
end

construct_id_dict(ids::Vector{T}, dates::ClosedInterval{Vector{Date}}) where {T} = construct_id_dict(ids, dates.left, dates.right)

function Base.iterate(iter::IterateTimelineTable{T, MNames, FNames}, state=1) where {T, MNames, FNames}
    if state > length(iter)
        return nothing
    end
    (iter[state], state+1)
end

# function Base.eltype(::IterateTimelineTable{T, MNames, FNames, N1, N2}) where {T, MNames, FNames, N1, N2}
#     Tuple{T, Vector{IterateOutput}}
# end
function Base.length(iter::IterateTimelineTable)
    length(iter.key_vec)
end

function Base.getindex(data::IterateTimelineTable, i::Int)
    #1 <= i <= length(data) || throw(BoundsError(data, i))
    id = iter_id(data)[i]
    r = iter_range(data)[i]
    mssngs = iter_missings(data)[id]
    if mssngs !== nothing
        mssngs = mssngs[r]
    end
    parent(data)[id, r, iter_cols(data), mssngs]
end

Base.firstindex(data::IterateTimelineTable) = 1
Base.lastindex(data::IterateTimelineTable) = length(data) 

function validate_iterator(
    data::IterateTimelineTable,
    out_vector::Vector
)
    all_ids = iter_index.(vcat(values(data.index_dict)...)) |> sort
    all_ids == 1:length(out_vector)
end

function Base.getindex(data::MarketData{T}, ids::Vector{T}, dates::ClosedInterval{Vector{Date}}, cols) where {T}
    u_ids = unique(ids)
    rs = date_range(calendar(data), dates)
    mssngs = Dict(
        u_ids .=>
        [
            combine_missing_bdays(
                get_missing_bdays.(
                    Ref(data),
                    id,
                    Ref(maximin(interval.(Ref(data), id, cols)...)),
                    cols)...
                    ) for id in u_ids]
    )
    IterateTimelineTable(
        data,
        cols,
        ids,
        rs,
        mssngs
    )
end

function Base.show(io::IO, data::IterateTimelineTable)
    println(io, "Iterable set of FixedTable with $(length(data)) unique datapoints")
    show(io, parent(data))
end