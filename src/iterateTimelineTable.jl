
struct IterateTimelineTable{T, MNames, FNames, N1, N2, CL}
    parent::MarketData{T, MNames, FNames, N1, N2}
    cols::CL
    key_vec::Vector{T}
    ranges::Vector{UnitRange{Int}}
    missing_vecs::Dict{T, Union{Nothing, SparseVector{Bool, Int}}}
    function IterateTimelineTable(
        data::MarketData{T, MNames, FNames, N1, N2},
        cols::CL,
        key_vec,
        ranges,
        missing_vecs
    ) where {T, MNames, FNames, N1, N2, CL}
        @assert Set(key_vec) ⊆ Set(keys(missing_vecs)) "Some Keys are Missing"
        @assert Set(key_vec) ⊆ Set(keys(data.firmdata)) "Some Keys are not in the Parent Data"
        new{T, MNames, FNames, N1, N2, CL}(data, cols, key_vec, ranges, missing_vecs)
    end
end

parent(data::IterateTimelineTable) = data.parent
iter_id(data::IterateTimelineTable) = data.key_vec
iter_range(data::IterateTimelineTable) = data.ranges
iter_missings(data::IterateTimelineTable) = data.missing_vecs
iter_cols(data::IterateTimelineTable) = data.cols

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
    parent(data)[id, r, iter_cols(data)..., missing_bdays=mssngs]
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

function Base.getindex(
    data::MarketData{T},
    ids::Vector{T},
    dates::ClosedInterval{Vector{Date}},
    cols::Union{Symbol, TimelineColumn, AbstractTerm}...
) where {T}
    u_ids = unique(ids)
    unique_cols = combine_columns(cols)
    rs = date_range(calendar(data), dates)
    range_limits = Dict(
        u_ids .=> [
            maximin(interval.(Ref(data), id, unique_cols)...) for id in u_ids
        ]
    )
    for i in eachindex(ids, rs)
        @inbounds rs[i] = maximin(rs[i], range_limits[ids[i]])
    end

    
    mssngs = Dict(
        u_ids .=>
        [
            combine_missing_bdays(
                get_missing_bdays.(
                    Ref(data),
                    id,
                    Ref(range_limits[id]),
                    unique_cols)...
                    ) for id in u_ids]
    )

    new_cols = convert_cols.(cols)
    IterateTimelineTable(
        data,
        new_cols,
        ids,
        rs,
        mssngs
    )
end

function Base.getindex(data::MarketData{T}, ids::Vector{T}, dates::ClosedInterval{Vector{Date}}, cols::FormulaTerm) where {T}
    sch = apply_schema(cols, schema(cols, data))
    data[ids, dates, sch.lhs, sch.rhs.terms...]
end

function Base.show(io::IO, data::IterateTimelineTable)
    println(io, "Iterable set of FixedTable with $(length(data)) unique datapoints")
    show(io, parent(data))
end