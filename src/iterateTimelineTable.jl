
struct IterateOutput
    index::Int
    dates::ClosedInterval{Date}
end

iter_index(x::IterateOutput) = x.index
iter_dates(x::IterateOutput) = x.dates

struct IterateTimelineTable{T, MNames, FNames, N1, N2}
    parent::MarketData{T, MNames, FNames, N1, N2}
    index_dict::Dict{T, Vector{IterateOutput}}
    key_vec::Vector{T}
    function IterateTimelineTable(data::MarketData{T, MNames, FNames, N1, N2}, index_dict, key_vec) where {T, MNames, FNames, N1, N2}
        @assert Set(key_vec) ⊆ Set(keys(index_dict)) "Some Keys are Missing"
        @assert Set(keys(index_dict)) ⊆ Set(keys(data.firmdata)) "Some Keys are not in the Parent Data"
        new{T, MNames, FNames, N1, N2}(data, index_dict, key_vec)
    end
end

parent(data::IterateTimelineTable) = data.parent

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
    if state > length(iter.key_vec)
        return nothing
    end
    (iter[state], state+1)
end

function Base.eltype(::IterateTimelineTable{T, MNames, FNames, N1, N2}) where {T, MNames, FNames, N1, N2}
    Tuple{T, Vector{IterateOutput}}
end
function Base.length(iter::IterateTimelineTable)
    length(iter.key_vec)
end

total_length(data::IterateTimelineTable) = sum(length.(values(data.index_dict)))

function Base.getindex(data::IterateTimelineTable, i::Int)
    1 <= i <= length(data) || throw(BoundsError(data, i))
    id = data.key_vec[i]
    (id, data.index_dict[id])
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

function Base.getindex(data::MarketData{T}, ids::Vector{T}, dates::ClosedInterval{Vector{Date}}) where {T}
    id_dict = construct_id_dict(ids, dates)
    IterateTimelineTable(
        data,
        id_dict,
        unique(ids)
    )
end

function Base.show(io::IO, data::IterateTimelineTable)
    println(io, "Iterable set of TimelineTable with $(length(data)) unique firms")
    println(io, "and a total number of iterations of $(total_length(data)) and a parent of")
    show(io, parent(data))
end