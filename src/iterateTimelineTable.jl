
struct IterateFixedTable{T, N, CL<:Union{Symbol, String}, COL<:Union{MatrixTerm, SVector{N, Symbol}}}
    parent::MarketData{T}
    col_names::SVector{N, CL}
    cols::COL
    key_vec::Vector{T}
    ranges::Vector{UnitRange{Int}}
    missing_vecs::Dict{T, Union{Nothing, OffsetVector{Bool, SparseVector{Bool, Int}}}}
    req_lengths::Vector{Int}
    function IterateFixedTable(
        data::MarketData{T},
        col_names::SVector{N, CL},
        cols::COL,
        key_vec,
        ranges,
        missing_vecs,
        req_lengths=fill(0, length(key_vec))
    ) where {T, N, CL, COL}
        @assert Set(key_vec) ⊆ Set(keys(missing_vecs)) "Some Keys are Missing"
        @assert Set(key_vec) ⊆ Set(keys(data.firmdata)) "Some Keys are not in the Parent Data"
        @assert length(key_vec) == length(ranges) == length(req_lengths) "Vectors must be same length"
        new{T, N, CL, COL}(data, col_names, cols, key_vec, ranges, missing_vecs, req_lengths)
    end
end

parent(data::IterateFixedTable) = data.parent
iter_id(data::IterateFixedTable) = data.key_vec
iter_range(data::IterateFixedTable) = data.ranges
iter_missings(data::IterateFixedTable) = data.missing_vecs
iter_cols(data::IterateFixedTable) = data.cols
iter_col_names(data::IterateFixedTable) = data.col_names
iter_req_lengths(data::IterateFixedTable) = data.req_lengths

function Base.iterate(iter::IterateFixedTable{T}, state=1) where {T}
    if state > length(iter)
        return nothing
    end
    (iter[state], state+1)
end

# function Base.eltype(::IterateFixedTable{T, MNames, FNames, N1, N2}) where {T, MNames, FNames, N1, N2}
#     Tuple{T, Vector{IterateOutput}}
# end
function Base.length(iter::IterateFixedTable)
    length(iter.key_vec)
end

function Base.getindex(data::IterateFixedTable, i::Int)
    #1 <= i <= length(data) || throw(BoundsError(data, i))
    id = iter_id(data)[i]
    r = iter_range(data)[i]
    mssngs = iter_missings(data)[id]
    if mssngs !== nothing
        mssngs = mssngs[r]
        if nnz(mssngs) == 0
            mssngs = nothing
        end
    end
    parent(data)[id, r, iter_cols(data), mssngs, col_names=iter_col_names(data), req_length=iter_req_lengths(data)[i]]
end

Base.firstindex(data::IterateFixedTable) = 1
Base.lastindex(data::IterateFixedTable) = length(data) 

function joined_missings(data, id, range_limits, cols)
    r = range_limits[id]
    out = combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(r), cols)...)
    if out !== nothing
        OffsetVector(out, first(r)-1)
    else
        nothing
    end
end

function Base.getindex(
    data::MarketData{T},
    ids::Vector{T},
    dates::ClosedInterval{Vector{Date}},
    cols
) where {T}
    @assert length(ids) == length(dates.left) == length(dates.right) "Vectors are not the same length"
    u_ids = unique(ids)

    rs = date_range(calendar(data), dates)
    r_lengths = length.(rs)
    range_limits = Dict(
        u_ids .=> [
            maximin(interval.(Ref(data), id, cols)...) for id in u_ids
        ]
    )
    for i in eachindex(ids, rs)
        @inbounds rs[i] = maximin(rs[i], range_limits[ids[i]])
    end

    
    mssngs = Dict(u_ids .=> joined_missings.(Ref(data), u_ids, Ref(range_limits), Ref(cols)))

    if eltype(cols) <: AbstractTerm
        col_names = SVector{length(cols)}(coefnames(cols))
        final_cols = MatrixTerm(cols)
    else
        col_names = SVector{length(cols)}(cols)
        final_cols = SVector{length(cols)}(Symbol.(cols))
    end
    IterateFixedTable(
        data,
        col_names,
        final_cols,
        ids,
        rs,
        mssngs,
        r_lengths
    )
end

function Base.getindex(
    data::MarketData{T},
    ids::Vector{T},
    dates::ClosedInterval{Vector{Date}},
    f::FormulaTerm;
    check_intercept=true
) where {T}
    f = adjust_formula(f; check_intercept)
    sch = apply_schema(f, schema(f, data))
    out = (sch.lhs, sch.rhs.terms...)
    data[ids, dates, out]
end


function Base.show(io::IO, data::IterateFixedTable)
    println(io, "Iterable set of FixedTable with $(length(data)) unique datapoints")
    show(io, parent(data))
end

function StatsModels.schema(f::FormulaTerm, d::MarketData)
    StatsModels.Schema(
        Dict(
            term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0)
        )
    )
end

function StatsModels.schema(f::FormulaTerm, d::FixedTable)
    StatsModels.Schema(
        Dict(
            term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0)
        )
    )
end

function Base.getindex(
    data::MarketData{T},
    ids::Vector{T},
    dates::ClosedInterval{Vector{Date}},
) where {T}
    (data, ids, dates)
end