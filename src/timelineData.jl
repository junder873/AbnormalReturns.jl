"""
```@setup general
data_dir = joinpath("..", "..", "test", "data") # hide
using CSV, DataFramesMeta, Dates, AbnormalReturns

df_firm = CSV.File(joinpath(data_dir, "daily_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath(data_dir, "mkt_ret.csv")) |> DataFrame
df_mkt[!, :mkt] = df_mkt.mktrf .+ df_mkt.rf
df_events = CSV.File(joinpath(data_dir, "firm_earnings_announcements.csv")) |> DataFrame
mkt_data = MarketData(
    df_mkt,
    df_firm
)
```
"""

"""
    struct DataVector
        data::OffsetVector{Float64, Vector{Float64}}
        missing_bdays::Union{Nothing, OffsetVector{Bool, SparseVector{Bool, Int}}}
        interval::UnitRange{Int}
        function DataVector(data, missing_bdays, interval)
            if missing_bdays !== nothing
                @assert length(data) == length(missing_bdays) "Data does not match length of dates or missings"
                @assert data.offsets == missing_bdays.offsets
            end
            @assert length(interval) == length(data)
            new(data, missing_bdays, interval)
        end
    end

    DataVector(data::AbstractVector{T}, offset::Int) where {T}

    DataVector(data::AbstractVector, d::Date, hc::MarketCalendar)
"""
struct DataVector
    data::OffsetVector{Float64, Vector{Float64}}
    missing_bdays::Union{Nothing, OffsetVector{Bool, SparseVector{Bool, Int}}}
    interval::UnitRange{Int}
    function DataVector(data, missing_bdays, interval)
        if missing_bdays !== nothing
            @assert length(data) == length(missing_bdays) "Data does not match length of dates or missings"
            @assert data.offsets == missing_bdays.offsets
        end
        @assert length(interval) == length(data)
        new(data, missing_bdays, interval)
    end
end

struct MarketData{T}
    calendar::MarketCalendar
    marketdata::Dict{Symbol, DataVector} # column names as symbols
    firmdata::Dict{T, Dict{Symbol, DataVector}} # data stored by firm id and then by column name as symbol
end

struct AllowMissing{mssng} end

"""
    struct FixedTable{N, T<:AbstractFloat, AV <: AbstractVector{T}, CL<:Union{Symbol, String}} <: AbstractMatrix{T}
        data::SVector{N, AV}
        cols::SVector{N, CL}
        req_length::Int
        function FixedTable(xs::SVector{N, AV}, cols::SVector{N, CL}, req_length=0) where {T, N, AV<:AbstractVector{T}, CL}
            new{N, T, AV, CL}(xs, cols, req_length)
        end
    end


This provides a fixed-width interface that is designed to allow quick access
(either through accessing a column number `data[:, 1]` or accessing a column
name `data[:, :a]`). `req_length` is an optional parameter that specifies the
length that a user originally requested, which is used in later functions
to determine if the FixedTable has too few rows.
"""
struct FixedTable{N, T<:AbstractFloat, AV <: AbstractVector{T}, CL<:Union{Symbol, String}} <: AbstractMatrix{T}
    data::SVector{N, AV}
    cols::SVector{N, CL}
    req_length::Int
    function FixedTable(xs::SVector{N, AV}, cols::SVector{N, CL}, req_length=0) where {T, N, AV<:AbstractVector{T}, CL}
        new{N, T, AV, CL}(xs, cols, req_length)
    end
end

function all_equal_length(xs)
    l1 = length(xs[1])
    for x in xs[2:end]
        if length(x) != l1
            return false
        end
    end
    true
end

#Base.length(data::FixedTable{N}) where {N} = N * length(data[1])
Base.size(data::FixedTable{N}) where {N} = (length(data[:, 1]), N)
Base.size(data::Adjoint{X, FixedTable{N}}) where {X, N} = (N, data.parent[1])

LinearAlgebra.adjoint(x::FixedTable{N, T}) where {N, T} = Adjoint{T, typeof(x)}(x)

"""
Checks whether each firm_id-date pair is unique, assumes that vectors are sorted by firm_id then date

Returns true if there is at least one firm_id-date pair repeated, false if all are unique
"""
function all_unique_obs(firm_ids::AbstractVector, dates::AbstractVector)
    @assert length(firm_ids) == length(dates) "Length of vectors are not the same"
    for i in 2:length(firm_ids)
        @inbounds if firm_ids[i] == firm_ids[i-1] && dates[i] == dates[i-1]
            return true
        end
    end
    return false
end

"""
    function MarketData(
        df_market,
        df_firms;
        date_col_market=:date,
        date_col_firms=:date,
        id_col=:permno,
        add_intercept_col=true,
        valuecols_market=nothing,
        valuecols_firms=nothing
    )

## Arguments
- `df_market`: A Tables.jl compatible source that stores market data, indexed by date.
    The dates must be a unique set. The column name for the date column is specified
    by the keyword argument "date_col_market"
- `df_firms`: A Tables.jl compatible source that stores firm data. Each firm must have a
    unique set of dates. The column name for the date column is specified
    by the keyword argument "date_col_firms" and the firm ID column is specified
    by the keyword argument "id_col"
- `valuecols_market=nothing`: If left as nothing, all other columns in `df_market` are
    used as the value columns. These are the columns that are stored in the resulting
    dataset. Otherwise a vector of Symbol or String specifying column names.
- `valuecols_firms=nothing`: Same as above
- `id_col=:permno`: The column corresponding to the set of firm IDs in `df_firms`
- `add_intercept_col=true`: Whether to add a column to the data for an intercept (which is
    always equal to 1)

MarketData is the main data storage structure. Data is stored for each firm in
a Dict, where the data itself is a NamedTuple (names corresponding to column names,
such as "ret"), and the keys for the Dict corresponding to firm IDs. The MarketData
struct also stores overall market data and a calendar of dates.

Any firm data must have a corresponding market data date, so there cannot be a
firm return if there is not a market return on that date.

## Example

```@example general
df_firm = CSV.File(joinpath(data_dir, "daily_ret.csv"))
df_mkt = CSV.File(joinpath(data_dir, "mkt_ret.csv"))
df_mkt[!, :mkt] = df_mkt.mktrf .+ df_mkt.rf
mkt_data = MarketData(
    df_mkt,
    df_firm
)
```
"""
function MarketData(
    df_market,
    df_firms;
    date_col_market=:date,
    date_col_firms=:date,
    id_col=:permno,
    add_intercept_col=true,
    valuecols_market=nothing,
    valuecols_firms=nothing
)
    df_market = DataFrame(df_market)
    df_firms = DataFrame(df_firms)
    if add_intercept_col
        df_market[!, :intercept] .= 1.0
        if valuecols_market !== nothing && :intercept ∉ valuecols_market
            push!(valuecols_market, :intercept)
        end
    end
    if valuecols_market === nothing
        valuecols_market = Symbol.([n for n in Symbol.(names(df_market)) if n ∉ [date_col_market]])
    end
    if valuecols_firms === nothing
        valuecols_firms = Symbol.([n for n in Symbol.(names(df_firms)) if n ∉ [date_col_firms, id_col]])
    end

    select!(df_market, vcat([date_col_market], valuecols_market))
    dropmissing!(df_market, date_col_market)
    #dropmissing!(df_market)
    sort!(df_market)

    if any(nonunique(df_market, [date_col_market]))
        @error("There are duplicate date rows in the market data")
    end

    select!(df_firms, vcat([id_col, date_col_firms], valuecols_firms))
    dropmissing!(df_firms, [id_col, date_col_firms])
    sort!(df_firms, [id_col, date_col_firms])

    # since this is sorted, just a single iteration is enough to check
    if all_unique_obs(df_firms[:, id_col], df_firms[:, date_col_firms])
        @error("There are duplicate id-date rows in the firm data")
    end

    cal = MarketCalendar(df_market[:, date_col_market])

    market_data = Dict(
        valuecols_market .=>
        DataVector.(
            Tables.columns(df_market[:, valuecols_market]),
            0
            # Ref(df_market[:, date_col_market]),
            # Ref(cal)
        )
    )

    check_all_businessdays(unique(df_firms[:, date_col_firms]), cal)

    gdf = groupby(df_firms, id_col)
    df_temp = combine(
        gdf,
        date_col_firms => (x -> add_missing_bdays(x, cal)) => date_col_firms
    )
    if nrow(df_temp) > 0
        insertcols!(df_temp, [col => missing for col in valuecols_firms]...)
        df_firms = vcat(
            df_firms,
            df_temp
        )
        sort!(df_firms, [id_col, date_col_firms])
        gdf = groupby(df_firms, id_col)
    end

    col_tab = columntable(gdf)


    firm_data = Dict{typeof(col_tab[id_col][1][1]), Dict{Symbol, DataVector}}()
    sizehint!(firm_data, length(col_tab[id_col]))

    for i in 1:length(col_tab[id_col])
        firm_data[col_tab[id_col][i][1]] = Dict(
            valuecols_firms .=> (
                DataVector(
                    col_tab[col][i],
                    col_tab[date_col_firms][i][1],
                    cal
                ) for col in valuecols_firms
            )
        )
    end

    MarketData(
        cal,
        market_data,
        firm_data
    )
end

function check_all_businessdays(dates, cal)
    bday_list = isbday(cal, dates)
    if !all(bday_list)
        bday_list_inv = (!).(bday_list)
        if sum(bday_list_inv) <= 3
            @error("Dates $(dates[bday_list_inv]) are not in the MARKET_DATA_CACHE")
        else
            d1_test = findfirst(bday_list_inv)
            d2_test = findlast(bday_list_inv)
            @error("Dates $(dates[d1_test]) ... $(dates[d2_test]) are not in the MARKET_DATA_CACHE")
        end
    end
end

function setup_calendar_data(f, data, offset)
    if any(ismissing.(data))
        missing_days = OffsetVector(sparsevec(any(ismissing, data, dims=2)), offset)
    else
        missing_days = nothing#OffsetVector(spzeros(Bool, size(data, 1)), offsets)
    end
    data = OffsetVector(coalesce.(data, zero(nonmissingtype(eltype(data)))), offset)
    f(
        data,
        missing_days,
        offset+1:offset+length(data)
    )
end

function DataVector(data::AbstractVector{T}, offset::Int) where {T}
    if any(ismissing.(data))
        if all(ismissing.(data))
            return DataVector(
                OffsetVector(zeros(nonmissingtype(T), 1), offset),
                OffsetVector(sparsevec([1], [true]), offset),
                offset+1:offset+1
            )
        end
        i = findfirst(!ismissing, data)
        j = findlast(!ismissing, data)
        setup_calendar_data(
            DataVector,
            data[i:j],
            offset+i-1
        )
    else
        setup_calendar_data(
            DataVector,
            data,
            offset
        )
    end
end

function DataVector(data::AbstractVector, d::Date, hc::MarketCalendar)
    DataVector(data, date_pos(hc, d)-1)
end

function add_missing_bdays(dates, cal)
    out = Date[]
    if length(dates) == bdayscount(cal, dates[1], dates[end]) + 1
        return out
    end
    dates_counter = 1
    for d in listbdays(cal, dates[1], dates[end])
        if dates[dates_counter] > d
            push!(out, d)
        else
            dates_counter += 1
        end
    end
    out
end

function get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::Union{Symbol, String, ContinuousTerm, Term}) where {T}
    mssngs = data_missing_bdays(data[id, col])
    if mssngs === nothing
        nothing
    else
        temp = mssngs[r]
        if nnz(temp) == 0
            nothing
        else
            temp
        end
    end
end

get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::InterceptTerm{true}) where {T} = get_missing_bdays(data, id, r, :intercept)

function get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::FunctionTerm) where {T}
    combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(r), col.args)...)
end
function get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::InteractionTerm) where {T}
    combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(r), col.terms)...)
end
function get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::StatsModels.LeadLagTerm{<:Any, typeof(lead)}) where {T}
    mssngs = data_missing_bdays(data[id, col.term])
    if mssngs === nothing
        nothing
    else
        temp = mssngs[r .+ col.nsteps]
        if nnz(temp) == 0
            nothing
        else
            temp
        end
    end
end
function get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::StatsModels.LeadLagTerm{<:Any, typeof(lag)}) where {T}
    mssngs = data_missing_bdays(data[id, col.term])
    if mssngs === nothing
        nothing
    else
        temp = mssngs[r .- col.nsteps]
        if nnz(temp) == 0
            nothing
        else
            temp
        end
    end
end

combine_missing_bdays(vals::Nothing...) = nothing
function combine_missing_bdays(vals::Union{Nothing, SparseVector{Bool, Int}}...)::SparseVector{Bool, Int}
    i = findfirst(!isnothing, vals)
    out = vals[i]
    for j in i+1:length(vals)
        if vals[j] !== nothing
            out = out .| vals[j]
        end
    end
    out
end

function Base.view(data::DataVector, r::UnitRange)
    view(data.data, r)
end

function Base.view(data::DataVector, r::AbstractVector)
    view(data.data, r)
end

function maximin(rs::UnitRange{Int}...)
    max(first.(rs)...):min(last.(rs)...)
end

##################################################
# Basic get functions that return a DataVector
##################################################

@inline function Base.getindex(
    data::MarketData{T},
    id::T,
    col::Symbol
) where {T}
    if col ∈ keys(data.marketdata)
        data.marketdata[col]
    else
        data.firmdata[id][col]
    end
end

Base.getindex(data::MarketData{T}, id::T, col::String) where {T} = data[id, Symbol(col)]
Base.getindex(data::MarketData{T}, id::T, col::Union{Term, ContinuousTerm}) where {T} = data[id, col.sym]
Base.getindex(data::MarketData{T}, id::T, col::StatsModels.LeadLagTerm) where {T} = data[id, col.term]

##################################################
# Simple function to return single value
##################################################
function Base.getindex(
    data::DataVector,
    loc::Int
)
    if loc ∉ interval(data)
        missing
    elseif data.missing_bdays !== nothing && data.missing_bdays[loc]
        missing
    else
        data.data[loc]
    end
end
function Base.getindex(
    data::MarketData{T},
    id::T,
    dt::Date,
    col::Union{AbstractString, Symbol}
) where {T}
    col = data[id, col]
    l = date_pos(calendar(data), dt)
    col[l]
end

##################################################
# Basic get functions that return a view of the
# underlying data
##################################################
function index_values(r, x)
    out = zeros(Int, length(x) - sum(x))
    j = 0
    for i in eachindex(r, x)
        if !x[i]
            j += 1
            @inbounds out[j] = r[i]
        end
        
    end
    out
end
Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::Symbol, missing_bdays::AbstractVector{Bool}) where {T} =
    view(data[id, col], index_values(r, missing_bdays))

Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::Symbol, ::Nothing=nothing) where {T} =
    view(data[id, col], r)

Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::String, missing_bdays=nothing) where {T} = 
    data[id, r, Symbol(col), missing_bdays]

Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::Union{Term, ContinuousTerm}, missing_bdays=nothing) where {T} =
    data[id, r, col.sym, missing_bdays]

Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::StatsModels.LeadLagTerm{<:Any, typeof(lead)}, missing_bdays=nothing) where {T} =
    data[id, r .+ col.nsteps, col.term, missing_bdays]


Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::StatsModels.LeadLagTerm{<:Any, typeof(lag)}, missing_bdays=nothing) where {T} =
    data[id, r .- col.nsteps, col.term, missing_bdays]

Base.getindex(data::MarketData{T}, id::T, r::UnitRange, col::InterceptTerm{true}, missing_bdays=nothing) where {T} =
    data[id, r, :intercept, missing_bdays]

##################################################
# Slightly more complex get functions to deal
# with interaction and function terms
# (these greatly slow down theprocess and are
# not recommended)
##################################################

@inline function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange,
    col::InteractionTerm,
    missing_bdays=nothing
) where {T}
    l = length(r)
    if missing_bdays !== nothing
        l -= count(missing_bdays)
    end
    out = fill(1.0, l)
    for t in col.terms
        out .*= data[id, r, t, missing_bdays]
    end
    out = OffsetVector(out, r[1]-1)
    if missing_bdays === nothing
        view(out, r)
    else
        view(out, r[1]:r[end]-count(missing_bdays))
    end
end

@inline function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange,
    col::FunctionTerm,
    missing_bdays=nothing
) where {T}
    out = OffsetVector(col.f.((data[id, r, t, missing_bdays] for t in col.args)...), r[1]-1)
    if missing_bdays === nothing
        view(out, r)
    else
        view(out, r[1]:r[end]-count(missing_bdays))
    end
end

##################################################
# Get functions that return a Fixed Table
# Starting with ones where the range is already
# defined (r), the assumption is this r
# is already checked and an error is returned
# if it is out of bounds anywhere
##################################################

@inline function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange{Int},
    cols::SVector{N, CL},
    missing_bdays::Union{Nothing, AbstractVector{Bool}};
    req_length::Int=0,
    col_names::SVector{N}=cols
) where {T, N, CL <: Union{String, Symbol}}
    if missing_bdays !== nothing
        @assert length(missing_bdays) == length(r) "Missing days is wrong length"
    end
    FixedTable(
        getindex.(Ref(data), id, Ref(r), cols, Ref(missing_bdays)),
        col_names,
        req_length
    )
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange{Int},
    cols::CL,
    missing_bdays::Union{Nothing, AbstractVector{Bool}};
    req_length::Int=0,
    col_names=nothing
) where {T, CL <: MatrixTerm}
    if missing_bdays !== nothing
        @assert length(missing_bdays) == length(r) "Missing days is wrong length"
    end
    N = length(cols.terms)
    FixedTable(
        SVector{N}(data[id, r, col, missing_bdays] for col in cols.terms),
        col_names === nothing ? SVector{N}(coefnames(cols.terms)) : col_names,
        req_length
    )
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange{Int},
    f::FormulaTerm,
    missing_bdays=nothing;
    req_length::Int=0,
    check_intercept=true
) where {T}
    f = adjust_formula(f; check_intercept)
    sch = apply_schema(f, schema(f, data))
    cols = MatrixTerm((sch.lhs, sch.rhs.terms...))
    data[
        id,
        r,
        cols,
        missing_bdays=missing_bdays,
        req_length=req_length
    ]
end

##################################################
# Get functions that return a Fixed Table
# these allow input for a date range which is
# converted to unitranges
##################################################

function Base.getindex(
    data::MarketData{T},
    id::T,
    dates::ClosedInterval{Date},
    cols::SVector{N};
) where {T, N}
    r1 = date_range(calendar(data), dates)
    r = maximin(r1, interval.(Ref(data), id, cols)...)
    mssngs = combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(r), cols)...)
    data[id, r, cols, mssngs, req_length=length(r1)]
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    dates::ClosedInterval{Date},
    cols::AbstractVector
) where {T}
    cols = SVector{length(cols)}(cols)
    data[id, dates, cols]
end


function Base.getindex(
    data::MarketData{T},
    id::T,
    dates::ClosedInterval{Date},
    f::FormulaTerm;
    check_intercept=true
) where {T}
    # this has to deal with schema here since in order to get missing bdays
    # and other items schema is necessary
    f = adjust_formula(f; check_intercept)
    sch = apply_schema(f, schema(f, data))
    out = (sch.lhs, sch.rhs.terms...)
    r1 = date_range(calendar(data), dates)
    r = maximin(r1, interval.(Ref(data), id, out)...)
    mssngs = combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(r), out)...)
    data[id, r, MatrixTerm(out), mssngs, req_length=length(r1)]
end

##################################################
# getindex that simply returns a tuple if columns
# are not provided
##################################################

Base.getindex(data::MarketData{T}, id::T, r::Union{UnitRange, ClosedInterval{Date}}) where {T} = (data, id, r)

##################################################
# functions related to FixedTable
##################################################

function Base.getindex(data::FixedTable{N}, ::Colon, i::Int) where {N}
    @assert 1 <= i <= N "Column not in data"
    @inbounds data.data[i]
end

function Base.getindex(data::FixedTable{N}, i::Int, j::Int) where {N}
    @assert 1 <= j <= N "Column not in data"
    data.data[j][i]
end

function Base.getindex(
    data::FixedTable{N, T, AV, CL},
    ::Colon,
    col::CL
) where {N, T, AV, CL}
    @assert col ∈ names(data) "Column is not in the data"
    i = findfirst(==(col), names(data))
    data[:, i]
end

function Base.getindex(data::FixedTable{N, T, AV, CL}, ::Colon, col) where {N, T, AV, CL}
    data[:, CL(col)]
end

Base.names(x::FixedTable) = x.cols




# function Base.length(data::TimelineTable{false})
#     real_dates = dates_min_max(data_dates(data), norm_dates(data))
#     c = bdayscount(calendar(data), dt_min(real_dates), dt_max(real_dates)) + isbday(calendar(data), dt_max(real_dates))
#     new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), real_dates)
#     return c - nnz(new_mssngs)
# end


function Base.length(x::DataVector)
    if missing_bdays === nothing
        length(raw_values(x))
    else
        length(data_missing_bdays(x)) - nnz(data_missing_bdays(x))
    end
end

# DataFrames.dropmissing(data::TimelineTable{false}) = data
# DataFrames.allowmissing(data::TimelineTable{true}) = data

# function DataFrames.dropmissing(data::TimelineTable{true})
#     TimelineTable(
#         parent(data),
#         data_id(data),
#         AllowMissing{false},
#         data_dates(data),
#         DictIndex(names(data)),
#         data_missing_bdays(data),
#         norm_dates(data)
#     )
# end

# function DataFrames.allowmissing(data::TimelineTable{false})
#     TimelineTable(
#         parent(data),
#         data_id(data),
#         AllowMissing{true},
#         data_dates(data),
#         DictIndex(names(data)),
#         data_missing_bdays(data),
#         norm_dates(data)
#     )
# end

raw_values(x::DataVector) = x.data
data_missing_bdays(x::DataVector) = x.missing_bdays


interval(x::DataVector) = x.interval
interval(data::MarketData{T}, id::T, col::Union{Symbol, String, ContinuousTerm, Term}) where {T} = interval(data[id, col])
function interval(data::MarketData{T}, id::T, col::FunctionTerm) where {T}
    maximin(interval.(Ref(data), id, col.args)...)
end
function interval(data::MarketData{T}, id::T, col::InteractionTerm) where {T}
    maximin(interval.(Ref(data), id, col.terms)...)
end
function interval(data::MarketData{T}, id::T, col::StatsModels.LeadLagTerm{<:Any, typeof(lead)}) where {T}
    r = interval(data[id, col.term])
    maximin(r, r .- col.nsteps)
end
function interval(data::MarketData{T}, id::T, col::StatsModels.LeadLagTerm{<:Any, typeof(lag)}) where {T}
    r = interval(data[id, col.term])
    maximin(r, r .+ col.nsteps)
end

interval(data::MarketData{T}, id::T, col::InterceptTerm{true}) where {T} = interval(data[id, :intercept])

calendar(x::MarketData) = x.calendar

data_dates(x::MarketData) = x.dates


function Base.show(io::IO, data::MarketData{T}) where {T}
    println(io, "MarketData with ID type $T with $(length(getfield(data, :firmdata))) unique firms")
    println(io, data.calendar)
    println(io, "Market Columns: $(join(keys(data.marketdata), ", "))")
    println(io, "Firm Columns: $(join(keys(data.firmdata[first(keys(data.firmdata))]), ", "))")
end

dt_min(x::ClosedInterval{Date}) = x.left
dt_max(x::ClosedInterval{Date}) = x.right


cal_dt_min(x::MarketData) = cal_dt_min(calendar(x))
cal_dt_max(x::MarketData) = cal_dt_max(calendar(x))

##################################
# transformation function to add new column to data
# this is not necessarily optimized, but is
# faster for columns that are in marketdata
##################################

in_market_data(data::MarketData, t::Symbol) = t ∈ keys(data.marketdata)
in_market_data(data::MarketData, t::ContinuousTerm) = in_market_data(data, t.sym)
in_market_data(data::MarketData, t::InteractionTerm) = all(in_market_data.(Ref(data), t.terms))
in_market_data(data::MarketData, t::StatsModels.LeadLagTerm) = in_market_data(data, t.term)
in_market_data(data::MarketData, t::FunctionTerm) = all(in_market_data.(Ref(data), t.args))

function DataFrames.transform!(
    data::MarketData,
    f::FormulaTerm
)
    sch = apply_schema(f, schema(f, data))
    coef_new = coefnames(sch.lhs) |> Symbol
    if all(in_market_data.(Ref(data), sch.rhs.terms))
        id = first(keys(data.firmdata))
        rs = maximin((interval(data, id, t) for t in sch.rhs.terms)...)
        mssngs = combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(rs), sch.rhs.terms)...)
        vals = fill(0.0, length(rs))
        for t in sch.rhs.terms
            vals = vals .+ data[id, rs, t, nothing] # no missings since the vector should equal 0 for those and the mssngs data will be for that
        end
        data.marketdata[coef_new] = DataVector(
            OffsetVector(vals, first(rs)-1),
            mssngs,
            rs
        )
    else
        for id in keys(data.firmdata)
            rs = maximin((interval(data, id, t) for t in sch.rhs.terms)...)
            mssngs = combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(rs), sch.rhs.terms)...)
            vals = fill(0.0, length(rs))
            for t in sch.rhs.terms
                vals = vals .+ data[id, rs, t, nothing] # no missings since the vector should equal 0 for those and the mssngs data will be for that
            end
            data.firmdata[id][coef_new] = DataVector(
                OffsetVector(vals, first(rs)-1),
                mssngs,
                rs
            )
        end
    end
end

function adjust_formula(f; check_intercept=true)
    if !check_intercept
        f
    elseif !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
        FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs);
    elseif StatsModels.omitsintercept(f)
        FormulaTerm(f.lhs, Tuple(r for r in f.rhs if !isa(r, ConstantTerm)))
    else
        f
    end
end