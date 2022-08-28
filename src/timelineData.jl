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
    struct DataVector <: CalendarData
        data::Vector{Float64}
        missing_bdays::SparseVector{Bool, Int}
        dates::ClosedInterval{Date}
        calendar::MarketCalendar
        function DataVector(data, missing_bdays, dates, calendar)
            @assert length(data) == length(missing_bdays) == bdayscount(calendar, dates.left, dates.right) + 1 "Data does not match length of dates or missings"
            new(data, missing_bdays, dates, calendar)
        end
    end

    DataVector(data::AbstractVector, dates::ClosedInterval{Date}, cal::MarketCalendar)

    DataVector(data::AbstractVector, dates::AbstractVector{Date}, cal::MarketCalendar)
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

struct MarketData{T, MNames, FNames, N1, N2}
    calendar::MarketCalendar
    marketdata::NamedTuple{MNames, NTuple{N1, DataVector}} # column names as symbols
    firmdata::Dict{T, NamedTuple{FNames, NTuple{N2, DataVector}}} # data stored by firm id and then by column name as symbol
end

struct AllowMissing{mssng} end

"""
    mutable struct TimelineTable{Mssng, T, MNames, FNames, N1, N2} <: Tables.AbstractColumns
        parent::MarketData{T, MNames, FNames, N1, N2}
        "The current ID"
        id::T
        "Whether the functions return Vector{Union{Missing, Float64}} or Vector{Float64}"
        allow_missing::Type{AllowMissing{Mssng}}
        "Actual returnable dates that guarentees a square matrix made of the minimum and maximum of the DataVectors for this ID and set of columns"
        dates::ClosedInterval{Date}
        "Index of column names and any lags/leads"
        cols::DictIndex
        "missing days between the dates"
        missing_bdays::SparseVector{Bool, Int}
        "dates requested, which provides what data is automatically returned by functions"
        req_dates::ClosedInterval{Date}
    end


This type provides [Tables.jl](https://github.com/JuliaData/Tables.jl) access to
`DataVector`s.

Functions to update some values are:
- `AbnormalReturns.update_id!`: changes the `id` (which updates `dates` and `missing_bdays`)
- `select!`: changes the `cols` (also updates `dates` and `missing_bdays` if `cols` is different)
- `AbnormalReturns.update_dates!`: changes `req_dates` (does not update anything else)
"""
struct FixedTable{N, T, AV <: AbstractVector{T}} <: Tables.AbstractColumns
    data::NTuple{N, AV}
    cols::SVector{N, TimelineColumn}
    function FixedTable(xs::NTuple{N, AV}, cols::SVector{N, TimelineColumn}) where {T, N, AV<:AbstractVector{T}}
        @assert all_equal_length(xs) "Not all the same length"
        new{N, T, AV}(xs, cols)
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
    valuecols_market=nothing,
    valuecols_firms=nothing
)
    df_market = DataFrame(df_market)
    df_firms = DataFrame(df_firms)
    if valuecols_market === nothing
        valuecols_market = Symbol.([n for n in Symbol.(names(df_market)) if n ∉ [date_col_market]])
    end
    if valuecols_firms === nothing
        valuecols_firms = Symbol.([n for n in Symbol.(names(df_firms)) if n ∉ [date_col_firms, id_col]])
    end

    df_market = select(df_market, vcat([date_col_market], valuecols_market))
    dropmissing!(df_market, date_col_market)
    #dropmissing!(df_market)
    sort!(df_market)

    if any(nonunique(df_market, [date_col_market]))
        @error("There are duplicate date rows in the market data")
    end

    df_firms = select(df_firms, vcat([id_col, date_col_firms], valuecols_firms))
    dropmissing!(df_firms, [id_col, date_col_firms])
    sort!(df_firms, [id_col, date_col_firms])

    # since this is sorted, just a single iteration is enough to check
    if all_unique_obs(df_firms[:, id_col], df_firms[:, date_col_firms])
        @error("There are duplicate id-date rows in the firm data")
    end

    cal = MarketCalendar(df_market[:, date_col_market])

    market_data = NamedTuple(
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


    firm_data = Dict{typeof(col_tab[id_col][1][1]), NamedTuple{Tuple(valuecols_firms), NTuple{length(valuecols_firms), DataVector}}}()
    sizehint!(firm_data, length(col_tab[id_col]))

    for i in 1:length(col_tab[id_col])
        firm_data[col_tab[id_col][i][1]] = NamedTuple(
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
        data = OffsetVector(coalesce.(data, zero(nonmissingtype(eltype(data)))), offset)
    else
        missing_days = nothing#OffsetVector(spzeros(Bool, size(data, 1)), offsets)
    end
    f(
        data,
        missing_days,
        offset:offset+length(data)-1
    )
end

function DataVector(data::AbstractVector{T}, offset::Int) where {T}
    if any(ismissing.(data))
        if all(ismissing.(data))
            return DataVector(
                OffsetVector(zeros(nonmissingtype(T), 1), offset),
                OffsetVector(sparsevec([1], [true]), offset),
                offset:offset
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

function get_missing_bdays(data::MarketData{T}, id::T, r::UnitRange{Int}, col::TimelineColumn) where {T}
    mssngs = data_missing_bdays(data[id, col])
    if mssngs === nothing
        nothing
    else
        mssngs[r .- shift_count(col)]
    end
end

combine_missing_bdays(vals::Nothing...) = nothing
function combine_missing_bdays(vals::Union{Nothing, SparseVector{Bool, Int}}...)::SparseVector{Bool, Int}
    (|).(vals...)
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

# function Base.getindex(
#     data::MarketData{T, MNames, FNames},
#     id::T,
#     cols::Vector{TimelineColumn}=TimelineColumn.([FNames..., MNames...]),
#     allow_mssng::Type{AllowMissing{Mssng}}=AllowMissing{false},
# ) where {T, MNames, FNames, Mssng}
#     real_dates = range_available(get_all_ranges(data, id, cols))
#     data[id, real_dates, cols, allow_mssng]
# end

Base.getindex(data::NamedTuple, col::TimelineColumn) = data[Symbol(col)]

function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    col::Symbol
) where {T, MNames, FNames}
    if col ∈ MNames
        data.marketdata[col]
    else
        data.firmdata[id][col]
    end
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    col::TimelineColumn
) where {T}
    data[id, Symbol(col)]
end

function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    r::UnitRange,
    col::TimelineColumn,
    missing_bdays::AbstractVector{Bool}
) where {T, MNames, FNames}
    if Symbol(col) ∈ MNames
        view(data[col], (r .+ col.shifts)[Not(missing_bdays)])
    else
        view(data[id][col], (r .+ col.shifts)[Not(missing_bdays)])
    end
end

function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    r::UnitRange,
    col::TimelineColumn,
) where {T, MNames, FNames}
    if Symbol(col) ∈ MNames
        view(data.marketdata[col], r .+ col.shifts)
    else
        view(data.firmdata[id][col], r .+ col.shifts)
    end
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange{Int},
    cols::SVector{N, TimelineColumn},
    ::Nothing
) where {T, N}
    FixedTable(
        NTuple{N}((data[id, r, col] for col in cols)),
        cols
    )
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    r::UnitRange{Int},
    cols::SVector{N, TimelineColumn},
    missing_days::AbstractVector{Bool}
) where {T, N}
    @assert length(missing_days) == length(r) "Missing days is wrong length"
    FixedTable(
        NTuple{N}((data[id, r, col, missing_days] for col in cols)),
        cols
    )
end


function Base.getindex(
    data::MarketData{T},
    id::T,
    dates::ClosedInterval{Date},
    cols::SVector{N, TimelineColumn},
    #missing_days::Union{Nothing, AbstractVector{Bool}}=nothing
) where {T, N}
    r = maximin(date_range(calendar(data), dates), interval.(Ref(data), id, cols)...)
    mssngs = combine_missing_bdays(get_missing_bdays.(Ref(data), id, Ref(r), cols)...)
    data[id, r, cols, mssngs]
end

function Base.getindex(
    data::MarketData{T},
    id::T,
    dates::ClosedInterval{Date},
    cols::AbstractVector,
    #missing_days::Union{Nothing, AbstractVector{Bool}}=nothing
) where {T}
    cols = SVector{length(cols)}(TimelineColumn.(cols))
    data[id, dates, cols]
end

function Base.getindex(data::FixedTable{N}, ::Colon, i::Int) where {N}
    @assert 1 <= i <= N "Column not in data"
    data.data[i]
end

function Base.getproperty(data::FixedTable, key::Symbol)
    getfield(data, key)
end

function Base.getindex(
    data::FixedTable,
    ::Colon,
    col::TimelineColumn
)
    @assert col ∈ names(data) "Column is not in the data"
    i = findfirst(==(col), names(data))
    data[:, i]
end


function Base.getindex(
    data::FixedTable,
    ::Colon,
    col
)
    data[:, TimelineColumn(col)]
end


Base.names(x::FixedTable) = x.cols

Tables.istable(::Type{<:FixedTable}) = true
Tables.columnaccess(::Type{<:FixedTable}) = true
function Tables.schema(obj::FixedTable{N, T}) where {N, T}
    Tables.Schema(Tables.columnnames(obj), fill(T, N))
end

Tables.columns(x::FixedTable) = x

Tables.getcolumn(x::FixedTable, nm::Symbol) = x[:, nm]
Tables.getcolumn(x::FixedTable, nm::Int) = x[:, nm]
Tables.getcolumn(x::FixedTable, nm::TimelineColumn) = x[:, nm]




function Tables.columnnames(x::FixedTable)
    Symbol.(String.(names(x)))
end

# function Base.length(data::TimelineTable{false})
#     real_dates = dates_min_max(data_dates(data), norm_dates(data))
#     c = bdayscount(calendar(data), dt_min(real_dates), dt_max(real_dates)) + isbday(calendar(data), dt_max(real_dates))
#     new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), real_dates)
#     return c - nnz(new_mssngs)
# end


function Base.length(x::DataVector)
    return length(data_missing_bdays(x)) - nnz(data_missing_bdays(x))
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

function interval(data::MarketData{T}, id::T, col::TimelineColumn) where {T}
    r = interval(data[id, col])
    if shift_count(col) == 0
        r
    else
        maximin(r, r .+ shift_count(col))
    end
end




calendar(x::MarketData) = x.calendar

firmdata(data::MarketData{T}, i::T, sym::Symbol) where {T} = data.firmdata[i][sym]
marketdata(data::MarketData, sym::Symbol) = data.marketdata[sym]
pick_data(data::MarketData{T, MNames, FNames}, id::T, sym::Symbol) where {T, MNames, FNames} = sym ∈ FNames ? firmdata(data, id, sym) : marketdata(data, sym)

pick_data(data::MarketData{T}, id::T, col::TimelineColumn) where {T} = pick_data(data, id, Symbol(col))

data_dates(x::MarketData) = x.dates

data_loc(data::MarketData{T, MNames, FNames}, sym::Symbol) where {T, MNames, FNames} = sym ∈ FNames ? :firmdata : :marketdata
function data_loc(data::MarketData{T, MNames, FNames}, syms::Vector{Symbol}) where {T, MNames, FNames}
    all_locs = data_loc.(Ref(data), syms)
    (
        syms[all_locs .== :firmdata],
        syms[all_locs .== :marketdata]
    )
end


function Base.show(io::IO, data::MarketData{T, MNames, FNames}) where {T, MNames, FNames}
    println(io, "MarketData with ID type $T with $(length(getfield(data, :firmdata))) unique firms")
    println(io, data.calendar)
    println(io, "Market Columns: $(join(MNames, ", "))")
    println(io, "Firm Columns: $(join(FNames, ", "))")
end

# function Base.show(io::IO, data::TimelineTable{Mssng, T, MNames, FNames}) where {Mssng, T, MNames, FNames}
#     println(io, "TimelineTable for $(data_id(data)) and available columns $(join([MNames..., FNames...], ","))")
#     println(io, "Maximum dates given current column selection: $(data_dates(data))")
#     println(io, "Currently selected dates: $(norm_dates(data))")
#     pretty_table(io, data; header = vcat(["Date"], names(data)), tf=tf_simple)
# end

dt_min(x::ClosedInterval{Date}) = x.left
dt_max(x::ClosedInterval{Date}) = x.right


cal_dt_min(x::MarketData) = cal_dt_min(calendar(x))
cal_dt_max(x::MarketData) = cal_dt_max(calendar(x))