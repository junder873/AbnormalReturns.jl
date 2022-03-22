
abstract type CalendarData end

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

struct DataMatrix <: CalendarData
    data::Matrix{Float64}
    missing_bdays::SparseVector{Bool, Int}# corresponds to each row with a missing value
    dates::ClosedInterval{Date}
    calendar::MarketCalendar
    function DataMatrix(data, missing_bdays, dates, calendar)
        @assert size(data, 1) == length(missing_bdays) == bdayscount(calendar, dates.left, dates.right) + 1 "Data does not match length of dates or missings"
        new(data, missing_bdays, dates, calendar)
    end
end

struct MarketData{T, MNames, FNames, N1, N2}
    calendar::MarketCalendar
    marketdata::NamedTuple{MNames, NTuple{N1, DataVector}} # column names as symbols
    firmdata::Dict{T, NamedTuple{FNames, NTuple{N2, DataVector}}} # data stored by firm id and then by column name as symbol
end

struct AllowMissing{mssng} end

mutable struct TimelineTable{Mssng, T, MNames, FNames, N1, N2} <: Tables.AbstractColumns
    parent::MarketData{T, MNames, FNames, N1, N2}
    id::T
    allow_missing::Type{AllowMissing{Mssng}}
    dates::ClosedInterval{Date}# actual returnable dates that guarentees a square matrix
    cols::DictIndex
    missing_bdays::SparseVector{Bool, Int}
    req_dates::ClosedInterval{Date}# dates requested, matters for calculating missing obs
end

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

    if any(nonunique(df_firms, [id_col, date_col_firms]))
        @error("There are duplicate id-date rows in the firm data")
    end

    cal = MarketCalendar(df_market[:, date_col_market])

    market_data = NamedTuple(
        valuecols_market .=>
        DataVector.(
            Tables.columns(df_market[:, valuecols_market]),
            Ref(df_market[:, date_col_market]),
            Ref(cal)
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
                    col_tab[date_col_firms][i],
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

function setup_calendar_data(f, data, dates::ClosedInterval{Date}, cal::MarketCalendar)
    if any(ismissing.(data))
        missing_days = sparsevec(any(ismissing, data, dims=2))
        data = coalesce.(data, zero(nonmissingtype(eltype(data))))
    else
        missing_days = spzeros(Bool, size(data, 1))
    end
    f(
        data,
        missing_days,
        dates,
        cal
    )
end

function DataVector(data::AbstractVector, dates::ClosedInterval{Date}, cal::MarketCalendar)
    setup_calendar_data(
        DataVector,
        data,
        dates,
        cal
    )
end

function DataVector(data::AbstractVector, dates::AbstractVector{Date}, cal::MarketCalendar)
    if any(ismissing.(data))
        if all(ismissing.(data))
            return DataVector(
                zeros(nonmissingtype(eltype(data)), 1),
                SparseVector([1], [true]),
                dates[1] .. dates[1],
                cal
            )
        end
        i = findfirst(!ismissing, data)
        j = findlast(!ismissing, data)
        setup_calendar_data(
            DataVector,
            data[i:j],
            dates[i] .. dates[j],
            cal
        )
    else
        setup_calendar_data(
            DataVector,
            data,
            dates[1] .. dates[end],
            cal
        )
    end
end

function DataMatrix(data::AbstractMatrix, dates::ClosedInterval{Date}, cal::MarketCalendar)
    setup_calendar_data(
        DataMatrix,
        data,
        dates,
        cal
    )
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

function get_missing_bdays(cal::MarketCalendar, missing_bdays::SparseVector{Bool, Int}, dates_missings::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    missing_bdays[date_range(cal, dates_missings, new_dates)]
end

function get_missing_bdays(data, new_dates::ClosedInterval{Date})
    get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), new_dates)
end

function date_range(cal::MarketCalendar, timeline_dates::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    dates = dates_min_max(timeline_dates, new_dates)
    s = bdayscount(cal, dt_min(timeline_dates), dt_min(dates)) + 1
    e = s + bdayscount(cal, dt_min(dates), dt_max(dates)) - !isbday(cal, dt_max(dates))
    s:e
end

function date_range(cal::MarketCalendar, data::CalendarData, dt_min::Date, dt_max::Date)
    date_range(cal, data_dates(data), dt_min .. dt_max)
end

function date_range(cal::MarketCalendar, data::CalendarData, dates::ClosedInterval{Date})
    date_range(cal, data_dates(data), dates)
end

# , mssngs::SparseVector{Bool, Int}=spzeros(Bool, size(raw_values(data), 1))
function Base.getindex(data::DataMatrix, dates::ClosedInterval{Date})
    new_dates = dates_min_max(data_dates(data), dates)
    r = date_range(calendar(data), data_dates(data), new_dates)
    data.data[r, :]
end

function Base.getindex(data::DataVector, dates::ClosedInterval{Date})
    new_dates = dates_min_max(data_dates(data), dates)
    r = date_range(calendar(data), data_dates(data), new_dates)
    data.data[r]
end

function combine_missing_bdays(vals::SparseVector{Bool, Int}...)::SparseVector{Bool, Int}
    (|).(vals...)
end

function combine_missing_bdays(vals::Vector{SparseVector{Bool, Int}})::SparseVector{Bool, Int}
    out = sparsevec(zeros(Bool, length(vals[1])))
    for v in vals
        out = out .| v
    end
    out
end

function Base.view(data::DataVector, dates::ClosedInterval{Date}, mssngs::SparseVector{Bool, Int})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    r = date_range(calendar(data), data_dates(data), dates)
    @assert length(mssngs) == length(r) "Missing data vector is the wrong length"
    if nnz(mssngs) == 0
        view(data.data, r)
    else
        view(data.data, r[(!).(mssngs)])
    end
end

function Base.view(data::DataMatrix, dates::ClosedInterval{Date}, mssngs::SparseVector{Bool, Int})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    r = date_range(calendar(data), data_dates(data), dates)
    @assert length(mssngs) == length(r) "Missing data vector is the wrong length"
    if nnz(mssngs) == 0
        view(data.data, r, :)
    else
        view(data.data, r[(!).(mssngs)], :)
    end
end

function Base.view(data::DataVector, dates::ClosedInterval{Date})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    r = date_range(calendar(data), data_dates(data), dates)
    view(data.data, r)
end

function Base.view(data::DataMatrix, dates::ClosedInterval{Date})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    r = date_range(calendar(data), data_dates(data), dates)
    view(data.data, r, :)
end

Base.getindex(data::NamedTuple{Names, NTuple{N, DataVector}}, col::TimelineColumn) where {Names, N} = shift(data[Symbol(col)], col.shifts)

function Base.getindex(data::MarketData{T, MNames, FNames}, id::T, col::TimelineColumn) where {T, MNames, FNames}
    shift(pick_data(data, id, Symbol(col)), col.shifts)
end

function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    cols::Vector{TimelineColumn}=TimelineColumn.([FNames..., MNames...]),
    allow_mssng::Type{AllowMissing{Mssng}}=AllowMissing{false},
) where {T, MNames, FNames, Mssng}
    real_dates = dates_min_max(get_all_dates(data, id, cols))
    data[id, real_dates, cols, allow_mssng]
end

function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    dates::ClosedInterval{Date},
    cols::Vector{TimelineColumn}=TimelineColumn.([FNames..., MNames...]),
    allow_mssng::Type{AllowMissing{Mssng}}=AllowMissing{false}
) where {T, MNames, FNames, Mssng}
    out = TimelineTable(
        data,
        id,
        allow_mssng,
        dates_min_max(get_all_dates(data, id, cols)),
        DictIndex(cols),
        sparsevec(zeros(Bool, bdayscount(calendar(data), dt_min(dates), dt_max(dates)))),
        dates,
    )
    out.missing_bdays = combine_all_missing_bdays(out)
    out
end

function Base.getindex(data::TimelineTable, col::TimelineColumn)
    parent(data)[data.id, col]
end
Base.getindex(data::TimelineTable, col::Symbol) = data[TimelineColumn(col)]

function Base.getindex(
    data::TimelineTable,
    ::Colon,
    col::TimelineColumn
)
    ret_dates = dates_min_max(data_dates(data), norm_dates(data))
    data[ret_dates, col]
end

function Base.getindex(
    data::TimelineTable{false},
    dates::ClosedInterval{Date},
    col::TimelineColumn
)
    new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), dates)
    if nnz(new_mssngs) == 0
        data[col][dates]
    else
        data[col][dates][(!).(new_mssngs)]
    end
end

function Base.getindex(
    data::TimelineTable{true},
    dates::ClosedInterval{Date},
    col::TimelineColumn
)
    new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), dates)
    out = data[col][dates]
    out = allowmissing(out)
    out[new_mssngs] .= missing
    out
end

function Base.getindex(
    data::TimelineTable,
    ::Colon,
    col::Symbol
)
    data[:, TimelineColumn(col)]
end

function Base.getindex(
    data::TimelineTable,
    dates::ClosedInterval{Date},
    col::Symbol
)
    data[dates, TimelineColumn(col)]
end

function Base.getproperty(
    data::TimelineTable,
    sym::Symbol
)
    if sym == :calendar
        getfield(getfield(data, :parent), :calendar)
    elseif sym == :firmdata
        getfield(getfield(data, :parent), :firmdata)
    elseif sym == :marketdata
        getfield(getfield(data, :parent), :marketdata)
    else
        getfield(data, sym)
    end
end



function get_dates(data::TimelineTable{true})
    listbdays(calendar(data), dt_min(data), dt_max(data))
end

function get_dates(data::TimelineTable{false})
    real_dates = dates_min_max(data_dates(data), norm_dates(data))
    out = listbdays(calendar(data), dt_min(real_dates), dt_max(real_dates))
    new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), real_dates)
    if nnz(new_mssngs) == 0
        return out
    else
        return out[(!).(new_mssngs)]
    end
end

Base.names(x::TimelineTable) = x.cols.cols

Tables.istable(::Type{<:TimelineTable}) = true
Tables.columnaccess(::Type{<:TimelineTable}) = true
function Tables.schema(obj::TimelineTable{Mssngs}) where {Mssngs}
    col_sym = Tables.columnnames(obj)
    col_types = vcat(
        [Date],
        fill(Mssngs ? Union{Float64, Missing} : Float64, length(col_sym)-1)
    )
    Tables.Schema(col_sym, col_types)
end

Tables.columns(x::TimelineTable) = x
function Tables.getcolumn(x::TimelineTable, i::Int)
    if i == 1
        get_dates(x)
    else
        Tables.getcolumn(x, x.lookup[i-1])# subtract 1 since date is the first column produced
    end
end
function Tables.getcolumn(x::TimelineTable, nm::TimelineColumn)
    x[:, nm]
end
function Tables.getcolumn(x::TimelineTable, nm::Symbol)
    if nm == :date
        return get_dates(x)
    end
    if nm ∈ Tables.columnnames(x)
        pos = findfirst(nm .== Tables.columnnames(x))
        return Tables.getcolumn(x, pos)
    end
end

function Tables.columnnames(x::TimelineTable)
    vcat([:date], Symbol.(String.(names(x))))
end

function Base.length(data::TimelineTable{false})
    real_dates = dates_min_max(data_dates(data), norm_dates(data))
    c = bdayscount(calendar(data), dt_min(real_dates), dt_max(real_dates)) + 1
    new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), real_dates)
    return c - nnz(new_mssngs)
end

function Base.length(data::TimelineTable{true})
    bdayscount(calendar(data), dt_min(data), dt_max(data)) + 1
end

function Base.length(x::DataVector)
    return length(data_missing_bdays(x)) - nnz(data_missing_bdays(x))
end

DataFrames.dropmissing(data::TimelineTable{false}) = data
DataFrames.allowmissing(data::TimelineTable{true}) = data

function DataFrames.dropmissing(data::TimelineTable{true})
    TimelineTable(
        parent(data),
        data_id(data),
        AllowMissing{false},
        data_dates(data),
        DictIndex(names(data)),
        data_missing_bdays(data),
        norm_dates(data)
    )
end

function DataFrames.allowmissing(data::TimelineTable{false})
    TimelineTable(
        parent(data),
        data_id(data),
        AllowMissing{true},
        data_dates(data),
        DictIndex(names(data)),
        data_missing_bdays(data),
        norm_dates(data)
    )
end

function update_dates!(
    data::TimelineTable,
    dates::ClosedInterval{Date}
)
    data.req_dates = dates
    data
end

function update_id!(
    data::TimelineTable{Mssng, T},
    id::T
) where {Mssng, T}
    if id != data.id
        data.id = id
        data.dates = dates_min_max(get_all_dates(data))
        data.missing_bdays = combine_all_missing_bdays(data)
    end
    data
end

function DataFrames.select!(x::TimelineTable, cols::Vector)
    select!(x, TimelineColumn.(cols))
end
function DataFrames.select!(data::TimelineTable, cols::Vector{TimelineColumn})
    updates = sort(cols) != sort(names(data)) # if the columns are changing, then update dates and missing_bdays
    data.cols = DictIndex(cols)
    if updates
        data.dates = dates_min_max(get_all_dates(data))
        data.missing_bdays = combine_all_missing_bdays(data)
    end
    data
end

function get_all_dates(data::MarketData{T}, id::T, cols::Vector{TimelineColumn}) where {T}
    shift_dates.(
        pick_data.(Ref(data), id, Symbol.(cols)),
        shift_count.(cols)
    )
end

function get_all_dates(data::TimelineTable)
    get_all_dates(parent(data), data.id, names(data))
end

function combine_all_missing_bdays(data::TimelineTable)
    out = sparsevec(
        zeros(
            Bool,
            bdayscount(
                calendar(data),
                dt_min(data_dates(data)),
                dt_max(data_dates(data))
            ) + 1
        )
    )
    for col in names(data)
        mssngs = data_missing_bdays(data[col])
        if nnz(mssngs) > 0
            out = out .| get_missing_bdays(data[col], data_dates(data))
        end
    end
    out
end

raw_values(x::CalendarData) = x.data
data_dates(x::CalendarData) = x.dates
data_missing_bdays(x::CalendarData) = x.missing_bdays

data_dates(x::TimelineTable) = x.dates
data_missing_bdays(x::TimelineTable) = x.missing_bdays
norm_dates(x::TimelineTable) = x.req_dates
parent(data::TimelineTable) = data.parent

calendar(x::CalendarData) = x.calendar
calendar(x::MarketData) = x.calendar
calendar(x::TimelineTable) = calendar(parent(x))

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

dt_min(x::ClosedInterval{Date}) = x.left
dt_min(x::TimelineTable) = dt_min(norm_dates(x))
dt_min(x::CalendarData) = dt_min(data_dates(x))
dt_max(x::ClosedInterval{Date}) = x.right
dt_max(x::TimelineTable) = dt_max(norm_dates(x))
dt_max(x::CalendarData) = dt_max(data_dates(x))

cal_dt_min(x::CalendarData) = cal_dt_min(calendar(x))
cal_dt_max(x::CalendarData) = cal_dt_max(calendar(x))
cal_dt_min(x::TimelineTable) = cal_dt_min(calendar(x))
cal_dt_max(x::TimelineTable) = cal_dt_max(calendar(x))
cal_dt_min(x::MarketData) = cal_dt_min(calendar(x))
cal_dt_max(x::MarketData) = cal_dt_max(calendar(x))