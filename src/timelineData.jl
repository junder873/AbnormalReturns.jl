

struct DataVector
    data::Vector{Float64}
    missing_bdays::SparseVector{Bool, Int}
    dates::ClosedInterval{Date}
    calendar::MarketCalendar
    function DataVector(data, missing_bdays, dates, calendar)
        @assert length(data) == length(missing_bdays) == bdayscount(calendar, dates.left, dates.right) + 1 "Data does not match length of dates or missings"
        new(data, missing_bdays, dates, calendar)
    end
end

struct RegressionCache
    data::Matrix{Float64}
    dates::ClosedInterval{Date}
    calendar::MarketCalendar
end

struct MarketData{T, MNames, FNames, N1, N2}
    calendar::MarketCalendar
    marketdata::NamedTuple{MNames, NTuple{N1, DataVector}} # column names as symbols
    firmdata::Dict{T, NamedTuple{FNames, NTuple{N2, DataVector}}} # data stored by firm id and then by column name as symbol
end

mutable struct TimelineTable{Names, N}
    data::NamedTuple{Names, NTuple{N, DataVector}}
    dates::ClosedInterval{Date}# actual returnable dates that guarentees a square matrix
    cols::DictIndex
    missing_bdays::SparseVector{Bool, Int}
    req_dates::ClosedInterval{Date}# dates requested, matters for calculating missing obs
    calendar::MarketCalendar
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
        data = data[i:j]
        dates = dates[i:j]
        missing_days = sparsevec(findall(ismissing, data))
        data = coalesce.(data, zero(nonmissingtype(eltype(data))))
    else
        missing_days = sparsevec(zeros(Bool, length(data)))
    end
    
    
    DataVector(
        data,
        missing_days,
        dates[1] .. dates[end],
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

function dates_min_max(dates::Vector{ClosedInterval{Date}})::ClosedInterval{Date}
    maximum(dt_min.(dates)) .. minimum(dt_max.(dates))
end
function dates_min_max(date1::ClosedInterval{Date}, date2::ClosedInterval{Date})::ClosedInterval{Date}
    max(dt_min(date1), dt_min(date2)) .. min(dt_max(date1), dt_max(date2))
end

function all_dates(data::MarketData)
    shift_dates.(pick_data.(Ref(data), Symbol.(names(data))), shift_count.(names(data)))
end

function get_missing_bdays(cal::MarketCalendar, missing_bdays::SparseVector{Bool, Int}, dates_missings::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    missing_bdays[date_range(cal, dates_missings, new_dates)]
end

function get_missing_bdays(data::DataVector, new_dates::ClosedInterval{Date})
    get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), new_dates)
end

function date_range(cal::MarketCalendar, timeline_dates::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    dates = dates_min_max(timeline_dates, new_dates)
    s = bdayscount(cal, dt_min(timeline_dates), dt_min(dates)) + 1
    e = s + bdayscount(cal, dt_min(dates), dt_max(dates)) - !isbday(cal, dt_max(dates))
    s:e
end

function date_range(cal::MarketCalendar, data::DataVector, dt_min::Date, dt_max::Date)
    date_range(cal, data_dates(data), dt_min .. dt_max)
end

function date_range(cal::MarketCalendar, data::DataVector, dates::ClosedInterval{Date})
    date_range(cal, data_dates(data), dates)
end

function Base.getindex(data::RegressionCache, dates::ClosedInterval{Date}, mssngs::SparseVector{Bool, Int})
    new_dates = dates_min_max(data_dates(data), dates)
    new_mssngs = get_missing_bdays(calendar(data), mssngs, dates, new_dates)
    r = date_range(calendar(data), data_dates(data), new_dates)
    if nnz(new_mssngs) == 0
        return data.data[r, :]
    else
        data.data[r[(!).(new_mssngs)], :]
    end
end

function Base.getindex(data::DataVector, dates::ClosedInterval{Date})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    data.data[date_range(calendar(data), data_dates(data), dates)] 
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

function check_col(x::Symbol, g1, g2)
    x ∈ g1 || x ∈ g2
end

function check_col(x::Vector{Symbol}, g1, g2)
    all(check_col.(x, Ref(g1), Ref(g2)))
end

check_col(x::Symbol, g1) = x ∈ g1

function check_col(x::Vector{Symbol}, g1)
    any(check_col.(x, Ref(g1)))
end

function Base.view(data::DataVector, dates::ClosedInterval{Date})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    r = date_range(calendar(data), data_dates(data), dates)
    view(data.data, r)
end

function Base.view(data::RegressionCache, dates::ClosedInterval{Date})
    @assert dt_min(data) <= dt_min(dates) && dt_max(dates) <= dt_max(data) "Date range out of bounds"
    r = date_range(calendar(data), data_dates(data), dates)
    view(data.data, r, :)
end

###########################################################
# Access functions for multiple columns and range of dates
# returns another data matrix
###########################################################

Base.getindex(data::NamedTuple{Names, NTuple{N, DataVector}}, col::TimelineColumn) where {Names, N} = shift(data[Symbol(col)], col.shifts)

function Base.getindex(data::MarketData{T, MNames, FNames}, id::T, cols::Vector{TimelineColumn}=TimelineColumn.([FNames..., MNames...])) where {T, MNames, FNames}
    out_data = merge(
        data.firmdata[id][FNames],
        data.marketdata[MNames]
    )
    real_dates = dates_min_max([data_dates(out_data[col]) for col in cols])
    missing_bdays = combine_missing_bdays([get_missing_bdays(out_data[col], real_dates) for col in cols])
    TimelineTable(
        out_data,
        real_dates,
        DictIndex(cols),
        missing_bdays,
        real_dates,
        data.calendar
    )
end

function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    dates::ClosedInterval{Date},
    cols::Vector{TimelineColumn}=TimelineColumn.([FNames..., MNames...])
) where {T, MNames, FNames}
    out_data = merge(
        data.firmdata[id][FNames],
        data.marketdata[MNames]
    )
    real_dates = dates_min_max([data_dates(out_data[col]) for col in cols])
    missing_bdays = combine_missing_bdays([get_missing_bdays(out_data[col], real_dates) for col in cols])
    TimelineTable(
        out_data,
        dates,
        DictIndex(cols),
        missing_bdays,
        dates,
        data.calendar
    )
end

Base.getindex(data::TimelineTable, col::TimelineColumn) = col.shifts == 0 ? data.data[Symbol(col)] : shift(data.data[Symbol(col)], col.shifts)
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
    data::TimelineTable,
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

function Base.getindex(
    data::TimelineTable,
    dates::ClosedInterval{Date}
)
    TimelineTable(
        data.data,
        data.dates,
        DictIndex(names(data)),
        data.missing_bdays,
        dates,
        data.calendar
    )
end

function get_dates(data::TimelineTable)
    out = listbdays(data.calendar, dt_min(data), dt_max(data))
    if nnz(data.missing_bdays) == 0
        return out
    else
        return out[(!).(data.missing_bdays)]
    end
end

Base.names(x::TimelineTable) = x.cols.cols

Tables.istable(::Type{<:TimelineTable}) = true
Tables.columnaccess(::Type{<:TimelineTable}) = true
function Tables.schema(obj::TimelineTable) where {mssngs}
    col_sym = Tables.columnnames(obj)
    col_types = fill(Float64, length(col_sym))
    col_types[1] = Date
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
    if nm ∈ Tables.columnanmes(x)
        pos = findfirst(nm .== Tables.columnnames(x))
        return Tables.getcolumn(x, pos)
    end
end

function Tables.columnnames(x::MarketData)
    vcat([:date], Symbol.(String.(names(x))))
end


function Base.length(data::TimelineTable)
    real_dates = dates_min_max(data_dates(data), norm_dates(data))
    c = bdayscount(calendar(data), dt_min(real_dates), dt_max(real_dates)) + 1
    new_mssngs = get_missing_bdays(calendar(data), data_missing_bdays(data), data_dates(data), real_dates)
    return c - nnz(new_mssngs)
end

function Base.length(x::DataVector)
    return length(data_missing_bdays(x)) - nnz(data_missing_bdays(x))
end

function update_dates!(
    data::TimelineTable{Names, N},
    dates::ClosedInterval{Date}
) where {Names, N}
    data.req_dates = dates
    data
end



function DataFrames.select!(x::TimelineTable, cols::Vector)
    select!(x, TimelineColumn.(cols))
end
function DataFrames.select!(data::TimelineTable, cols::Vector{TimelineColumn})
    updates = sort(cols) != sort(names(data)) # if the columns are changing, then update dates and missing_bdays
    data.cols = DictIndex(cols)
    if updates
        data.dates = dates_min_max([data_dates(data[col]) for col in names(data)])
        data.missing_bdays = combine_missing_bdays([get_missing_bdays(data[col], data_dates(data)) for col in names(data)])
    end
    data
end


raw_values(x::DataVector) = x.data
data_dates(x::DataVector) = x.dates
data_missing_bdays(x::DataVector) = x.missing_bdays

raw_values(x::RegressionCache) = x.data
data_dates(x::RegressionCache) = x.dates
data_missing_bdays(x::RegressionCache) = x.missing_bdays

data_dates(x::TimelineTable) = x.dates
data_missing_bdays(x::TimelineTable) = x.missing_bdays
norm_dates(x::TimelineTable) = x.req_dates

calendar(x::DataVector) = x.calendar
calendar(x::RegressionCache) = x.calendar
calendar(x::MarketData) = x.calendar
calendar(x::TimelineTable) = x.calendar

firmdata(x::MarketData{T}, i::T, sym::Symbol) where {T} = x.firmdata[i][sym]
marketdata(x::MarketData, sym::Symbol) = x.marketdata[sym]
pick_data(data::MarketData{T, MNames, FNames}, id::T, sym::Symbol) where {T, MNames, FNames} = sym ∈ FNames ? firmdata(data, id, sym) : marketdata(data, sym)

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
dt_min(x::TimelineTable) = x.req_dates.left
dt_min(x::DataVector) = x.dates.left
dt_min(x::RegressionCache) = x.dates.left
dt_max(x::ClosedInterval{Date}) = x.right
dt_max(x::TimelineTable) = x.req_dates.right
dt_max(x::DataVector) = x.dates.right
dt_max(x::RegressionCache) = x.dates.right

cal_dt_min(x::DataVector) = cal_dt_min(calendar(x))
cal_dt_max(x::DataVector) = cal_dt_max(calendar(x))
cal_dt_min(x::TimelineTable) = cal_dt_min(calendar(x))
cal_dt_max(x::TimelineTable) = cal_dt_max(calendar(x))
cal_dt_min(x::MarketData) = cal_dt_min(calendar(x))
cal_dt_max(x::MarketData) = cal_dt_max(calendar(x))