

struct DataVector
    data::Vector{Float64}
    missing_bdays::Union{Nothing, Set{Int}}
    dates::ClosedInterval{Date}
end

mutable struct RegressionCache
    terms::MatrixTerm
    data::Matrix{Float64}
    dates::ClosedInterval{Date}
    missing_bdays::Union{Nothing, Vector{Int}}
    calendar::MarketCalendar
end

mutable struct MarketData
    calendar::MarketCalendar
    colmapping::Dict{Symbol, Symbol} # data is stored either in marketdata or in firmdata
    marketdata::Dict{Symbol, DataVector} # column names as symbols
    firmdata::Dict{Int, Dict{Symbol, DataVector}} # data stored by firm (int) and then by column name as symbol
    regression_cache::Union{Nothing, RegressionCache}
end

abstract type TimelineTable <: Tables.AbstractColumns end

mutable struct TimelineTableNoMissing <: TimelineTable
    parent::MarketData
    firmdata::Union{Nothing, Dict{Symbol, DataVector}}
    dates::ClosedInterval{Date}
    cols::DictIndex
    missing_bdays::Union{Nothing, Vector{Int}}# nothing if no missing bdays, empty vector if never calculated
end

mutable struct TimelineTableWithMissing <: TimelineTable
    parent::MarketData
    firmdata::Union{Nothing, Dict{Symbol, DataVector}}
    dates::ClosedInterval{Date}
    cols::DictIndex
    missing_bdays::Union{Nothing, Vector{Int}}# nothing if no missing bdays, empty vector if never calculated
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
    dropmissing!(df_market)
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

    market_data = Dict{Symbol, DataVector}()
    dates = minimum(df_market[:, date_col_market]) .. maximum(df_market[:, date_col_market])
    for col in valuecols_market
        market_data[col] = DataVector(
            df_market[:, col],
            nothing,
            dates
        )
    end

    col_map = Dict{Symbol, Symbol}()
    for col in Symbol.(valuecols_market)
        col_map[col] = :marketdata
    end
    for col in Symbol.(valuecols_firms)
        col_map[col] = :firmdata
    end

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

    firm_matrix = Matrix(df_firms[:, valuecols_firms])
    firm_dates = df_firms[:, date_col_firms]
    firm_ids = unique(df_firms[:, id_col])

    firm_data = Dict{Int, Dict{Symbol, DataVector}}()
    sizehint!(firm_data, length(firm_ids))

    for firm_id in firm_ids
        idx = gdf.keymap[(firm_id,)]
        s = gdf.starts[idx]
        e = gdf.ends[idx]
        firm_data[firm_id] = Dict{Symbol, DataVector}()
        for (j, col) in enumerate(valuecols_firms)
            firm_data[firm_id][col] = DataVector(
                firm_matrix[s:e, j],
                firm_dates[s:e]
            )
        end
    end

    MarketData(
        cal,
        col_map,
        market_data,
        firm_data,
        nothing
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

function DataVector(data::AbstractVector, dates::Vector{Date})
    if any(ismissing.(data))
        if all(ismissing.(data))
            return DataVector(
                zeros(nonmissingtype(eltype(data)), 1),
                Set([1]),
                dates[1] .. dates[1]
            )
        end
        i = findfirst(!ismissing, data)
        j = findlast(!ismissing, data)
        data = data[i:j]
        dates = dates[i:j]
        missing_days = findall(ismissing, data)
        data = coalesce.(data, zero(nonmissingtype(eltype(data))))
    else
        missing_days = Int[]
    end
    
    
    DataVector(
        data,
        length(missing_days) == 0 ? nothing : Set(missing_days),
        dates[1] .. dates[end]
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

dates_min_max(dates::Vector{ClosedInterval{Date}}) = dates_min_max(dates...)
function dates_min_max(dates::ClosedInterval{Date}...)
    maximum([d.left for d in dates]) .. minimum([d.right for d in dates])
end

function adjust_missing_bdays(cal::MarketCalendar, missing_bdays::Vector{Int}, dates_missings::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    if new_dates.right < dates_missings.left || new_dates.left > dates_missings.right
        return nothing
    end
    #if new_dates.left <= dates_missings.left
    missing_bdays = missing_bdays .+ (bdayscount(cal, new_dates.left, dates_missings.left))
    e = bdayscount(cal, new_dates.left, new_dates.right) + 1
    return missing_bdays[1 .<= missing_bdays .<= e]
end

function adjust_missing_bdays(cal::MarketCalendar, missing_bdays::Set{Int}, dates_missings::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    if new_dates.right < dates_missings.left || new_dates.left > dates_missings.right
        return nothing
    end
    return adjust_missing_bdays(cal, missing_bdays |> collect, dates_missings, new_dates)
end

adjust_missing_bdays(cal::MarketCalendar, missing_bdays::Nothing, dates_missings::ClosedInterval{Date}, new_dates::ClosedInterval{Date}) = nothing

function date_range(cal::MarketCalendar, timeline_dates::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    dates = dates_min_max(timeline_dates, new_dates)
    s = bdayscount(cal, timeline_dates.left, dates.left) + 1
    e = s + bdayscount(cal, dates.left, dates.right) - !isbday(cal, dates.right)
    s:e
end

function date_range(cal::MarketCalendar, data::DataVector, dt_min::Date, dt_max::Date)
    data_range(cal, data.dates, dt_min .. dt_max)
end

# function Base.getindex(data::RegressionData, dates::ClosedInterval{Date})
#     new_dates = dates_min_max(data.dates, dates)
#     new_missings = adjust_missing_bdays(data, new_dates)
#     new_missings2 = new_missings === nothing ? nothing : Set(new_missings)
#     DataVector(raw_values(data)[date_range(data, new_dates)], new_missings2, new_dates, data.calendar)
# end

function Base.getindex(data::RegressionCache, dates::ClosedInterval{Date})
    new_dates = dates_min_max(data.dates, dates)
    r = date_range(data.calendar, data.dates, new_dates)
    new_missings = adjust_missing_bdays(data.calendar, data.dates, new_dates)
    if new_missings === nothing
        data.data[r, :]
    else
        data.data[r[Not(new_missings)], :]
    end
end

function combine_missing_bdays(vals::Union{Nothing, Set{Int}}...)
    out = Set{Int}[]
    for x in vals
        x === nothing && continue
        out = union(out, x)
    end
    if length(out) == 0
        return nothing
    else
        return out
    end
end



###########################################################
# Access functions for multiple columns and range of dates
# returns another data matrix
###########################################################

function Base.getindex(data::MarketData, col::TimelineColumn)
    if !haskey(data.colmapping, Symbol(col))
        throw("Column is not in the data")
    end
    if data.colmapping[Symbol(col)] == :firmdata
        throw("An ID must be supplied to access columns of firm data")
    end

    out = data.marketdata[Symbol(col)]
    if col.shifts == 0
        return out
    end
    return shift(out, col.shifts, data.calendar)
end

function Base.getindex(data::MarketData, id::Int, col::TimelineColumn)
    if !haskey(data.colmapping, Symbol(col))
        throw("Column is not in the data")
    end

    if data.colmapping[Symbol(col)] == :marketdata
        return data[col]
    end
    out = data.firmdata[id][Symbol(col)]
    if col.shifts == 0
        return out
    end
    return shift(out, col.shifts, data.calendar)
end



function Base.getindex(
    data::MarketData,
    id::Int,
    dates::ClosedInterval{Date},
    cols::Vector
    )
    data[id, dates, TimelineColumn.(cols)]
end
function Base.getindex(
    data::MarketData,
    id::Int,
    dates::ClosedInterval{Date},
    cols::Vector{TimelineColumn}
    )
    if !haskey(data.firmdata, id)
        throw("Data for firm id $id is not stored in the data")
    end
    if any([!haskey(data.colmapping, x) for x in Symbol.(cols) if x != :date])
        throw("Not all columns are in the data")
    end

    return TimelineTableNoMissing(
        data,
        data.firmdata[id],
        dates,
        DictIndex(cols),
        Int[]
    )
end

function Base.getindex(data::MarketData, id::Int, ::Colon, cols::Vector)
    data[id, :, TimelineColumn.(cols)]
end
function Base.getindex(data::MarketData, id::Int, ::Colon, cols::Vector{TimelineColumn})
    if !haskey(data.firmdata, id)
        throw("Data for firm id $id is not stored in the data")
    end
    if any([!haskey(data.colmapping, x) for x in Symbol.(cols) if x != :date])
        throw("Not all columns are in the data")
    end

    firmdata = data.firmdata[id]
    dates = dates_min_max(data_dates.([data[id, col] for col in cols if Symbol(col) != :date]))
    TimelineTableNoMissing(
        data,
        firmdata,
        dates,
        DictIndex(cols),
        Int[]
    )
end

function Base.getindex(data::MarketData, id::Int, dates::ClosedInterval{Date}, ::Colon)
    data[id, dates, TimelineColumn.(vcat([:date], Symbol.(keys(data.colmapping))))]
end

function Base.getindex(data::MarketData, id::Int, ::Colon, ::Colon)
    data[id, :, TimelineColumn.(vcat([:date], Symbol.(keys(data.colmapping))))]
end

function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, cols::Vector)
    data[dates, TimelineColumn.(cols)]
end
function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, cols::Vector{TimelineColumn})
    if any([!haskey(data.colmapping, x) for x in Symbol.(cols) if x != :date])
        throw("Not all columns are in the data")
    end
    if any([data.colmapping[x] for x in Symbol.(cols) if x != :date] .== :firmdata)
        throw("An ID must be supplied to access columns of firm data")
    end

    TimelineTableNoMissing(
        data,
        nothing,
        dates,
        DictIndex(cols),
        Int[]
    )
end

Base.getindex(data::MarketData, ::Colon, cols::Vector) = data[:, TimelineColumn.(cols)]
function Base.getindex(data::MarketData, ::Colon, cols::Vector{TimelineColumn})
    dates = dates_min_max(data_dates.([data[col] for col in cols if Symbol(col) != :date])...)
    data[dates, cols]
end

Base.getindex(data::MarketData, id::Int, ::Colon, cols::Vector) = data[:, TimelineColumn.(cols)]
function Base.getindex(data::MarketData, id::Int, ::Colon, cols::Vector{TimelineColumn})
    dates = dates_min_max(data_dates.([data[id, col] for col in cols if Symbol(col) != :date])...)
    data[dates, cols]
end

function Base.getindex(data::MarketData, ::Colon, ::Colon)
    cols = vcat([:date], [col for col in Symbol.(keys(data.colmapping)) if data.colmapping[col] == :marketdata])
    data[:, cols]
end

function Base.getindex(data::MarketData, id::Int, ::Colon, ::Colon)
    cols = vcat([:date], [col for col in Symbol.(keys(data.colmapping))])
    data[id, :, cols]
end

function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, cols::Vector)
    data[dates, TimelineColumn.(cols)]
end
function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, cols::Vector{TimelineColumn})
    if !(data.dates == dates && sort(cols) == sort(names(data)))
        data.missing_bdays = Int[]
    end
    data.dates = dates
    data.cols = DictIndex(cols)
    return data
end

function Base.getindex(data::TimelineTable, ::Colon, cols::Vector)
    data[:, TimelineColumn.(cols)]
end
function Base.getindex(data::TimelineTable, ::Colon, cols::Vector{TimelineColumn})
    data[data.dates, cols]
end

function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, ::Colon)
    data[dates, names(data)]
end

function Base.getproperty(obj::TimelineTable, sym::Symbol)
    if sym == :calendar
        obj.parent.calendar
    elseif sym == :marketdata
        obj.parent.marketdata
    elseif sym == :colmapping
        obj.parent.colmapping
    elseif sym == :regression_cache
        obj.parent.regression_cache
    else
        getfield(obj, sym)
    end
end

function Base.getindex(obj::TimelineTable, nm::Symbol)
    obj[TimelineColumn(nm)]
end

function Base.getindex(obj::TimelineTable, nm::TimelineColumn)
    out = getproperty(obj, obj.colmapping[nm.name])[nm.name]
    if nm.shifts == 0
        return out
    end
    return shift(out, nm.shifts, obj.calendar)
end

function Base.getindex(obj::TimelineTableNoMissing, ::Colon, nm::Symbol)::Union{Vector{Date}, Vector{Float64}}
    obj[:, TimelineColumn(nm)]
end
function Base.getindex(obj::TimelineTableNoMissing, ::Colon, nm::TimelineColumn)::Union{Vector{Date}, Vector{Float64}}
    dates = dates_min_max(obj.dates, [obj[col].dates for col in names(obj) if Symbol(col) != :date]...)

    if Symbol(nm) == :date
        return listbdays(obj.calendar, dates.left, dates.right)
    end

    vec = obj[nm]
    local_dates = dates_min_max(dates, vec.dates)
    s = bdayscount(obj.calendar, vec.dates.left, local_dates.left) + 1
    e = s + bdayscount(obj.calendar, local_dates.left, local_dates.right) - !isbday(obj.calendar, local_dates.right)
    out = raw_values(vec)[s:e]
    if obj.missing_bdays !== nothing && length(obj.missing_bdays) == 0
        calculate_missing_bdays!(obj)
    end
    if obj.missing_bdays === nothing
        return out
    else
        return out[Not(adjust_missing_bdays(obj.calendar, obj.missing_bdays, obj.dates, dates))]
    end
end

function Base.getindex(obj::TimelineTableWithMissing, ::Colon, nm::Symbol)::Union{Vector{Date}, Vector{Union{Missing, Float64}}}
    obj[:, TimelineColumn(nm)]
end
function Base.getindex(obj::TimelineTableWithMissing, ::Colon, nm::TimelineColumn)::Union{Vector{Date}, Vector{Union{Missing, Float64}}}

    if Symbol(nm) == :date
        return listbdays(obj.calendar, obj.dates.left, obj.dates.right)
    end
    vec = obj[nm]
    local_dates = dates_min_max(obj.dates, vec.dates)
    s = bdayscount(obj.calendar, vec.dates.left, local_dates.left) + 1
    e = s + bdayscount(obj.calendar, local_dates.left, local_dates.right) - !isbday(obj.calendar, local_dates.right)
    out = raw_values(vec)[s:e]
    if obj.missing_bdays !== nothing && length(obj.missing_bdays) == 0
        calculate_missing_bdays!(obj)
    end
    to_add_pre = ifelse(
        obj.dates.left < local_dates.left,
        bdayscount(obj.calendar, obj.dates.left, local_dates.left),
        0
    )
    to_add_post = ifelse(
        obj.dates.right > local_dates.right,
        bdayscount(obj.calendar, local_dates.right, obj.dates.right),
        0
    )
    if vec.missing_bdays !== nothing
        out = allowmissing(out)
        out[adjust_missing_bdays(obj.calendar, vec.missing_bdays, vec.dates, local_dates)] .= missing
    end
    out = vcat(
        repeat([missing], to_add_pre),
        out,
        repeat([missing], to_add_post)
    )
    return out
end

function calculate_missing_bdays!(obj::TimelineTable)
    col_dates = dates_min_max(obj.dates, [obj[col].dates for col in names(obj) if Symbol(col) != :date]...)
    out = Set{Int}()
    if obj.dates.left < col_dates.left # set 1:first available dates as missing
        union!(out, Set(1:bdayscount(obj.calendar, obj.dates.left, col_dates.left)))
    end
    if obj.dates.right > col_dates.right # set last available date:end as missing
        r = Set(
            1:bdayscount(obj.calendar, col_dates.right, obj.dates.right) .+
            bdayscount(obj.calendar, obj.dates.left, col_dates.right)
        )
        union!(out, r)
    end
    for col in names(obj)
        Symbol(col) == :date && continue
        obj[col].missing_bdays === nothing && continue
        temp = adjust_missing_bdays(
            obj.calendar,
            obj[col].missing_bdays,
            obj[col].dates,
            obj.dates
        )
        if temp !== nothing
            union!(out, temp)
        end
    end
    if length(out) == 0
        obj.missing_bdays = nothing
    else
        obj.missing_bdays = collect(out)
    end
end

Base.names(x::TimelineTable) = x.cols.cols

Tables.istable(::Type{<:TimelineTable}) = true
Tables.columnaccess(::Type{<:TimelineTable}) = true
function Tables.schema(m::TimelineTableNoMissing)
    col_types = Type[]
    for col in names(m)
        if Symbol(col) == :date
            push!(col_types, Date)
        else
            push!(col_types, Float64)
        end
    end
    Tables.Schema(Tables.columnnames(m), col_types)
end

function Tables.schema(m::TimelineTableWithMissing)
    col_types = Type[]
    for col in names(m)
        if Symbol(col) == :date
            push!(col_types, Date)
        else
            push!(col_types, Union{Missing, Float64})
        end
    end
    Tables.Schema(names(m), col_types)
end

Tables.columns(x::TimelineTable) = x
function Tables.getcolumn(x::TimelineTable, i::Int)
    Tables.getcolumn(x, x.lookup[i])
end
function Tables.getcolumn(x::TimelineTable, nm::TimelineColumn)
    x[:, nm]
end
function Tables.getcolumn(x::TimelineTable, nm::Symbol)
    p = Meta.parse(nm |> string)
    if typeof(p) <: Symbol
        Tables.getcolumn(x, TimelineColumn(nm))
    elseif typeof(p) <: Expr
        if p.args[1] ∉ (:lag, :lead)
            throw("Arg must be lag or lead")
        end
        v = length(p.args) == 3 ? p.args[3] : 1
        if p.args[1] == :lag
            Tables.getcolumn(x, lag(p.args[2], v))
        else p.args[1] == :lead
            Tables.getcolumn(x, lead(p.args[2], v))
        end
    else
        throw("Unknown error")
    end
end

Tables.columnnames(x::TimelineTable) = Symbol.(String.(names(x)))


function Base.length(x::TimelineTableNoMissing)
    if (x.missing_bdays !== nothing && length(x.missing_bdays) == 0)
        calculate_missing_bdays!(x)
    end
    sub = x.missing_bdays === nothing ? 0 : length(x.missing_bdays)
    return bdayscount(x.calendar, x.dates.left, x.dates.right) - sub + isbday(x.calendar, x.dates.right)
end

function Base.length(x::TimelineTableWithMissing)
    bdayscount(x.calendar, x.dates.left, x.dates.right) + isbday(x.calendar, x.dates.right)
end

function Base.length(x::DataVector)
    sub = x.missing_bdays === nothing ? 0 : length(x.missing_bdays)
    return length(raw_values(x)) - sub
end

DataFrames.dropmissing(x::TimelineTableNoMissing) = x
DataFrames.allowmissing(x::TimelineTableWithMissing) = x

function DataFrames.dropmissing(x::TimelineTableWithMissing)
    TimelineTableNoMissing(
        x.parent,
        x.firmdata,
        x.dates,
        x.cols,
        x.missing_bdays
    )
end

function DataFrames.allowmissing(x::TimelineTableNoMissing)
    TimelineTableWithMissing(
        x.parent,
        x.firmdata,
        x.dates,
        x.cols,
        x.missing_bdays
    )
end

function DataFrames.select!(x::TimelineTable, cols::Vector)
    select!(x, TimelineColumn.(cols))
end
function DataFrames.select!(x::TimelineTable, cols::Vector{TimelineColumn})
    #@assert all([col ∈ x.colnames for col in cols]) "Not all columns are in the table"
    if sort(cols) != sort(names(x))
        x.missing_bdays = Int[]
    end
    x.cols = DictIndex(cols)
    x
end


raw_values(x::DataVector) = x.data
data_dates(x::DataVector) = x.dates
data_missing_bdays(x::DataVector) = x.missing_bdays

raw_values(x::RegressionCache) = x.data
data_dates(x::RegressionCache) = x.dates
data_missing_bdays(x::RegressionCache) = x.missing_bdays
