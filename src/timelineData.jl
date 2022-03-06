

struct DataVector
    data::Vector{Float64}
    missing_bdays::Union{Nothing, Set{Int}}
    dates::ClosedInterval{Date}
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

struct AllowMissing{mssng} end

mutable struct TimelineTable{Mssng, T, MNames, FNames, N1, N2} <: Tables.AbstractColumns
    parent::MarketData{T, MNames, FNames, N1, N2}
    allow_missing::Type{AllowMissing{Mssng}}
    firmdata::Union{Nothing, NamedTuple{FNames, NTuple{N2, DataVector}}}
    dates::ClosedInterval{Date}
    cols::DictIndex
    missing_bdays::Union{Nothing, Vector{Int}}# nothing if no missing bdays, empty vector if never calculated
    produce_date::Bool# whether when converting to tables.jl interface to produce a date column
    # by assumption, this will be the first column produced
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

    market_data = NamedTuple(valuecols_market .=> DataVector.(Tables.columns(df_market[:, valuecols_market]), Ref(df_market[:, date_col_market])))

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
                    col_tab[date_col_firms][i]
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

function DataVector(data::AbstractVector, dates::AbstractVector{Date})
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
        missing_days = any(ismissing.(data)) ? Set(findall(ismissing, data)) : nothing
        data = coalesce.(data, zero(nonmissingtype(eltype(data))))
    else
        missing_days = nothing
    end
    
    
    DataVector(
        data,
        missing_days,
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
function dates_min_max(dates::ClosedInterval{Date}...)::ClosedInterval{Date}
    maximum([d.left for d in dates]) .. minimum([d.right for d in dates])
end

function adjust_missing_bdays(cal::MarketCalendar, missing_bdays::Vector{Int}, dates_missings::ClosedInterval{Date}, new_dates::ClosedInterval{Date})
    if new_dates.right < dates_missings.left || new_dates.left > dates_missings.right
        return nothing
    end
    #if new_dates.left <= dates_missings.left
    missing_bdays = missing_bdays .+ (bdayscount(cal, new_dates.left, dates_missings.left))
    e = bdayscount(cal, new_dates.left, new_dates.right) + isbday(cal, new_dates.right)
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
    date_range(cal, data.dates, dt_min .. dt_max)
end

function date_range(cal::MarketCalendar, data::DataVector, dates::ClosedInterval{Date})
    date_range(cal, data.dates, dates)
end

function Base.getindex(data::RegressionCache, dates::ClosedInterval{Date}, mssngs::Nothing)
    new_dates = dates_min_max(data.dates, dates)
    r = date_range(data.calendar, data.dates, new_dates)
    @view data.data[r, :]
end

function Base.getindex(data::RegressionCache, dates::ClosedInterval{Date}, mssngs::Vector{Int})
    new_dates = dates_min_max(data.dates, dates)
    new_mssngs = adjust_missing_bdays(data.calendar, mssngs, dates, new_dates)
    if new_mssngs === nothing
        return data[dates, nothing]
    end
    r = date_range(data.calendar, data.dates, new_dates)
    data.data[r[Not(new_mssngs)], :]
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

###########################################################
# Access functions for multiple columns and range of dates
# returns another data matrix
###########################################################

@inline function Base.getindex(data::MarketData{T, MNames, FNames}, col::TimelineColumn) where {T, MNames, FNames}
    if Symbol(col) ∈ FNames
        throw("An ID must be supplied to access columns of firm data")
    end

    out = data.marketdata[Symbol(col)]
    if col.shifts == 0
        return out
    end
    return shift(out, col.shifts, data.calendar)
end

@inline function Base.getindex(data::MarketData{T, MNames, FNames}, id::T, col::TimelineColumn) where {T, MNames, FNames}
    if Symbol(col) ∈ MNames
        return data[col]
    end
    out = data.firmdata[id][Symbol(col)]
    if col.shifts == 0
        return out
    end
    return shift(out, col.shifts, data.calendar)
end



@inline function Base.getindex(
    data::MarketData{T},
    id::T,
    dates::ClosedInterval{Date},
    cols::Vector
    ) where {T}
    data[id, dates, TimelineColumn.(cols)]
end
@inline @views function Base.getindex(
    data::MarketData{T, MNames, FNames},
    id::T,
    dates::ClosedInterval{Date},
    cols::Vector{TimelineColumn}
    ) where {T, MNames, FNames}
    if !haskey(data.firmdata, id)
        throw(ArgumentError("Data for firm id $id is not stored in the data"))
    end
    if !check_col(Symbol.(cols), MNames, FNames)
        throw(ArgumentError("Not all columns are in the data"))
    end

    return TimelineTable(
        data,
        AllowMissing{false},
        data.firmdata[id],
        dates,
        DictIndex(cols),
        Int[],
        true
    )
end

@inline function Base.getindex(data::MarketData{T}, id::T, ::Colon, cols::Vector) where {T}
    data[id, :, TimelineColumn.(cols)]
end
@inline function Base.getindex(data::MarketData{T, MNames, FNames}, id::T, ::Colon, cols::Vector{TimelineColumn}) where {T, MNames, FNames}
    if !haskey(data.firmdata, id)
        throw(ArgumentError("Data for firm id $id is not stored in the data"))
    end
    if !check_col(Symbol.(cols), MNames, FNames)
        throw(ArgumentError("Not all columns are in the data"))
    end

    firmdata = data.firmdata[id]
    dates = dates_min_max(data_dates.([data[id, col] for col in cols]))
    TimelineTable(
        data,
        AllowMissing{false},
        firmdata,
        dates,
        DictIndex(cols),
        Int[],
        true
    )
end

@inline function Base.getindex(data::MarketData{T, MNames, FNames}, id::T, dates::ClosedInterval{Date}, ::Colon) where {T, MNames, FNames}
    data[id, dates, TimelineColumn.([MNames..., FNames...])]
end

@inline function Base.getindex(data::MarketData{T, MNames, FNames}, id::T, ::Colon, ::Colon) where {T, MNames, FNames}
    data[id, :, TimelineColumn.([MNames..., FNames...])]
end

@inline function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, cols::Vector)
    data[dates, TimelineColumn.(cols)]
end
@inline function Base.getindex(data::MarketData{T, MNames, FNames}, dates::ClosedInterval{Date}, cols::Vector{TimelineColumn}) where {T, MNames, FNames}
    if !check_col(Symbol.(cols), MNames, FNames)
        throw(ArgumentError("Not all columns are in the data"))
    end
    if check_col(Symbol.(cols), FNames)
        throw(ArgumentError("An ID must be supplied to access columns of firm data"))
    end

    TimelineTable(
        data,
        AllowMissing{false},
        nothing,
        dates,
        DictIndex(cols),
        Int[],
        true
    )
end

@inline Base.getindex(data::MarketData, ::Colon, cols::Vector) = data[:, TimelineColumn.(cols)]
@inline function Base.getindex(data::MarketData, ::Colon, cols::Vector{TimelineColumn})
    dates = dates_min_max(data_dates.([data[col] for col in cols])...)
    data[dates, cols]
end

@inline function Base.getindex(data::MarketData{T, MNames, FNames}, ::Colon, ::Colon) where {T, MNames, FNames}
    data[:, [MNames...]]
end

@inline function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, cols::Vector)
    data[dates, TimelineColumn.(cols)]
end
@inline function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, cols::Vector{TimelineColumn})
    if !(data.dates == dates && sort(cols) == sort(names(data)))
        data.missing_bdays = Int[]
    end
    data.dates = dates
    data.cols = DictIndex(cols)
    return data
end

@inline function Base.getindex(data::TimelineTable, ::Colon, cols::Vector)
    data[:, TimelineColumn.(cols)]
end
@inline function Base.getindex(data::TimelineTable, ::Colon, cols::Vector{TimelineColumn})
    data[data.dates, cols]
end

@inline function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, ::Colon)
    data[dates, names(data)]
end

function Base.getproperty(obj::TimelineTable{Mssngs, T, MNames, FNames}, sym::Symbol) where {Mssngs, T, MNames, FNames}
    if sym == :calendar
        obj.parent.calendar
    elseif sym == :marketdata
        obj.parent.marketdata
    elseif sym == :regression_cache
        obj.parent.regression_cache
    elseif sym == :date
        get_dates(obj)
    elseif sym ∈ MNames
        getproperty(obj.parent.marketdata, sym)
    elseif sym ∈ FNames
        getproperty(obj.firmdata, sym)
    else
        getfield(obj, sym)
    end
end

function Base.getindex(obj::TimelineTable, nm::Symbol)
    obj[TimelineColumn(nm)]
end

function Base.getindex(obj::TimelineTable, nm::TimelineColumn)
    shift(getproperty(obj, Symbol(nm)), nm.shifts, obj.calendar)
end

function get_dates(obj::TimelineTable{false})
    out = listbdays(obj.calendar, obj.dates.left, obj.dates.right)
    if obj.missing_bdays !== nothing && length(obj.missing_bdays) == 0
        calculate_missing_bdays!(obj)
    end
    if obj.missing_bdays === nothing
        return out
    else
        return out[Not(obj.missing_bdays)]
    end
end

function get_dates(obj::TimelineTable{true})
    listbdays(obj.calendar, obj.dates.left, obj.dates.right)
end

function Base.getindex(obj::TimelineTable{false}, ::Colon, nm::Symbol)::Vector{Float64}
    obj[:, TimelineColumn(nm)]
end
@inline @views function Base.getindex(obj::TimelineTable{false}, ::Colon, nm::TimelineColumn)::Vector{Float64}
    dates = dates_min_max(obj.dates, [obj[col].dates for col in names(obj)]...)

    vec = obj[nm]
    local_dates = dates_min_max(dates, vec.dates)
    r = date_range(obj.calendar, vec, local_dates)
    out = vec.data[r]
    if obj.missing_bdays !== nothing && length(obj.missing_bdays) == 0
        calculate_missing_bdays!(obj)
    end
    if obj.missing_bdays === nothing
        return out
    else
        return out[Not(adjust_missing_bdays(obj.calendar, obj.missing_bdays, obj.dates, dates))]
    end
end

function Base.getindex(obj::TimelineTable{true}, ::Colon, nm::Symbol)::Vector{Union{Missing, Float64}}
    obj[:, TimelineColumn(nm)]
end
function Base.getindex(obj::TimelineTable{true}, ::Colon, nm::TimelineColumn)::Vector{Union{Missing, Float64}}

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
    col_dates = dates_min_max(obj.dates, [obj[col].dates for col in names(obj)]...)
    out = Set{Int}()
    if obj.dates.left < col_dates.left # set 1:first available dates as missing
        union!(out, Set(1:bdayscount(obj.calendar, obj.dates.left, col_dates.left)))
    end
    if obj.dates.right > col_dates.right # set last available date:end as missing
        r = Set(
            bdayscount(obj.calendar, obj.dates.left, col_dates.right) + 2:bdayscount(obj.calendar, obj.dates.left, obj.dates.right)+isbday(obj.calendar, obj.dates.right)
        )
        union!(out, r)
    end
    for col in names(obj)
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
function Tables.schema(obj::TimelineTable{mssngs}) where {mssngs}
    col_sym = Tables.columnnames(obj)
    col_types = fill(mssngs ? Union{Float64, Missing} : Float64, length(col_sym))
    if obj.produce_date
        col_types[1] = Date
    end
    Tables.Schema(col_sym, col_types)
end

Tables.columns(x::TimelineTable) = x
function Tables.getcolumn(x::TimelineTable, i::Int)
    if x.produce_date
        i -= 1
    end
    Tables.getcolumn(x, x.lookup[i])
end
function Tables.getcolumn(x::TimelineTable, nm::TimelineColumn)
    x[:, nm]
end
function Tables.getcolumn(x::TimelineTable{Mssngs, T, MNames, FNames}, nm::Symbol) where {Mssngs, T, MNames, FNames}
    if nm == :date
        return x.date
    end
    if nm ∈ MNames || nm ∈ FNames
        return x[:, TimelineColumn(nm)]
    end
    if nm ∈ Tables.columnnames(x)
        pos = findfirst(nm .== Tables.columnnames(x))
        if x.produce_date
            pos -= 1
        end
        return x[:, names(x)[pos]]
    end
    p = Meta.parse(nm |> string)
    if typeof(p) <: Symbol
        Tables.getcolumn(x, TimelineColumn(nm))
    elseif typeof(p) <: Expr
        if p.args[1] ∉ (:lag, :lead)
            throw(ArgumentError("Arg must be lag or lead"))
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

function Tables.columnnames(x::TimelineTable)
    out = Symbol.(String.(names(x)))
    x.produce_date ? vcat([:date], out) : out
end


function Base.length(x::TimelineTable{false})
    if (x.missing_bdays !== nothing && length(x.missing_bdays) == 0)
        calculate_missing_bdays!(x)
    end
    sub = x.missing_bdays === nothing ? 0 : length(x.missing_bdays)
    return bdayscount(x.calendar, x.dates.left, x.dates.right) - sub + isbday(x.calendar, x.dates.right)
end

function Base.length(x::TimelineTable{true})
    bdayscount(x.calendar, x.dates.left, x.dates.right) + isbday(x.calendar, x.dates.right)
end

function Base.length(x::DataVector)
    sub = x.missing_bdays === nothing ? 0 : length(x.missing_bdays)
    return length(raw_values(x)) - sub
end

DataFrames.dropmissing(x::TimelineTable{false}) = x
DataFrames.allowmissing(x::TimelineTable{true}) = x

function DataFrames.dropmissing(x::TimelineTable{true})
    TimelineTable(
        x.parent,
        AllowMissing{false},
        x.firmdata,
        x.dates,
        x.cols,
        x.missing_bdays,
        x.produce_date
    )
end

function DataFrames.allowmissing(x::TimelineTable{false})
    TimelineTable(
        x.parent,
        AllowMissing{true},
        x.firmdata,
        x.dates,
        x.cols,
        x.missing_bdays,
        x.produce_date
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

function Base.show(io::IO, data::MarketData{T, MNames, FNames}) where {T, MNames, FNames}
    println(io, "MarketData with ID type $T with $(length(data.firmdata)) unique firms")
    println(io, data.calendar)
    println(io, "Market Columns: $(join(MNames, ", "))")
    println(io, "Firm Columns: $(join(FNames, ", "))")
end