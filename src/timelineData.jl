

struct DataVector
    vec::Vector{Float64}
    missing_bdays::Union{Nothing, Set{Int}}
    dates::ClosedInterval{Date}
end

struct MarketData
    calendar::MarketCalendar
    colmapping::Dict{Symbol, Symbol}
    marketdata::Dict{Symbol, DataVector}
    firmdata::Dict{Int, Dict{Symbol, DataVector}}
end

abstract type TimelineTable end

mutable struct TimlineTableNoMissing <: TimelineTable
    parent::MarketData
    firmdata::Union{Nothing, Dict{Symbol, DataVector}}
    dates::ClosedInterval{Date}
    colnames::Vector{Symbol}
    lookup::Dict{Int, Symbol}
    missing_bdays::Union{Nothing, Vector{Int}}# nothing if no missing bdays, empty vector if never calculated
end

mutable struct TimelineTableWithMissing <: TimelineTable
    parent::MarketData
    firmdata::Union{Nothing, Dict{Symbol, DataVector}}
    dates::ClosedInterval{Date}
    colnames::Vector{Symbol}
    lookup::Dict{Int, Symbol}# reverse of other itmes since data is just stored in vectors
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

function date_range(cal::MarketCalendar, data::DataVector, dt_min::Date, dt_max::Date)
    dt_min = max(data.dates.left, dt_min)
    dt_max = min(data.dates.right, dt_max)
    s = bdayscount(cal, data.dt_min, dt_min) + 1
    e = s + bdayscount(cal, dt_min, dt_max) - !isbday(cal, dt_max)
    s:e
end

###########################################################
# Access functions for multiple columns and range of dates
# returns another data matrix
###########################################################


function Base.getindex(
    data::MarketData,
    id::Int,
    dates::ClosedInterval{Date},
    cols::Vector{Symbol}
    )
    if !haskey(data.firmdata, id)
        throw("Data for firm id $id is not stored in the data")
    end
    if any([!haskey(data.colmapping, x) for x in cols if x != :date])
        throw("Not all columns are in the data")
    end

    return TimelineTableNoMissing(
        data,
        data.firmdata[id],
        dates,
        cols,
        Dict{Int, Symbol}(1:length(cols) .=> cols),
        Int[]
    )
end

function Base.getindex(data::MarketData, id::Int, ::Colon, cols::Vector{Symbol})
    if !haskey(data.firmdata, id)
        throw("Data for firm id $id is not stored in the data")
    end
    if any([!haskey(data.colmapping, x) for x in cols if x != :date])
        throw("Not all columns are in the data")
    end

    firmdata = data.firmdata[id]
    dates = if data.colmapping[col] == :firmdata
        minimum([x.dates.left for x in values(firmdata)]) .. maximum([x.dates.right for x in values(firmdata)])
    else
        data.cal.dt_min .. data.cal.dt_max
    end
    TimelineTableNoMissing(
        data,
        firmdata,
        dates,
        cols,
        Dict{Int, Symbol}(1:length(cols) .=> cols),
        Int[]
    )
end

function Base.getindex(data::MarketData, id::Int, dates::ClosedInterval{Date}, ::Colon)
    data[id, dates, vcat([:date], Symbol.(keys(data.colmapping)))]
end

function Base.getindex(data::MarketData, id::Int, ::Colon, ::Colon)
    data[id, :, vcat([:date], Symbol.(keys(data.colmapping)))]
end

function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, cols::Vector{Symbol})
    if any([!haskey(data.colmapping, x) for x in cols if x != :date])
        throw("Not all columns are in the data")
    end
    if any([data.colmapping[x] for x in cols] .== :firmdata)
        throw("An ID must be supplied to access columns of firm data")
    end

    TimelineTable(
        data,
        nothing,
        dates,
        cols,
        Dict{Int, Symbol}(1:length(cols) .=> cols),
        Int[]
    )
end



function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, cols::Vector{Symbol})
    data.dates = dates
    data.colnames = cols
    data.lookup = Dict{Int, Symbol}(1:length(cols) .=> cols)
    if !(data.dates == dates && sort(cols) == sort(data.cols))
        data.missing_bdays = Int[]
    end
    return data
end

function Base.getindex(data::TimelineTable, ::Colon, cols::Vector{Symbol})
    data[data.dates, cols]
end

function Base.getindex(data::TimelineTable, dates::ClosedInterval{Date}, ::Colon)
    data[dates, data.colnames]
end

function Base.getproperty(obj::TimelineTable, sym::Symbol)
    if sym == :calendar
        obj.parent.calendar
    elseif sym == :marketdata
        obj.parent.marketdata
    elseif sym == :colmapping
        obj.parent.colmapping
    else
        getfield(obj, sym)
    end
end

function Base.getindex(obj::TimelineTableNoMissing, ::Colon, nm::Symbol)::Vector{Float64}
    dt_min = max(obj.dates.left, [x[col].dates.left for col in names(x) if x != :date]...)
    dt_max = min(obj.dates.right, [x[col].dates.right for col in names(x) if x != :date]...)

    if nm == :date
        return listbdays(obj.calendar, dt_min, dt_max)
    end
    vec = getproperty(obj, obj.colmapping[nm])[nm]
    local_dt_min = max(dt_min, vec.dates.left)
    local_dt_max = min(dt_max, vec.dates.right)
    s = bdayscount(obj.calendar, vec.dates.left, local_dt_min) + 1
    e = s + bdayscount(obj.calendar, local_dt_min, local_dt_max) - !isbday(obj.calendar, local_dt_max)
    out = vec.vec[s:e]
    if obj.missing_bdays !== nothing && length(obj.missing_bdays) == 0
        calculate_missing_bdays!(x)
    end
    if obj.missing_bdays === nothing
        return out
    else
        return out[Not(adjust_missing_bdays(obj.calendar, obj.missing_bdays, obj.dates, dt_min .. dt_max))]
    end
end

function Base.getindex(obj::TimelineTableWithMissing, ::Colon, nm::Symbol)::Vector{Union{Missing, Float64}}
    dt_min = obj.dates.left
    dt_max = obj.dates.right

    if nm == :date
        return listbdays(obj.calendar, dt_min, dt_max)
    end
    vec = getproperty(obj, obj.colmapping[nm])[nm]
    local_dt_min = max(dt_min, vec.dates.left)
    local_dt_max = min(dt_max, vec.dates.right)
    s = bdayscount(obj.calendar, vec.dates.left, local_dt_min) + 1
    e = s + bdayscount(obj.calendar, local_dt_min, local_dt_max) - !isbday(obj.calendar, local_dt_max)
    out = vec.vec[s:e]
    if obj.missing_bdays !== nothing && length(obj.missing_bdays) == 0
        calculate_missing_bdays!(x)
    end
    to_add_pre = ifelse(
        dt_min < local_dt_min,
        bdayscount(obj.calendar, dt_min, local_dt_min),
        0
    )
    to_add_post = ifelse(
        dt_max > local_dt_max,
        bdayscount(obj.calendar, local_dt_max, dt_max),
        0
    )
    if vec.missing_bdays !== nothing
        out = allowmissing(out)
        out[adjust_missing_bdays(obj.cal, vec.missing_bdays, vec.dates, local_dt_min .. local_dt_max)] .= missing
    end
    if to_add_pre > 0 || to_add_post > 0
        out = vcat(
            repeat([missing], to_add_pre),
            out,
            repeat([missing], to_add_post)
        )
    end
    return out
end

function calculate_missing_bdays!(obj::TimelineTable)
    col_dt_min = maximum([obj[col].dates.left for col in names(obj) if col != :date])
    col_dt_max = minimum([obj[col].dates.right for col in names(obj) if col != :date])
    out = Set{Int}()
    if obj.dates.left < col_dt_min
        union!(out, Set(1:bdayscount(obj.calendar, obj.dates.left, col_dt_min)))
    end
    if obj.dates.right > col_dt_max
        r = Set(
            1:bdayscount(obj.calendar, col_dt_max, obj.dates.right) .+
            bdayscount(obj.calendar, obj.dates.left, col_dt_max)
        )
        union!(out, r)
    end
    for col in names(obj)
        col == :date && continue
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

Base.names(x::TimelineTable) = x.colnames

Tables.istable(::Type{<:TimelineTable}) = true
Tables.columnaccess(::Type{<:TimelineTable}) = true
function Tables.schema(m::TimelineTableNoMissing)
    col_types = Type[]
    for col in names(m)
        if col == :date
            push!(col_types, Date)
        else
            push!(col_types, Float64)
        end
    end
    Tables.Schema(names(m), col_types)
end

function Tables.schema(m::TimelineTableWithMissing)
    col_types = Type[]
    for col in names(m)
        if col == :date
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
function Tables.getcolumn(x::TimelineTable, nm::Symbol)
    x[:, nm]
end
Tables.columnnames(x::TimelineTable) = names(x)

function Base.length(x::TimelineTableNoMissing)
    if (x.missing_bdays !== nothing && length(x.missing_bdays) == 0)
        calculate_missing_bdays!(x)
    end
    sub = x.missing_bdays === nothing ? 0 : length(x.missing_bdays)
    return bdayscount(x.calendar, x.dates.left, x.dates.right) - sub + 1
end

function Base.length(x::TimelineTableWithMissing)
    bdayscount(x.calendar, x.dates.left, x.dates.right) + 1
end

DataFrames.dropmissing(x::TimelineTableNoMissing) = x
DataFrames.allowmissing(x::TimelineTableWithMissing) = x

function DataFrames.dropmissing(x::TimelineTableWithMissing)
    TimelineTableNoMissing(
        x.parent,
        x.firmdata,
        x.dates,
        x.colnames,
        x.lookup,
        x.missing_bdays
    )
end

function DataFrames.allowmissing(x::TimelineTableNoMissing)
    TimelineTableWithMissing(
        x.parent,
        x.firmdata,
        x.dates,
        x.colnames,
        x.lookup,
        x.missing_bdays
    )
end

function DataFrames.select!(x::TimelineTable, cols::Vector{Symbol})
    @assert all([col ∈ x.colnames for col in cols]) "Not all columns are in the table"
    x.colnames = cols
    x.lookup = Dict{Int, Symbol}(1:length(cols) .=> cols)
    if sort(cols) != sort(x.colnames)
        x.missing_bdays = Int[]
    end
    x
end
