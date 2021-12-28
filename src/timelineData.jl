

struct DataMatrix{A <: Real}
    cols::Index
    matrix::Matrix{A}
    missing_bdays::Union{Nothing, Dict{Symbol, Set{Int}}}
    dt_min::Date
    dt_max::Date
    cal::HolidayCalendar
    date_col::Bool
end

struct MarketData
    calendar::MarketCalendar
    colmapping::Dict{Symbol, Symbol}
    marketdata::DataMatrix
    firmdata::Dict{Int, DataMatrix}
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
    start_time = now()
    df_market = DataFrame(df_market)
    df_firms = DataFrame(df_firms)
    if valuecols_market === nothing
        valuecols_market = Symbol.([n for n in Symbol.(names(df_market)) if n ∉ [date_col_market]])
    end
    if valuecols_firms === nothing
        valuecols_firms = Symbol.([n for n in Symbol.(names(df_firms)) if n ∉ [date_col_firms, id_col]])
    end

    println("1 ", Dates.canonicalize(now() - start_time))
    df_market = select(df_market, vcat([date_col_market], valuecols_market))
    dropmissing!(df_market)
    sort!(df_market)

    if any(nonunique(df_market, [date_col_market]))
        @error("There are duplicate date rows in the market data")
    end

    println("2 ", Dates.canonicalize(now() - start_time))
    df_firms = select(df_firms, vcat([id_col, date_col_firms], valuecols_firms))
    dropmissing!(df_firms, [id_col, date_col_firms])
    sort!(df_firms, [id_col, date_col_firms])

    println("3 ", Dates.canonicalize(now() - start_time))

    if any(nonunique(df_firms, [id_col, date_col_firms]))
        @error("There are duplicate id-date rows in the firm data")
    end

    cal = MarketCalendar(df_market[:, date_col_market])

    market_data = DataMatrix(
        Index(valuecols_market),
        Matrix(df_market[:, valuecols_market]),
        nothing,
        minimum(df_market[:, date_col_market]),
        maximum(df_market[:, date_col_market]),
        cal,
        false
    )

    println("4 ", Dates.canonicalize(now() - start_time))
    col_map = Dict{Symbol, Symbol}()
    for col in Symbol.(valuecols_market)
        col_map[col] = :marketdata
    end
    for col in Symbol.(valuecols_firms)
        col_map[col] = :firmdata
    end

    check_all_businessdays(unique(df_firms[:, date_col_firms]), cal)

    println("5 ", Dates.canonicalize(now() - start_time))

    gdf = groupby(df_firms, id_col)
    println("6 ", Dates.canonicalize(now() - start_time))
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

    println("7 ", Dates.canonicalize(now() - start_time))
    
    df_idx_base = combine(
        gdf,
        valuecols_firms .=> find_missings .=> valuecols_firms,
        date_col_firms .=> [minimum, maximum] .=> [:date_min, :date_max]
    )

    println("8 ", Dates.canonicalize(now() - start_time))
    for col in valuecols_firms
        df_firms[!, col] = coalesce.(
            df_firms[:, col],
            zeros(
                nonmissingtype(eltype(df_firms[:, col])),
                nrow(df_firms)
            )
        )
    end

    firm_matrix = Matrix(df_firms[:, valuecols_firms])



    println("9 ", Dates.canonicalize(now() - start_time))
    firm_data = Dict{Int, DataMatrix}()
    sizehint!(firm_data, nrow(df_idx_base))

    col_indx = Index(valuecols_firms)

    for i in 1:nrow(df_idx_base)
        idx = gdf.keymap[(df_idx_base[i, id_col],)]
        s = gdf.starts[idx]
        e = gdf.ends[idx]
        firm_data[df_idx_base[i, id_col]] = DataMatrix(
            col_indx,
            firm_matrix[s:e, :],
            row_to_dict(df_idx_base[i, valuecols_firms]),
            df_idx_base[i, :date_min],
            df_idx_base[i, :date_max],
            cal,
            false
        )
    end

    println("10 ", Dates.canonicalize(now() - start_time))

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

function find_missings(x)
    Set(
        findall(
            ismissing.(x)
        )
    )
end

function row_to_dict(df_row::DataFrameRow)
    if all(length.(values(df_row)) .== 0)
        nothing
    else
        Dict{Symbol, Set{Int}}(Symbol.(names(df_row)) .=> values(df_row))
    end
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


function adjust_missing_bdays(missing_bdays, s::Int, e::Int)
    Set(filter(x -> s <= x <= e, missing_bdays) .- (s - 1))
end

function adjust_missing_bdays(missing_bdays, s::Int)
    Set(missing_bdays .+ s)
end

function Base.vcat(data1::DataMatrix, data2::DataMatrix)
    @assert advancebdays(data1.cal, data1.dt_max, 1) == data2.dt_min "Dates in the two sets must be sequential"
    @assert data1.cols == data2.cols "Columns in the two sets must be the same"

    if data1.missing_bdays === nothing && data2.missing_bdays === nothing
        new_missings = nothing
    elseif data1.missing_bdays === nothing
        new_missings = data2.missing_bdays
        for key in keys(new_missings)
            new_missings[key] = adjust_missing_bdays(new_missings[key], length(data1))
        end
    elseif data2.missing_bdays === nothing
        new_missings = data1.missing_bdays
    else
        new_missings = data1.missing_bdays
        for key in keys(new_missings)
            new_missings[key] = union(new_missings[key], adjust_missing_bdays(data2.missing_bdays[key], length(data1)))
        end
    end
    DataMatrix(
        data1.cols,
        vcat(data1.matrix, data2.matrix),
        new_missings,
        data1.dt_min,
        data2.dt_max,
        data1.cal,
        data1.date_col || data2.date_col
    )
end

function date_range(data, dt_min, dt_max)
    dt_min = max(data.dt_min, dt_min)
    dt_max = min(data.dt_max, dt_max)
    s = bdayscount(data.cal, data.dt_min, dt_min) + 1
    e = s + bdayscount(data.cal, dt_min, dt_max) - !isbday(data.cal, dt_max)
    s:e
end

"""
    get_firm_data(id::Real, date_start::Date, date_end::Date, col::Symbol="ret")

Fetches a vector from the FIRM_DATA_CACHE for a specific firm over a date range.
"""

###########################################################
# Access functions for multiple columns and range of dates
# returns another data matrix
###########################################################
function Base.getindex(
    data::DataMatrix{T},
    dates::ClosedInterval{Date},
    cols::Vector{Symbol}
    ) where {T<:Real}
    if any([!haskey(data.cols, x) for x in cols if x != :date])
        println([x for x in cols if !haskey(data.cols, x)])
        throw("Not all columns are in the data")
    end
    if :date ∈ cols
        filter!(x -> x != :date, cols)
        date_col=true
    else
        date_col=false
    end

    dt_min = max(data.dt_min, dates.left)
    dt_max = min(data.dt_max, dates.right)
    s = bdayscount(data.cal, data.dt_min, dt_min) + 1
    e = s + bdayscount(data.cal, dt_min, dt_max) - !isbday(data.cal, dt_max)

    if data.missing_bdays !== nothing
        new_missings = Dict{Symbol, Set{Int}}()
        for col in cols
            new_missings[col] = adjust_missing_bdays(data.missing_bdays[col], s, e)
        end
        if all(length.(values(new_missings)) .== 0)
            new_missings = nothing
        end
    else
        new_missings = nothing
    end
    out = DataMatrix(
        Index(cols),
        data.matrix[s:e, data.cols[cols]],
        new_missings,
        dt_min,
        dt_max,
        data.cal,
        date_col
    )
    if data.dt_min > dates.left
        to_add_pre = bdayscount(data.cal, dates.left, data.dt_min)
        new_missings = Dict{Symbol, Set{Int}}(cols .=> repeat([Set(1:to_add_pre)], length(cols)))
        out = vcat(
            DataMatrix(
                out.cols,
                zeros(eltype(data.matrix), (to_add_pre, length(cols))),
                new_missings,
                dates.left,
                advancebdays(data.cal, dt_min, -1),
                data.cal,
                date_col
            ),
            out
        )
    end
    if data.dt_max < dates.right
        to_add_post = bdayscount(data.cal, data.dt_max, dates.right)
        new_missings = Dict{Symbol, Set{Int}}(cols .=> repeat([Set(1:to_add_post)], length(cols)))
        out = vcat(
            DataMatrix(
                out.cols,
                zeros(eltype(data.matrix), (to_add_post, length(cols))),
                new_missings,
                advancebdays(data.cal, dt_max, 1),
                dates.right,
                data.cal,
                date_col
            ),
            out
        )
    end
    out
end

function Base.getindex(data::DataMatrix, ::Colon, cols::Vector{Symbol})
    data[data.dt_min .. data.dt_max, cols]
end

function Base.getindex(data::DataMatrix{T}, dates::ClosedInterval{Date}, ::Colon) where {T<:Real}
    data[dates, data.cols.names]
end

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

    if :date ∈ cols
        filter!(x -> x != :date, cols)
        date_col=true
    else
        date_col=false
    end

    cols_market = Symbol[]
    cols_firm = Symbol[]
    for col in cols
        if data.colmapping[col] == :marketdata
            push!(cols_market, col)
        else
            push!(cols_firm, col)
        end
    end

    if length(cols_firm) == 0
        if date_col
            push!(cols_market, :date)
        end
        return data.marketdata[dates, cols_market]
    else
        firm_data = data.firmdata[id]
        if length(cols_market) == 0
            if date_col
                push!(cols_firm, :date)
            end
            return firm_data[dates, cols]
        end

        return hcat(
            firm_data[dates, cols_firm],
            data.marketdata[dates, cols_market];
            date_col
        )
    end
end

function Base.getindex(data::MarketData, id::Int, ::Colon, cols::Vector{Symbol})
    dates = if data.colmapping[col] == :firmdata
        data.firmdata[id].dt_min .. data.firmdata[id].dt_max
    else
        data.cal.dt_min .. data.cal.dt_max
    end
    data[id, dates, cols]
end

function Base.getindex(data::MarketData, id::Int, dates::ClosedInterval{Date}, ::Colon)
    data[id, dates, Symbol.(keys(data.colmapping))]
end

function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, cols::Vector{Symbol})
    if any([!haskey(data.colmapping, x) for x in cols if x != :date])
        throw("Not all columns are in the data")
    end
    if any([data.colmapping[x] for x in cols] .== :firmdata)
        throw("An ID must be supplied to access columns of firm data")
    end

    data.marketdata[dates, cols]
end

function Base.getindex(data::MarketData, ::Colon, cols::Vector{Symbol})
    data[data.cal.dt_min .. data.cal.dt_max, cols]
end

function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, ::Colon)
    data[dates, Symbol.(keys(data.colmapping))]
end

################################################
# Access functions for an individual columns
# returns a vector
################################################
function Base.getindex(
    data::DataMatrix{T},
    dates::ClosedInterval{Date},
    col::Symbol
    ) where {T<:Real}
    if !haskey(data.cols, col)
        throw("Column is not in the data")
    end

    r = date_range(data, dates.left, dates.right)
    out = Tables.getcolumn(data, col)[r]
    if data.dt_min > dates.left
        to_add_pre = bdayscount(data.cal, dates.left, data.dt_min)
        out = vcat(
            repeat([missing], to_add_pre),
            out
        )
    end
    if data.dt_max < dates.right
        to_add_post = bdayscount(data.cal, data.dt_max, dates.right)
        out = vcat(
            repeat([missing], to_add_post),
            out
        )
    end
    out
end

function Base.getindex(data::DataMatrix{T}, ::Colon, col::Symbol) where {T<:Real}
    data[data.dt_min .. data.dt_max, col]
end

function Base.getindex(data::MarketData,
    id::Int,
    dates::ClosedInterval{Date},
    col::Symbol
    )
    if !haskey(data.colmapping, col)
        throw("Column is not in the data")
    end
    if data.colmapping[col] == :firmdata
        return getindex(data.firmdata[id], dates, col)
    else
        return data.marketdata[dates, col]
    end
end

function Base.getindex(data::MarketData, dates::ClosedInterval{Date}, col::Symbol)
    if !haskey(data.colmapping, col)
        throw("Column is not in the data")
    end
    if data.colmapping[col] == :firmdata
        throw("To get firm data an ID must be supplied")
    else
        return data.marketdata[dates, col]
    end
end

function Base.getindex(data::MarketData, id::Int, ::Colon, col::Symbol)
    dates = if data.colmapping[col] == :firmdata
        data.firmdata[id].dt_min .. data.firmdata[id].dt_max
    else
        data.cal.dt_min .. data.cal.dt_max
    end
    data[id, dates, col]
end

function Base.getindex(data::MarketData, ::Colon, col::Symbol)
    data[data.cal.dt_min .. data.cal.dt_max, col]
end

#####################################################
# Access for single date/column
#####################################################

function Base.getindex(data::DataMatrix{T}, d::Date, col::Symbol) where {T<:Real}
    if !haskey(data.cols, col)
        throw("Column $col is not in the data")
    end

    if !isbday(data.cal, d) || d < data.dt_min || data.dt_max < d
        throw("$d is not a businessday in the data")
    end

    i = bdayscount(data.cal, data.dt_min, d)
    if data.missing_bdays === nothing || i ∉ data.missing_bdays[col]
        return data.matrix[i, data.cols[col]]
    else
        return missing
    end
end

function Base.getindex(data::MarketData, id::Int, d::Date, col::Symbol)
    if !haskey(data.colmapping, col)
        throw("Column $col is not in the data")
    end

    if data.colmapping[col] == :firmdata
        return data.firmdata[id][d, col]
    else
        return data.marketdata[d, col]
    end
end

function Base.getindex(data::MarketData, d::Date, col::Symbol)
    if !haskey(data.colmapping, col)
        throw("Column $col is not in the data")
    end

    if data.colmapping[col] == :firmdata
        throw("To get firm data, an ID must be supplied")
    else
        return data.marketdata[d, col]
    end
end



function Base.hcat(data1::DataMatrix{T}, data2::DataMatrix{T}; date_col=false) where {T<:Real}
    @assert data1.dt_min == data2.dt_min "Dates must match to merge"
    @assert data1.dt_max == data2.dt_max "Dates must match to merge"
    @assert size(data1.matrix, 1) == size(data2.matrix, 1) "Length of matrices must match"
    cols = merge(data1.cols, data2.cols)
    if data1.missing_bdays === nothing && data2.missing_bdays === nothing
        new_missings = nothing
    else
        new_missings1 = if data1.missing_bdays === nothing
            Dict{Symbol, Set{Int}}(data1.cols.names .=> repeat([Set{Int}()], length(data1.cols)))
        else
            data1.missing_bdays
        end
        new_missings2 = if data2.missing_bdays === nothing
            Dict{Symbol, Set{Int}}(data2.cols.names .=> repeat([Set{Int}()], length(data2.cols)))
        else
            data2.missing_bdays
        end
        new_missings = merge(new_missings1, new_missings2)
    end
    DataMatrix(
        cols,
        hcat(data1.matrix, data2.matrix),
        new_missings,
        data1.dt_min,
        data2.dt_max,
        data1.cal,
        date_col
    )
end

Base.values(x::DataMatrix) = x.matrix

Base.names(x::DataMatrix) = x.cols.names


Tables.istable(::Type{<:DataMatrix}) = true
Tables.columnaccess(::Type{<:DataMatrix}) = true
Tables.schema(m::DataMatrix{T}) where {T} = Tables.Schema(names(m), fill(eltype(T), size(values(m), 2)))

Tables.columns(x::DataMatrix) = x
function Tables.getcolumn(x::DataMatrix{T}, i::Int) where {T<:Real}
    col = x.cols.names[i]
    if x.missing_bdays === nothing || length(x.missing_bdays[col]) == 0
        x.matrix[:, i]
    else
        out = allowmissing(x.matrix[:, i])
        out[collect(x.missing_bdays[col])] .= missing
        out
    end
end
function Tables.getcolumn(x::DataMatrix{T}, nm::Symbol) where {T<:Real}
    if nm == :date
        listbdays(x.cal, x.dt_min, x.dt_max)
    else
        Tables.getcolumn(x, x.cols[nm])
    end
end
Tables.columnnames(x::DataMatrix) = x.date_col ? vcat([:date], names(x)) : names(x)

Base.length(x::DataMatrix) = size(x.matrix, 1)

function DataFrames.dropmissing(x::DataMatrix)
    if x.missing_bdays === nothing
        return x
    end
    drops = union(values(x.missing_bdays)...)
    new_data = x.matrix[Not(collect(drops)), :]
    DataMatrix(
        x.cols,
        new_data,
        nothing,
        x.dt_min,
        x.dt_max,
        x.cal,
        false
    )
end

# Tables.rowaccess(::Type{DataMatrix}) = true
# Tables.rows(x::DataMatrix) = x

# struct MatrixRow{T} <: Tables.AbstractRow
#     row::Int
#     source::DataMatrix{T}
#     date::Date
# end

# Base.eltype(x::DataMatrix{T}) where {T} = DataMatrixRow{T}