# The idea for this is heavily inspired by BusinessDays.jl
# Thank you to that package for making this possible

###################################################################
##
# Another idea to make this section better is to rewrite it relying
# on a single master calendar from BusienssDays (likely a custom
# calendar). This would have the advantage of being a more common
# interface and considerably reducing storage since the individual
# firms would just need date_start, date_end, and a vector of
# days that are missing relative to the busienssday calendar.
##
# The difficulty of ths is quickly getting the underlying data.
# There would be no loss of efficiency for firms that have no
# missing data, but for firms with missing dates then bdayscount
# would provide an innacurate position in the data and would need
# to be adjusted based on the missing dates, which could be slow.
##
###################################################################

struct DataMatrix
    cols::Index
    data::Matrix{<:Real}
    missing_bdays::Union{Nothing, Dict{Symbol, Set{Int}}}
    dt_min::Date
    dt_max::Date
    cal::HolidayCalendar
end

struct MarketData
    calendar::MarketCalendar
    colmapping::Dict{Symbol, Symbol}
    marketdata::DataMatrix
    firmdata::Dict{Int, DataMatrix}
end

struct CombinedData # this is the component that is tables.jl compatible
    cols::Index
    data::Matrix{<:Real}
    dates::Vector{Date}
    function Combinedata(cols, data)
        @assert length(cols) == size(data, 2) "Dimensions are mismatched"
        new(cols, data)
    end
    function Combinedata(cols, data, dates)
        @assert length(dates) == size(data, 1) "Dimensions are mismatched"
        @assert length(cols) == size(data, 2) "Dimensions are mismatched"
        new(cols, data, dates)
    end
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
        cal
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
            cal
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


# function Base.merge(data::FirmDataVector...)

"""
    get_firm_data(id::Real, date_start::Date, date_end::Date, col::Symbol="ret")

Fetches a vector from the FIRM_DATA_CACHE for a specific firm over a date range.
"""
function Base.getindex(data::DataMatrix, dates::StepRange{Date, <:DatePeriod}, cols::Vector{Symbol})
    if any([!haskey(data.cols, x) for x in cols])
        throw("Not all columns are in the data")
    end

    @assert data.dt_min <= dates[1] "Dates out of bounds"
    @assert data.dt_max >= dates[end] "Dates out of bounds"

    s = bdayscount(data.cal, data.dt_min, dates[1]) + 1
    e = s + bdayscount(data.cal, dates[1], dates[end]) - !isbday(data.cal, dates[end])
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
    DataMatrix(
        Index(cols),
        data.data[s:e, data.cols[cols]],
        new_missings,
        dates[1],
        dates[end],
        data.cal
    )
end

function Base.getindex(data::MarketData, id::Real, dates::StepRange{Date, <:DatePeriod}, cols::Vector{Symbol})
    if !haskey(data.firmdata, id)
        throw("Data for firm id $id is not stored in the data")
    end
    if any([!haskey(data.colmapping, x) for x in cols])
        throw("Not all columns are in the data")
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
        return data.marketdata[dates, cols_market]
    else
        firm_data = data.firmdata[id]
        dt_min = max(dates[1], firm_data.dt_min)
        dt_max = min(dates[end], firm_data.dt_max)
        if length(cols_market) == 0
            return firm_data[dt_min:Day(1):dt_max, cols_firm]
        end

        return merge(
            firm_data[dt_min:Day(1):dt_max, cols_firm],
            data.marketdata[dates, cols_market]
        )
    end
end

function Base.merge(data1::DataMatrix, data2::DataMatrix)
    @assert data1.dt_min == data2.dt_min "Dates must match to merge"
    @assert data1.dt_max == data2.dt_max "Dates must match to merge"
    @assert size(data1.data, 1) == size(data2.data, 1) "Length of matrices must match"
    cols = vcat(data1.cols.names, data2.cols.names)
    new_missings = if data1.missing_bdays === nothing && data2.missing_bdays === nothing
        nothing
    elseif data1.missing_bdays === nothing
        data2.missing_bdays
    elseif data2.missing_bdays === nothing
        data1.missing_bdays
    else
        merge(data1.missing_bdays, data2.missing_bdays)
    end
    DataMatrix(
        Index(cols),
        hcat(data1.data, data2.data),
        new_missings,
        data1.dt_min,
        data2.dt_max,
        data1.cal
    )
end
    
    
function Tables.MatrixTable(data::DataMatrix)
    data = if data.missing_bdays === nothing
        data.data
    else
        data.data[Not(union(values(data.missing_bdays))), :]
    end
    MatrixTable(
        data.cols.names,
        data.cols.lookup,
        data
    )
end
