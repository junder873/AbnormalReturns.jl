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



mutable struct MarketData
    date_start::Date
    date_end::Date
    dates::Vector{Date}
    cols::Vector{String}
    data::Matrix{<:Real}
    MarketData() = new()
end

struct CrspMarketCalendar <: HolidayCalendar end

function BusinessDays.isholiday(::CrspMarketCalendar, dt::Date; calendar_days::Vector{Date}=MARKET_DATA_CACHE.dates)
    dt ∉ calendar_days
end

function BusinessDays.isbday(hc::CrspMarketCalendar, dt::Date)::Bool
    if BusinessDays._getcachestate(hc)
        return isbday(BusinessDays._getholidaycalendarcache(hc), dt)
    else
        return !isholiday(hc, dt)
    end
end

struct FirmData
    date_start::Date
    date_end::Date
    missing_bday::Union{Nothing, Vector{UInt32}}
    data::Vector{<:Real}
    function FirmData(d1, d2, missing_bday, data)
        @assert d2 >= d1 "End Date must be greater than Start Date"
        @assert missing_bday === nothing || issorted(missing_bday) "Vector of missing days must be sorted"
        new(d1, d2, missing_bday, data)
    end
end

const FIRM_DATA_CACHE = Dict{Int, Dict{String, FirmData}}()
const MARKET_DATA_CACHE = MarketData()
const Warn_Settings = Dict{String, Bool}(
    "firm_dates_range" => false,
    "market_dates_range" => true
)

function update_market_data!(d1, d2, dates, cols, data)
    @assert d2 >= d1 "End Date must be greater than Start Date"
    @assert all(d1 .<= dates .<= d2) "All dates must be between start and end date"
    @assert length(cols) == size(data)[2] "Number of column names must be same as columns of matrix"
    @assert length(dates) == size(data)[1] "Number of rows must be same as dates"
    @assert unique(cols) == cols "Column names must be unique"

    MARKET_DATA_CACHE.date_start = d1
    MARKET_DATA_CACHE.date_end = d2
    MARKET_DATA_CACHE.dates = dates
    MARKET_DATA_CACHE.cols = cols
    MARKET_DATA_CACHE.data = data
end

"""
    FirmData(
        df::DataFrame;
        date_col="date",
        id_col="permno",
        valuecols=nothing
    )

Creates the cached data for firms. Does some initial cleaning (dropmissing and sort)
and then assigns the values to a Dictionary (typically with Permno as a key) with the values
as vectors in a second Dictionary.
"""
function FirmData(
    df::DataFrame;
    date_col="date",
    id_col="permno",
    valuecols=nothing
)
    if valuecols === nothing
        valuecols = [n for n in names(df) if n ∉ [date_col, id_col]]
    end
    if !isdefined(MARKET_DATA_CACHE, :data)
        temp = unique(df[:, [date_col]])
        MarketData(temp; date_col)
    end
    df = dropmissing(df, [id_col, date_col])
    check_all_businessdays(unique(df[:, date_col]))
    if any(nonunique(df, [id_col, date_col]))
        @error("There are duplicate id-date rows in the dataframe")
    end

    df = select(df, vcat([id_col, date_col], string.(valuecols)))
    dropmissing!(df, [id_col, date_col])
    sort!(df)


    for col in string.(valuecols)
        temp = select(df, [id_col, date_col, col])
        dropmissing!(temp)
        
        gdf = groupby(temp, id_col)
        for (i, g) in enumerate(gdf)
            if !haskey(FIRM_DATA_CACHE, keys(gdf)[i][1])
                FIRM_DATA_CACHE[keys(gdf)[i][1]] = Dict{String, FirmData}()
            end
            FIRM_DATA_CACHE[keys(gdf)[i][1]][col] = FirmData(
                g[1, date_col],
                g[end, date_col],
                find_missing_bday(g[:, date_col]; perform_checks=false),
                g[:, col]
            )
        end
    end
end

"""
    MarketData(
        df::DataFrame;
        date_col="date",
        valuecols=nothing,
        add_intercept=true,
        force_update=false
    )

Creates the cached data for the market. Does some initial cleaning and creates a timeline
of the existing data.
"""
function MarketData(
    df::DataFrame;
    date_col="date",
    valuecols=nothing,
    add_intercept=true,
    force_update=true
)
    if (
        !isdefined(MARKET_DATA_CACHE, :data) ||
        force_update ||
        (
            MARKET_DATA_CACHE.date_start <= minimum(df[:, date_col]) &&
            MARKET_DATA_CACHE.date_end >= maximum(df[:, date_col])
        )
    )
        
        if valuecols === nothing
            valuecols = [n for n in names(df) if n ≠ date_col]
        end
        df = dropmissing(df, vcat([date_col], valuecols))
        sort!(df, [date_col])
        if add_intercept
            df[!, :intercept] = ones(Int, nrow(df))
            valuecols=vcat(["intercept"], valuecols)
        end
        select!(df, vcat([date_col], valuecols))

        update_market_data!(
            minimum(df[:, date_col]),
            maximum(df[:, date_col]),
            df[:, date_col],
            string.(valuecols),
            Matrix(df[:, valuecols])
        )
        if BusinessDays._getcachestate(BusinessDays.symtocalendar(:CrspMarketCalendar))
            BusinessDays.cleancache("CrspMarketCalendar")
        end

        BusinessDays.initcache("CrspMarketCalendar", minimum(df[:, date_col]), maximum(df[:, date_col]))
    end
end

"""
    data_range(timed_data::TimelineData, d1::Date, d2::Date)

    data_range(firm_data::FirmData, mkt_data::MarketData, d1::Date, d2::Date)

Provides the data range based on the dates provided. If only one TimelineData type is provided,
returns a range (i.e., 15:60) that corresponds to the data stored in the type. If firm and market
data are provided, then returns a vector of integers to make sure the number of rows provided
in the market matrix matches the length of the data vector for the firm.
"""
function data_range(data_start::Date, d1::Date, d2::Date, not_range::Union{Nothing, AbstractVector}=nothing)
    s = bdayscount("CrspMarketCalendar", data_start, d1) + 1
    e = s + bdayscount("CrspMarketCalendar", d1, d2) - !isbday("CrspMarketCalendar", d2)
    if not_range !== nothing
        s -= sum(s .> not_range)
        e -= sum(e .>= not_range)
    end
    s:e
end

function data_range(d1::Date, d2::Date, firm_match::Union{Nothing, FirmData}=nothing)
    x = data_range(MARKET_DATA_CACHE.date_start, d1, d2)
    # need to fix, not is just based off of a different range than the date vector
    if firm_match !== nothing && firm_match.missing_bday !== nothing
        not_range = firm_match.missing_bday .+ bdayscount("CrspMarketCalendar", MARKET_DATA_CACHE.date_start, firm_match.date_start)
        #println(bdayscount("CrspMarketCalendar", MARKET_DATA_CACHE.date_start, firm_match.date_start))
        #println(x)
        not_range = not_range[x[1] .<= not_range .<= x[end]]
        not_range .-= (x[1] - 1)
        #println(not_range)
        x[Not(not_range)]
    else
        x
    end
end


function check_all_businessdays(dates)
    bday_list = isbday("CrspMarketCalendar", dates)
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

function find_missing_bday(dates::Vector{Date}; perform_checks=true)
    if perform_checks
        if !issorted(dates)
            dates = sort(dates)
        end
        check_all_businessdays(dates)
    end

    # since all days are unique to a firm, if there are no missing dates skip this whole process
    if length(dates) == bdayscount("CrspMarketCalendar", dates[1], dates[end]) + 1
        return nothing
    end
    out = UInt32[]
    dates_counter = 1
    for (i, d) in enumerate(MARKET_DATA_CACHE.dates[data_range(MARKET_DATA_CACHE.date_start, dates[1], dates[end])])
        if dates[dates_counter] > d
            push!(out, i)
        else
            dates_counter += 1
        end
    end
    out
end


"""
    get_firm_data(id::Real, date_start::Date, date_end::Date, col::String="ret")

Fetches a vector from the FIRM_DATA_CACHE for a specific firm over a date range.
"""
function get_firm_data(id::Real, date_start::Date, date_end::Date, col::String="ret")
    if !haskey(FIRM_DATA_CACHE, id)
        @warn("Data for firm id $id is not stored in the cached data")
        return [missing]
    end
    firm_data = FIRM_DATA_CACHE[id][col]
    if Warn_Settings["firm_dates_range"]
        date_start < firm_data.date_start && @warn "Minimum Date is less than Cached Firm Date Start, this will be adjusted."
        date_end > firm_data.date_end && @warn "Maximum Date is greater than Cached Firm Date End, this will be adjusted."
    end
    date_start = max(date_start, firm_data.date_start)
    date_end = min(date_end, firm_data.date_end)

    if date_end < date_start
        return [missing]
    end

    firm_data.data[data_range(firm_data.date_start, date_start, date_end, firm_data.missing_bday)]
end

col_pos(x::String, cols::Vector{String}) = findfirst(isequal(x), cols)

# only market data
"""
    get_market_data(date_start::Date, date_end::Date, cols_market::String...)

    get_market_data(id::Real, date_start::Date, date_end::Date, cols_market::Union{Nothing, Vector{String}}=nothing)

Fetches a Matrix of market data between two dates, if an id (Integer) is provided, then the rows of the matrix will be the same
length as the length of the vector for the firm betwee those two dates.
"""
function get_market_data(date_start::Date, date_end::Date, cols_market::String...)
    if Warn_Settings["market_dates_range"]
        date_start < MARKET_DATA_CACHE.date_start && @warn "Minimum Date is less than Cached Market Date Start, this will be adjusted."
        date_end > MARKET_DATA_CACHE.date_end && @warn "Maximum Date is greater than Cached Market Date End, this will be adjusted."
    end
    date_start = max(date_start, MARKET_DATA_CACHE.date_start)
    date_end = min(date_end, MARKET_DATA_CACHE.date_end)

    if length(cols_market) == 0
        pos = 1:length(MARKET_DATA_CACHE.cols)
    else
        @assert all([c in MARKET_DATA_CACHE.cols for c in cols_market]) "Not all columns are in the data"
        pos = [col_pos(c, MARKET_DATA_CACHE.cols) for c in cols_market]
    end
    MARKET_DATA_CACHE.data[data_range(MARKET_DATA_CACHE.date_start, date_start, date_end), pos]
end

# market data with same length as a firm data

function get_market_data(id::Real, date_start::Date, date_end::Date, cols_market::String...)
    if !haskey(FIRM_DATA_CACHE, id)
        @warn("Data for firm id $id is not stored in the cached data")
        return [missing]
    end
    firm_data = FIRM_DATA_CACHE[id]
    if Warn_Settings["firm_dates_range"]
        date_start < firm_data.date_start && @warn "Minimum Date is less than Cached Firm Date Start, this will be adjusted."
        date_end > firm_data.date_end && @warn "Maximum Date is greater than Cached Firm Date End, this will be adjusted."
    end
    if Warn_Settings["market_dates_range"]
        date_start < MARKET_DATA_CACHE.date_start && @warn "Minimum Date is less than Cached Market Date Start, this will be adjusted."
        date_end > MARKET_DATA_CACHE.date_end && @warn "Maximum Date is greater than Cached Market Date End, this will be adjusted."
    end
    date_start = max(date_start, firm_data.date_start, MARKET_DATA_CACHE.date_start)
    date_end = min(date_end, firm_data.date_end, MARKET_DATA_CACHE.date_end)

    if date_end < date_start
        return [missing]
    end

    if length(cols_market) == 0
        pos = 1:length(MARKET_DATA_CACHE.cols)
    elseif length(cols_market) == 1
        @assert cols_market[1] ∈ MARKET_DATA_CACHE.cols "$(cols_market[1]) is not in the market data cache"
        pos = col_pos(cols_market[1], MARKET_DATA_CACHE.cols)
    else
        @assert all([c in MARKET_DATA_CACHE.cols for c in cols_market]) "Not all columns are in the data"
        pos = [col_pos(c, MARKET_DATA_CACHE.cols) for c in cols_market]
    end
    MARKET_DATA_CACHE.data[data_range(date_start, date_end, firm_data), pos]
end

# market data and firm data with same length
"""
    get_firm_market_data(
        id::Real,
        date_start::Date,
        date_end::Date;
        cols_market::Union{Nothing, Vector{String}, String}=nothing,
        col_firm::String
    )

Returns a Tuple of a vector of firm data and a matrix (or vector if only a String is passed for
cols_market) of market data (matrix has same number of rows as vector length).
"""
function get_firm_market_data(
    id::Real,
    date_start::Date,
    date_end::Date;
    cols_market::Union{Nothing, Vector{String}, String}=nothing,
    col_firm::String="ret"
)
    if !haskey(FIRM_DATA_CACHE, id)
        @warn("Data for firm id $id is not stored in the cached data")
        return ([missing], [missing])
    end
    firm_data = FIRM_DATA_CACHE[id][col_firm]
    if Warn_Settings["firm_dates_range"]
        date_start < firm_data.date_start && @warn "Minimum Date is less than Cached Firm Date Start, this will be adjusted."
        date_end > firm_data.date_end && @warn "Maximum Date is greater than Cached Firm Date End, this will be adjusted."
    end
    if Warn_Settings["market_dates_range"]
        date_start < MARKET_DATA_CACHE.date_start && @warn "Minimum Date is less than Cached Market Date Start, this will be adjusted."
        date_end > MARKET_DATA_CACHE.date_end && @warn "Maximum Date is greater than Cached Market Date End, this will be adjusted."
    end
    date_start = max(date_start, firm_data.date_start, MARKET_DATA_CACHE.date_start)
    date_end = min(date_end, firm_data.date_end, MARKET_DATA_CACHE.date_end)

    if date_end < date_start
        return ([missing], [missing])
    end

    if cols_market === nothing
        pos = 1:length(MARKET_DATA_CACHE.cols)
    elseif typeof(cols_market) <: Vector
        @assert all([c in MARKET_DATA_CACHE.cols for c in cols_market]) "Not all columns are in the data"
        pos = [col_pos(c, MARKET_DATA_CACHE.cols) for c in cols_market]
    else
        @assert cols_market ∈ MARKET_DATA_CACHE.cols "$cols_market is not in the MARKET_DATA_CACHE"
        pos = col_pos(cols_market, MARKET_DATA_CACHE.cols)
    end
    (
        firm_data.data[data_range(firm_data.date_start, date_start, date_end, firm_data.missing_bday)],
        MARKET_DATA_CACHE.data[data_range(date_start, date_end, firm_data), pos]
    )
end

function clear_firm_cached_data!()
    for key in keys(FIRM_DATA_CACHE)
        delete!(FIRM_DATA_CACHE, key)
    end
end

firm_in_cache(id::Real) = haskey(FIRM_DATA_CACHE, id)