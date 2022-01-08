
Statistics.var(rr::RegressionModel) = rss(rr) / dof_residual(rr)
Statistics.std(rr::RegressionModel) = sqrt(var(rr))

# it seems like I want some kind of formula with a similar minobs option for a lot
# of these functions as well.
function Statistics.var(
    data::TimelineTable,
    col_firm::Symbol,
    col_market::Symbol;
    minobs::Real=0.8
)
    if minobs < 1
        minobs = bdayscount(data.cal, data.dt_min, data.dt_max) * minobs
    end
    m = dropmissing(data[:, [col_firm, col_market]]).matrix
    if size(m, 1) < minobs
        return missing
    end
    var(m[:, 1] .- m[:, 2])
end

function Statistics.std(
    data::TimelineTable,
    col_firm::Symbol,
    col_market::Symbol;
    minobs::Real=0.8
)
    sqrt(var(data, col_firm, col_market; minobs))
end

Statistics.var(rr::Missing) = missing
Statistics.std(rr::Missing) = missing

bh_return(x) = prod(1 .+ skipmissing(x)) - 1
bhar(x, y) = bh_return(x) - bh_return(y)

# for firm data
"""
    bh_return(id::Int, date_start::Date, date_end::Date, col_firm::String="ret")
    bh_return(date_start::Date, date_end::Date, cols_market::String="vwretd")

Calculates the buy and hold returns (also called geometric return) for TimelineData. If an Integer
is passed, then it is calculated based on the FIRM_DATA_CACHE (for the integer provided), otherwise
is calculated for the MARKET_DATA_CACHE.

These functions treat missing returns in the period implicitly as a zero return.
"""
function bh_return(data::TimelineTable, col)
    bh_return(data[:, col])
end


"""
    bhar(id::Int, date_start::Date, date_end::Date, cols_market::String="vwretd", col_firm::String="ret")
    bhar(id::Int, date_start::Date, date_end::Date, rr::Union{Missing, RegressionModel})

Calculates the difference between buy and hold returns relative to the market. If a RegressionModel type is passed, then
the expected return is estimated based on the regression (Fama French abnormal returns). Otherwise, the value is
based off of the value provided (typically a market wide return).

These functions treat missing returns in the period implicitly as a zero return.
"""
function bhar(
    data::TimelineTable,
    firm_col::Symbol,
    mkt_col::Symbol;
    minobs::Real=0.8
)
    if minobs < 1
        minobs = bdayscount(data.cal, data.dt_min, data.dt_max) * minobs
    end
    select!(data, [firm_col, mkt_col])
    dropmissing!(data)
    if size(m, 1) < minobs
        return missing
    end
    bhar(m[:, firm_col], m[:, market_col])
end

function bhar(
    data::TimelineTable,
    rr::RegressionModel;
    minobs::Real=0.8
)
    if minobs < 1
        minobs = bdayscount(data.calendar, data.dates.left, data.dates.right) * minobs
    end
    select!(data, StatsModels.termvars(rr.formula))
    dropmissing!(data)
    if length(data) < minobs || length(data) <= length(data.colnames)
        return missing
    end
    sch = apply_schema(rr.formula, schema(rr.formula, data))
    resp, pred = modelcols(sch, data)
    bhar(resp, predict(rr, pred))

end
function bhar(
    data::TimelineTable,
    rr::Missing;
    minobs::Real=0.8
)
    return missing
end

"""
    car(id::Int, date_start::Date, date_end::Date, cols_market::String="vwretd", col_firm::String="ret")
    car(id::Int, date_start::Date, date_end::Date, rr::Union{Missing, RegressionModel})

Calculates the difference between cumulative returns relative to the market. If a RegressionModel type is passed, then
the expected return is estimated based on the regression (Fama French abnormal returns). Otherwise, the value is
based off of the value provided (typically a market wide return).

Cumulative returns are the simple sum of returns, they are often used due to their ease to calculate but
undervalue extreme returns compared to buy and hold returns (bh_return or bhar).

These functions treat missing returns in the period implicitly as a zero return.
"""
function car(
    data::TimelineTable,
    firm_col::Symbol,
    mkt_col::Symbol;
    minobs::Real=0.8
)
    if minobs < 1
        minobs = bdayscount(data.cal, data.dt_min, data.dt_max) * minobs
    end
    m = dropmissing(data[[firm_col, mkt_col]]).matrix
    if size(m, 1) < minobs
        return missing
    end
    sum(m[:, 1] .- m[:, 2])
end

function car(
    data::TimelineTable,
    rr::RegressionModel;
    minobs::Real=.8
)
    if minobs < 1
        minobs = bdayscount(data.cal, data.dt_min, data.dt_max) * minobs
    end
    resp, pred = dropmissing_modelcols(rr, data)
    if length(resp) < minobs
        return missing
    end
    sum(resp .- predict(rr, pred))
end


function get_coefficient_val(rr::RegressionModel, coefname::String...)
    for x in coefname
        if x ∈ coefnames(rr)
            return coef(rr)[findfirst(x .== coefnames(rr))]
        end
    end
    @error("None of $(coefname) is in the RegressionModel model.")
end
"""
    alpha(rr::RegressionModel, coefname::String...="intercept")

"alpha" in respect to the the CAPM model, i.e., the intercept in the model.
This is the alpha from the estimation period.

This function finds the position of the coefficient name provided, defaults to "intercept".
If the coefname is not in the regression, then this function returns an error.
"""
alpha(rr::RegressionModel, coefname::String...="(Intercept)") = get_coefficient_val(rr, coefname...)


"""
    beta(rr::RegressionModel, coefname::String...=["mkt", "mktrf", "vwretd", "ewretd"])

"beta" in respect to the CAPM model, i.e., the coefficient on the market return minus the risk free rate.
This is the beta from the estimation period.

This function finds the position of the coefficient name provided, defaults to several common market returns.
If the coefname is not in the regression, then this function returns an error.
"""
beta(rr::RegressionModel, coefname::String...=["mkt", "mktrf", "vwretd", "ewretd"]...) = get_coefficient_val(rr, coefname...)

alpha(rr::Missing, coefname::String...="(Intercept)") = missing
beta(rr::Missing, coefname::String...="error") = missing