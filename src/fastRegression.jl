struct BasicReg{L,R} <: RegressionModel
    nobs::Int
    formula::FormulaTerm{L,R}
    coef::Vector{Float64}
    coefnames::Vector{String}
    yname::String
    tss::Float64
    rss::Float64
    residuals::Union{Vector{Float64}, Nothing}
    function BasicReg(nobs, formula::FormulaTerm{L,R}, coef, coefnames, yname, tss, rss, residuals) where {L,R}
        @assert rss >= 0 "Residual sum of squares must be greater than 0"
        @assert tss >= 0 "Total sum of squares must be greater than 0"
        @assert nobs >= 0 "Observations must be greater than 0"
        @assert length(coef) == length(coefnames) "Number of coefficients must be same as number of coefficient names"
        new{L,R}(nobs, formula, coef, coefnames, yname, tss, rss, residuals)
    end
    BasicReg(x::Int, f::FormulaTerm{L,R}) where {L, R} = new{L,R}(x, f)
end

function create_pred_matrix(data::TimelineTable{false}, sch)
    if data.regression_cache === nothing || data.regression_cache.terms != sch.rhs
        temp_data = data.parent[:, internal_termvars(sch.rhs)]
        mat = modelcols(sch.rhs, temp_data)
        data.parent.regression_cache = RegressionCache(
            sch.rhs,
            mat,
            temp_data.dates,
            nothing,
            data.calendar
        )
    end
end


function quick_reg(
    data::TimelineTable{false},
    f::FormulaTerm;
    minobs::Real=0.8,
    save_residuals::Bool=false,
    save_prediction_matrix::Bool=true
)

    if minobs < 1
        minobs = bdayscount(data.calendar, data.dates.left, data.dates.right) * minobs
    end

    if !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
        f = FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs)
    end
    select!(data, internal_termvars(f))

    if length(data) < minobs
        #throw("Too few observations")
        return BasicReg(length(data), f)
    end
    # a Schema is normally built by running schema(f, data)
    # but doing that repeatedly is quite slow and, in this case, does not
    # provide any extra use since all of the columns are already known to be continuous
    # and the other values (mean, var, min, max) are not used later
    sch = apply_schema(
        f,
        StatsModels.Schema(
            Dict(
                term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0))
            ),
    )

    if save_prediction_matrix
        create_pred_matrix(data, sch)
    end

    resp, pred = modelcols(sch, data)
    #resp = modelcols(sch.lhs, data)
    #pred = modelcols(sch.rhs, data; save_matrix=save_prediction_matrix)

    if length(resp) < minobs
        return BasicReg(0, f)
    end

    coef = cholesky!(Symmetric(pred' * pred)) \ (pred' * resp)
    resid = resp .- pred * coef
    rss = sum(abs2, resid)
    tss = sum(abs2, (resp .- mean(resp)))
    yname, xnames = coefnames(sch)

    BasicReg(
        length(resp),
        f,
        coef,
        xnames,
        yname,
        tss,
        rss,
        save_residuals ? resid : nothing
    )
end

function StatsBase.fit(
    ::Type{BasicReg},
    formula::FormulaTerm,
    data::TimelineTable{false}
)
    quick_reg(data, formula)
end




"""

    cache_reg(
        id::Int,
        est_min::Date,
        est_max::Date;
        cols_market::Union{Nothing, Vector{String}}=nothing,
        col_firm::String="ret",
        minobs::Real=.8,
        calendar="CrspMarketCalendar"
    )

An intentionally minamilistic linear regression of a vector of firm data on a matrix
of market data, where the firm data is provided by an Integer id and the range
of dates. Designed to quickly estimate Fama French predicted returns.

## Arguments

- `id::Int`: The firm identifier
- `est_min::Date`: The start of the estimation period
- `est_max::Date`: The end of the estimation period
- `cols_market::Union{Nothing, Vector{String}}=nothing`: A vector of columns that are stored
    in the MARKET_DATA_CACHE, should not be repeated and it is recommended to include `"intercept"`
    as the intercept (typically first). If `nothing` (default) then all columns saved will be used.
    Note that if all columns are used there could be errors in the regression (if the risk free rate
    is stored) since over short periods the risk free rate is often constant and conflicts with
    the intercept
- `col_firm::String="ret"`: The column for the firm vector, typically this is the return, it must be
    in the FIRM_DATA_CACHE data
- `minobs::Real=.8`: Minimum number of observations to run the regression, if the number provided
    is less than 1, then it is assumed to be a ratio (i.e., minimum observations is number of
    businessdays times minobs)
- `calendar="CrspMarketCalendar"`: calendar to use if minobs is less than 1, should be initialized cache if
    used.
- `save_residuals::Bool=false`: Whether or not to save the residuals in the regression
"""

StatsBase.predict(mod::BasicReg, x) = x * coef(mod)

StatsBase.coef(x::BasicReg) = isdefined(x, :coef) ? x.coef : missing
StatsBase.coefnames(x::BasicReg) = isdefined(x, :coefnames) ? x.coefnames : missing
StatsBase.responsename(x::BasicReg) = isdefined(x, :yname) ? x.yname : missing
StatsBase.nobs(x::BasicReg) = x.nobs
StatsBase.dof_residual(x::BasicReg) = isdefined(x, :coef) ? nobs(x) - length(coef(x)) : missing
StatsBase.r2(x::BasicReg) = 1 - (rss(x) / deviance(x))
StatsBase.adjr2(x::BasicReg) = 1 - rss(x) / deviance(x) * (nobs(x) - 1) / dof_residual(x)
StatsBase.islinear(x::BasicReg) = true
StatsBase.deviance(x::BasicReg) = isdefined(x, :tss) ? x.tss : missing
StatsBase.rss(x::BasicReg) = isdefined(x, :rss) ? x.rss : missing
function StatsBase.residuals(x::BasicReg)
    if x.residuals === nothing
        @error("To access residuals, run `cache_reg` with the option `save_residuals=true`")
    else
        x.residuals
    end
end

"""
    predict(rr::RegressionModel, date_start::Date, date_end::Date)

Uses a provided RegressionModel model and the cached saved market data to predict returns
between the two dates provided.

This will work with any RegressionModel provided as long as it provides methods for `coefnames`
that correspond to names stored in the MARKET_DATA_CACHE and a `coef` method. The model should be linear.
"""
# this might be too generalized...
function StatsBase.predict(rr::RegressionModel, data::TimelineTable{false})
    f = rr.formula
    select!(data, internal_termvars(f))
    

    # a Schema is normally built by running schema(f, data)
    # but doing that repeatedly is quite slow and, in this case, does not
    # provide any extra use since all of the columns are already known to be continuous
    # and the other values (mean, var, min, max) are not used later
    sch = apply_schema(
        f,
        StatsModels.Schema(
            Dict(
                term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0))
            ),
    )

    pred = modelcols(sch.rhs, data)
    predict(rr, pred)
end

