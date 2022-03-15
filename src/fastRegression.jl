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
        # @assert rss >= 0 "Residual sum of squares must be greater than 0"
        # @assert tss >= 0 "Total sum of squares must be greater than 0"
        # @assert nobs >= 0 "Observations must be greater than 0"
        # @assert length(coef) == length(coefnames) "Number of coefficients must be same as number of coefficient names"
        new{L,R}(nobs, formula, coef, coefnames, yname, tss, rss, residuals)
    end
    BasicReg(x::Int, f::FormulaTerm{L,R}) where {L, R} = new{L,R}(x, f)
end

function BasicReg(
    resp::AbstractVector{Float64},
    pred::AbstractMatrix{Float64},
    yname::String,
    xnames::Vector{String},
    f::FormulaTerm{L,R};
    save_residuals=false
)::BasicReg{L,R} where {L,R}
    if length(resp) <= size(pred, 2)
        return BasicReg(length(resp), f)
    end
    # display(pred)
    # display(resp)

    coef = cholesky!(Symmetric(pred' * pred)) \ (pred' * resp)
    resid = resp .- pred * coef
    rss = sum(abs2, resid)
    tss = sum(abs2, (resp .- mean(resp)))

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

function BasicReg(
    data::TimelineTable,
    cache::RegressionCache,
    sch::FormulaTerm{L2,R2},
    f::FormulaTerm{L,R};
    minobs::V=0.8,
    save_residuals::Bool=false
)::BasicReg{L,R} where {L,R,L2,R2,V<:Real}

    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    if length(data) < minobs
        return BasicReg(length(data), f)
    end
    resp = modelcols(int_lhs(sch), data)
    # resp = view(data[sch.lhs.sym], data.dates)
    # pred = view(cache, data.dates)
    pred = cache[data.dates, data.missing_bdays]
    # BasicReg(
    #     resp,
    #     pred,
    #     sch,
    #     f;
    #     save_residuals
    # )
    yname, xnames = coefnames(sch)
    BasicReg(
        resp,
        pred,
        yname,
        xnames,
        f;
        save_residuals
    )
end

function BasicReg(
    data::TimelineTable,
    sch::FormulaTerm{L2,R2},
    f::FormulaTerm{L,R};
    minobs::V=0.8,
    save_residuals::Bool=false
)::BasicReg{L,R,L2,R2} where {L,R,L2,R2,V<:Real}
    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    if length(data) < minobs
        return BasicReg(length(data), f)
    end
    resp, pred = modelcols(sch, data)
    BasicReg(
        resp,
        pred,
        sch,
        f;
        save_residuals
    )
end

function create_pred_matrix(data::TimelineTable, sch)
    # This allows missing to make sure all dates are in the data
    # it is not necessary to include missing bdays since the
    # timelinetable already has the same information

    #temp_data = allowmissing(data[:, internal_termvars(sch.rhs)])
    mat = modelcols(sch.rhs, data)
    #mat = coalesce.(mat, 0.0)
    #update_request!(data, cols=old_cols, dates=old_dates)
    RegressionCache(
            mat,
            data.dates,
            data.calendar
        )
end


function quick_reg(
    data::TimelineTable,
    f::FormulaTerm{L,R};
    minobs::V=0.8,
    save_residuals::Bool=false
) where {L,R,V<:Real}

    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
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
    sch = apply_schema(f, schema(f, data))
    resp, pred = modelcols(sch, data)
    yname, xnames = coefnames(sch)
    BasicReg(
        resp,
        pred,
        yname,
        xnames,
        f;
        save_residuals
    )
end

function fill_vector_reg(
    data::IterateMarketData{T},
    out_vector::Vector{BasicReg{L,R}},
    cache::RegressionCache,
    sch::FormulaTerm{L2,R2},
    cols::Vector{TimelineColumn},
    f::FormulaTerm{L,R};
    minobs::V=0.8,
    save_residuals::Bool=false
) where {T, L, R, L2, R2, V<:Real}
    @assert validate_iterator(data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"
    Threads.@threads for (u_id, iters, local_data) in data
        select!(local_data, cols)
        @inbounds for (i, d1, d2) in iters
            update_dates!(local_data, d1 .. d2)
            out_vector[i] = BasicReg(
                local_data,
                cache,
                sch,
                f;
                minobs,
                save_residuals
            )
        end
    end
    out_vector
end

function vector_reg(
    parent_data::MarketData{T},
    ids::AbstractVector{T},
    date_mins::AbstractVector{Date},
    date_maxs::AbstractVector{Date},
    f::FormulaTerm{L,R};
    minobs::V=0.8,
    save_residuals::Bool=false
) where {T,L,R,V<:Real}
    if !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
        f = FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs)
    end
    cols = internal_termvars(f)
    sch = apply_schema(f, schema(f, parent_data))
    cache = create_pred_matrix(parent_data[ids[1]], sch)
    out = fill(BasicReg(0, f), length(ids))
    fill_vector_reg(
        IterateMarketData(parent_data, ids, date_mins, date_maxs),
        out,
        cache,
        sch,
        cols,
        f;
        minobs,
        save_residuals
    )
end


"""

    cache_reg(
        id::Int,
        est_min::Date,
        est_max::Date;
        cols_market::Union{Nothing, Vector{String}}=nothing,
        col_firm::String="ret",
        minobs=.8,
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
- `minobs=.8`: Minimum number of observations to run the regression, if the number provided
    is less than 1, then it is assumed to be a ratio (i.e., minimum observations is number of
    businessdays times minobs)
- `calendar="CrspMarketCalendar"`: calendar to use if minobs is less than 1, should be initialized cache if
    used.
- `save_residuals::Bool=false`: Whether or not to save the residuals in the regression
"""

StatsBase.predict(mod::BasicReg, x) = x * coef(mod)

StatsBase.coef(x::BasicReg) = isdefined(x, :coef) ? x.coef : missing
StatsBase.coefnames(x::BasicReg) = isdefined(x, :coef) ? x.coefnames : missing
StatsBase.responsename(x::BasicReg) = isdefined(x, :coef) ? x.yname : missing
StatsBase.nobs(x::BasicReg) = x.nobs
StatsBase.dof_residual(x::BasicReg) = isdefined(x, :coef) ? nobs(x) - length(coef(x)) : missing
StatsBase.r2(x::BasicReg) = 1 - (rss(x) / deviance(x))
StatsBase.adjr2(x::BasicReg) = 1 - rss(x) / deviance(x) * (nobs(x) - 1) / dof_residual(x)
StatsBase.islinear(x::BasicReg) = true
StatsBase.deviance(x::BasicReg) = isdefined(x, :coef) ? x.tss : missing
StatsBase.rss(x::BasicReg) = isdefined(x, :coef) ? x.rss : missing
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
# # this might be too generalized...
# function StatsBase.predict(rr::RegressionModel, data::TimelineTable{false})
#     f = rr.formula
#     select!(data, internal_termvars(f))
    

#     # a Schema is normally built by running schema(f, data)
#     # but doing that repeatedly is quite slow and, in this case, does not
#     # provide any extra use since all of the columns are already known to be continuous
#     # and the other values (mean, var, min, max) are not used later
#     sch = apply_schema(f, schema(f, data))

#     pred = modelcols(sch.rhs, data)
#     predict(rr, pred)
# end

function rhs_str(nms, vals; intercept = "(Intercept)")
    out = String[]
    for (nm, val) in zip(nms, vals)
        val_str = string(round(val, digits=3))
        if nm == intercept
            push!(out, val_str)
        else
            push!(out, val_str * "*" * nm)
        end
    end
    join(out, " + ")
end

function Base.show(io::IO, rr::BasicReg)
    if !isdefined(rr, :coef)
        print(io, "Obs: $(nobs(rr)), $(rr.formula)")
    else
        print(io, "Obs: $(nobs(rr)), $(responsename(rr)) ~ $(rhs_str(coefnames(rr), coef(rr))), AdjR2: ", round(adjr2(rr) * 100, digits=3), "%")
        # mat = hcat(
        #     coefnames(rr),
        #     string.(round.(coef(rr), digits=3))
        # )
        # PrettyTables.pretty_table(io, mat)
    end
end