
Statistics.var(rr::RegressionModel) = rss(rr) / dof_residual(rr)
Statistics.std(rr::RegressionModel) = sqrt(var(rr))

# it seems like I want some kind of formula with a similar minobs option for a lot
# of these functions as well.
function Statistics.var(
    data::TimelineTable,
    col_firm,
    col_market;
    minobs=0.8
)
    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    select!(data, [col_firm, col_market])
    if length(data) < minobs
        return missing
    end
    var(data[:, col_firm] .- data[:, col_market])
end

function Statistics.std(
    data::TimelineTable,
    col_firm,
    col_market;
    minobs=0.8
)
    sqrt(var(data, col_firm, col_market; minobs))
end

function bh_return(vals::AbstractVector{Float64})
    out = 0.0
    @simd for x in vals
        out *= (1 + x)
    end
    out - 1
end
function bh_return(vals::AbstractVector{Union{Missing, Float64}})
    out = 0.0
    @simd for x in vals
        if ismissing(x)
            out *= 1
        else
            out *= (1+x)
        end
    end
    out - 1
end

function bh_return(pred::AbstractMatrix, coef)
    out = 0.0
    @simd for i in 1:size(pred, 2)
        @inbounds out *= (fast_pred(pred, coef, i) + 1)
    end
    out - 1
end

bhar(resp, pred, coef) = bh_return(resp) - bh_return(pred, coef)

bhar(x::AbstractVector, y::AbstractVector) = bh_return(x) - bh_return(y)

function cumulative_return(pred::AbstractMatrix, coef)
    out = 0.0
    @simd for i in 1:size(pred, 2)
        @inbounds out += fast_pred(pred, coef, i)
    end
    out
end
car(resp, pred, coef) = sum(resp) - cumulative_return(pred, coef)
car(x::AbstractVector, y::AbstractVector) = sum(skipmissing(x)) - sum(skipmissing(y))

# for firm data
"""
    bh_return(id::Int, date_start::Date, date_end::Date, col_firm::String="ret")
    bh_return(date_start::Date, date_end::Date, cols_market::String="vwretd")

Calculates the buy and hold returns (also called geometric return) for TimelineData. If an Integer
is passed, then it is calculated based on the FIRM_DATA_CACHE (for the integer provided), otherwise
is calculated for the MARKET_DATA_CACHE.

These functions treat missing returns in the period implicitly as a zero return.
"""
function bh_return(data::TimelineTable, col; minobs=0.8)
    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    t = data[:, col]
    if length(t) < minobs
        missing
    else
        bh_return(data[:, col])
    end
end

function sum_return(data::TimelineTable, col; minobs=0.8)
    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    t = data[:, col]
    if length(t) < minobs
        missing
    else
        sum(data[:, col])
    end
end

function simple_diff(
    data::TimelineTable,
    firm_col,
    mkt_col,
    fun;
    minobs=0.8
)
    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    select!(data, [firm_col, mkt_col])
    fun(data, firm_col; minobs) - fun(data, mkt_col; minobs)
end

function pred_diff(
    data::TimelineTable,
    rr::RegressionModel,
    fun;
    minobs::Float64=0.8
)
    if minobs < 1
        minobs = bdayscount(data.calendar, dt_min(data), dt_max(data)) * minobs
    end
    if !isdefined(rr, :coef)
        return missing
    end
    f = rr.formula
    sch = apply_schema(f, schema(f, data))
    select!(data, internal_termvars(sch))
    #data = dropmissing(data)
    if length(data) < minobs || length(data) <= length(names(data))
        return missing
    end
    # since the values in a continuous term of mean/var/min/max are not used here,
    # this just creates a schema from the available values without
    sch = apply_schema(f, schema(f, data))
    resp, pred = modelcols(sch, data)
    fun(resp, predict(rr, pred))
end

function fill_vector_pred!(
    out_vector::Vector{Union{Missing, Float64}},
    parent_data::MarketData{T},
    iters::Dict{T, Vector{IterateOutput}},
    cache::DataMatrix,
    sch,
    rrs::Vector{BasicReg{L,R}},
    fun
) where {T, L, R}
    @assert validate_iterator(iters, out_vector) "Length of out_vector does not match the number of indexes in the iterator"
    cols = internal_termvars(sch)
    Threads.@threads for u_id in keys(iters) |> collect
        data = parent_data[u_id, cols, AllowMissing{true}]
        resp = datavector_modelcols(int_lhs(sch), data)
        for iter in iters[u_id]
            isdefined(rrs[iter_index(iter)], :coef) || continue
            cur_dates = dates_min_max(data_dates(data), iter_dates(iter))
            if nnz(data_missing_bdays(data)) == 0
                x = view(resp, cur_dates)
                y = view(cache, cur_dates)
            else
                x = view(resp, cur_dates, new_mssngs)
                y = view(cache, cur_dates, new_mssngs)
            end
            
            length(x) < iter_minobs(iter) && continue

            @inbounds out_vector[iter_index(iter)] = fun(
                x,
                y,
                coef(rrs[iter_index(iter)])
                #predict(rrs[iter_index(iter)], y)
            )

        end
    end
    out_vector
end

function vector_pred_diff(
    parent_data::MarketData{T},
    ids::AbstractVector{T},
    date_mins::AbstractVector{Date},
    date_maxs::AbstractVector{Date},
    rrs::AbstractVector{BasicReg{L, R}},
    fun;
    minobs=0.8
) where {T, L, R}
    f = rrs[1].formula
    sch = apply_schema(f, schema(f, parent_data))
    cols = internal_termvars(sch)
    
    cache = create_pred_matrix(
        parent_data[ids[1], internal_termvars(sch.rhs), AllowMissing{true}],
        sch
    )
    out = Vector{Union{Missing, Float64}}(missing, length(ids))
    iter_dict = construct_id_dict(ids, date_mins, date_maxs, calendar(parent_data); minobs)
    fill_vector_pred!(
        out,
        parent_data,
        iter_dict,
        cache,
        sch,
        rrs,
        fun
    )
    if any(ismissing.(out))
        out
    else
        disallowmissing(out)
    end
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
    firm_col,
    mkt_col;
    minobs=0.8
)
    simple_diff(data, firm_col, mkt_col, bh_return; minobs)
end

function bhar(
    data::TimelineTable,
    rr::RegressionModel;
    minobs=0.8
)
    pred_diff(data, rr, bhar; minobs)
end

function bhar(
    parent_data::MarketData{T},
    ids::AbstractVector{T},
    date_mins::AbstractVector{Date},
    date_maxs::AbstractVector{Date},
    rrs::AbstractVector{<:BasicReg};
    minobs=0.8
) where {T}
    vector_pred_diff(parent_data, ids, date_mins, date_maxs, rrs, bhar; minobs)
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
    firm_col,
    mkt_col;
    minobs=0.8
)
    simple_diff(data, firm_col, mkt_col, sum_return; minobs)
end

function car(
    data::TimelineTable,
    rr::RegressionModel;
    minobs=0.8
)
    pred_diff(data, rr, car; minobs)
end

function car(
    parent_data::MarketData{T},
    ids::AbstractVector{T},
    date_mins::AbstractVector{Date},
    date_maxs::AbstractVector{Date},
    rrs::AbstractVector{<:BasicReg};
    minobs=0.8
) where {T}
    vector_pred_diff(parent_data, ids, date_mins, date_maxs, rrs, car; minobs)
end


function get_coefficient_val(rr::RegressionModel, coefname::String...)
    ismissing(coefnames(rr)) && return missing
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
