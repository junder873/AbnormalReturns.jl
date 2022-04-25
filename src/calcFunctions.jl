
Statistics.var(rr::RegressionModel) = rss(rr) / dof_residual(rr)
Statistics.std(rr::RegressionModel) = sqrt(var(rr))




function bh_return(vals::AbstractVector{Float64})
    out = 1.0
    @simd for x in vals
        out *= (1 + x)
    end
    out - 1
end
function bh_return(vals::AbstractVector{Union{Missing, Float64}})
    out = 1.0
    @simd for x in vals
        if ismissing(x)
            out *= 1
        else
            out *= (1+x)
        end
    end
    out - 1
end

function bh_return(pred::FixedWidthMatrix{N}, coef::SVector{N}) where {N}
    #@assert size(pred, 2) == length(coef) "Got Matrix of size $(size(pred)) and coefficients of $coef $pred"
    out = 1.0
    @simd for i in 1:size(pred)[1]
        out *= (point_pred(pred, coef, i) + 1)
    end
    out - 1
end

bhar(resp, pred, coef) = bh_return(resp) - bh_return(pred, coef)

bhar(x::AbstractVector, y::AbstractVector) = bh_return(x) - bh_return(y)

function cumulative_return(pred::FixedWidthMatrix{N}, coef::SVector{N}) where {N}
    #@assert size(pred, 2) == length(coef) "Got Matrix of size $(size(pred)) and coefficients of $coef $pred"
    out = 0.0
    @simd for i in 1:size(pred)[1]
        out += point_pred(pred, coef, i)
    end
    out
end
car(resp, pred, coef) = sum(resp) - cumulative_return(pred, coef)
car(x::AbstractVector, y::AbstractVector) = sum(skipmissing(x)) - sum(skipmissing(y))

var_diff(x, y) = var(x) + var(y) - 2 * cov(x, y)
# for firm data
"""

Calculates the buy and hold returns (also called geometric return) for TimelineData. If an Integer
is passed, then it is calculated based on the FIRM_DATA_CACHE (for the integer provided), otherwise
is calculated for the MARKET_DATA_CACHE.

These functions treat missing returns in the period implicitly as a zero return.
"""
function bh_return(data::TimelineTable, col; minobs=0.8)
    col = TimelineColumn(col)
    select!(data, [col])
    dates = dates_min_max(norm_dates(data), data_dates(data))
    if nnz(data_missing_bdays(data)) == 0
        x = view(data[col], dates)
    else
        new_missings = get_missing_bdays(data, dates)
        x = view(data[col], dates, new_missings)
    end
    if length(x) < adjust_minobs(minobs, calendar(data), norm_dates(data))
        missing
    else
        bh_return(x)
    end
end

function fill_vector_bh_return!(
    out_vector::Vector{Union{Missing, Float64}},
    iter_data::IterateTimelineTable,
    col::TimelineColumn;
    minobs=0.8
)
    @assert validate_iterator(iter_data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"

    Threads.@threads for (u_id, iter_indexes) in iter_data
        data = parent(iter_data)[u_id, [col], AllowMissing{true}]
        resp = data[col]
        for iter in iter_indexes
            cur_dates = dates_min_max(data_dates(data), iter_dates(iter))
            if nnz(data_missing_bdays(data)) == 0
                x = view(resp, cur_dates)
            else
                new_mssngs = get_missing_bdays(data, cur_dates)
                x = view(resp, cur_dates, new_mssngs)
            end
            
            length(x) < adjust_minobs(minobs, calendar(parent(iter_data)), iter_dates(iter)) && continue

            @inbounds out_vector[iter_index(iter)] = bh_return(x)

        end
    end
    out_vector
end

function bh_return(
    data::IterateTimelineTable,
    col;
    minobs=0.8
)
    col = TimelineColumn(col)
    out = Vector{Union{Missing, Float64}}(missing, total_length(data))
    fill_vector_bh_return!(
        out,
        data,
        col;
        minobs
    )
    if any(ismissing.(out))
        out
    else
        disallowmissing(out)
    end
end

function sum_return(data::TimelineTable, col; minobs=0.8)
    col = TimelineColumn(col)
    select!(data, [col])
    dates = dates_min_max(norm_dates(data), data_dates(data))
    if nnz(data_missing_bdays(data)) == 0
        x = view(data[col], dates)
    else
        new_missings = get_missing_bdays(data, dates)
        x = view(data[col], dates, new_missings)
    end
    if length(x) < adjust_minobs(minobs, calendar(data), norm_dates(data))
        missing
    else
        sum(x)
    end
end

function simple_diff(
    data::TimelineTable,
    firm_col,
    mkt_col,
    fun;
    minobs=0.8
)

    select!(data, [firm_col, mkt_col])
    dates = dates_min_max(norm_dates(data), data_dates(data))
    if nnz(data_missing_bdays(data)) == 0
        x = view(data[firm_col], dates)
        y = view(data[mkt_col], dates)
    else
        new_missings = get_missing_bdays(data, dates)
        x = view(data[firm_col], dates, new_missings)
        y = view(data[mkt_col], dates, new_missings)
    end
    if length(x) < adjust_minobs(minobs, calendar(data), norm_dates(data))
        missing
    else
        fun(x, y)
    end
end

function pred_diff(
    data::TimelineTable,
    rr::RegressionModel,
    fun;
    minobs::Float64=0.8
)
    if !isdefined(rr, :coef)
        return missing
    end
    f = rr.formula
    sch = apply_schema(f, schema(f, data))
    select!(data, internal_termvars(sch))
    if length(data) < adjust_minobs(minobs, calendar(data), norm_dates(data))
        return missing
    end

    resp, pred = modelcols(sch, data)
    fun(resp, pred, coef(rr))
end

function fill_vector_simple!(
    out_vector::Vector{Union{Missing, Float64}},
    iter_data::IterateTimelineTable,
    cache::DataVector,
    col::TimelineColumn,
    fun;
    minobs=0.8
)
    @assert validate_iterator(iter_data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"

    Threads.@threads for (u_id, iter_indexes) in iter_data
        data = parent(iter_data)[u_id, [col], AllowMissing{true}]
        resp = data[col]
        for iter in iter_indexes
            cur_dates = dates_min_max(data_dates(data), iter_dates(iter))
            if nnz(data_missing_bdays(data)) == 0
                x = view(resp, cur_dates)
                y = view(cache, cur_dates)
            else
                new_mssngs = get_missing_bdays(data, cur_dates)
                x = view(resp, cur_dates, new_mssngs)
                y = view(cache, cur_dates, new_mssngs)
            end
            
            length(x) < adjust_minobs(minobs, calendar(parent(iter_data)), iter_dates(iter)) && continue

            @inbounds out_vector[iter_index(iter)] = fun(
                x,
                y,
            )

        end
    end
    out_vector
end

function vector_simple_diff(
    data::IterateTimelineTable,
    col1::TimelineColumn,
    col2::TimelineColumn,
    fun;
    minobs=0.8
)
    cache = parent(data)[first(data)[1], [col2], AllowMissing{true}][col2]

    out = Vector{Union{Missing, Float64}}(missing, total_length(data))
    fill_vector_simple!(
        out,
        data,
        cache,
        col1,
        fun;
        minobs
    )
    if any(ismissing.(out))
        out
    else
        disallowmissing(out)
    end
end


function fill_vector_pred!(
    out_vector::Vector{Union{Missing, Float64}},
    iter_data::IterateTimelineTable{T},
    cache::NTuple{N, DataVector},
    sch,
    rrs::Vector{BasicReg{L, R, N}},
    fun;
    minobs=0.8
) where {T, L, R, N}
    @assert validate_iterator(iter_data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"
    @assert length(out_vector) == length(rrs) "Length of out_vector is not the same as number of regressions"
    cols = internal_termvars(sch)
    Threads.@threads for (u_id, iter_indexes) in iter_data
        data = parent(iter_data)[u_id, cols, AllowMissing{true}]
        resp = datavector_modelcols(int_lhs(sch), data)
        for iter in iter_indexes
            isdefined(rrs[iter_index(iter)], :coef) || continue
            cur_dates = dates_min_max(data_dates(data), iter_dates(iter))
            if nnz(data_missing_bdays(data)) == 0
                x = view(resp, cur_dates)
                y = FixedWidthMatrix(cache, cur_dates)
            else
                new_mssngs = get_missing_bdays(data, cur_dates)
                x = view(resp, cur_dates, new_mssngs)
                y = FixedWidthMatrix(cache, cur_dates, new_mssngs)
            end
            
            length(x) < adjust_minobs(minobs, calendar(parent(iter_data)), iter_dates(iter)) && continue

            @inbounds out_vector[iter_index(iter)] = fun(
                x,
                y,
                coef(rrs[iter_index(iter)])
            )

        end
    end
    out_vector
end

function vector_pred_diff(
    data::IterateTimelineTable,
    rrs::AbstractVector{BasicReg{L, R, N}},
    fun;
    minobs=0.8
) where {L, R, N}
    f = rrs[1].formula
    sch = apply_schema(f, schema(f, parent(data)))
    cache = create_pred_matrix(data, sch)

    out = Vector{Union{Missing, Float64}}(missing, total_length(data))
    fill_vector_pred!(
        out,
        data,
        cache,
        sch,
        rrs,
        fun;
        minobs
    )
    if any(ismissing.(out))
        out
    else
        disallowmissing(out)
    end
end


"""
    bhar(
        data::TimelineTable{Mssng, T, MNames, FNames},
        firm_col=FNames[1],
        mkt_col=MNames[1];
        minobs=0.8
    ) where {Mssng, T, MNames, FNames}

    bhar(
        data::TimelineTable,
        rr::RegressionModel;
        minobs=0.8
    )

    bhar(
        data::IterateTimelineTable{T, MNames, FNames},
        firm_col=FNames[1],
        mkt_col=MNames[1];
        minobs=0.8
    ) where {T, MNames, FNames}

    bhar(
        data::IterateTimelineTable,
        rrs::AbstractVector{<:BasicReg};
        minobs=0.8
    )

Calculates the difference between buy and hold returns (also referred to as geometric returns) for a firm and a benchmark.
If a regression is passed, then the benchmark is based on the coefficients from that regression and the performance of the benchmarks
in the regression. These are sometimes called Fama-French abnormal returns. Simple abnormal returns use a market index as the benchmark
(such as the S&P 500 or a value weighted return of all firms).

Similar to constructing the regression, passing an `IterateTimelineTable` will return a Vector and uses a more optimized method.
"""
function bhar(
    data::TimelineTable{Mssng, T, MNames, FNames},
    firm_col=FNames[1],
    mkt_col=MNames[1];
    minobs=0.8
) where {Mssng, T, MNames, FNames}
    simple_diff(data, TimelineColumn(firm_col), TimelineColumn(mkt_col), bhar; minobs)
end

function bhar(
    data::TimelineTable,
    rr::RegressionModel;
    minobs=0.8
)
    pred_diff(data, rr, bhar; minobs)
end

function bhar(
    data::IterateTimelineTable{T, MNames, FNames},
    firm_col=FNames[1],
    mkt_col=MNames[1];
    minobs=0.8
) where {T, MNames, FNames}
    vector_simple_diff(data, TimelineColumn(firm_col), TimelineColumn(mkt_col), bhar; minobs)
end

function bhar(
    data::IterateTimelineTable,
    rrs::AbstractVector{<:BasicReg};
    minobs=0.8
)
    vector_pred_diff(data, rrs, bhar; minobs)
end

"""
    car(
        data::TimelineTable{Mssng, T, MNames, FNames},
        firm_col=FNames[1],
        mkt_col=MNames[1];
        minobs=0.8
    ) where {Mssng, T, MNames, FNames}

    car(
        data::TimelineTable,
        rr::RegressionModel;
        minobs=0.8
    )

    car(
        data::IterateTimelineTable{T, MNames, FNames},
        firm_col=FNames[1],
        mkt_col=MNames[1];
        minobs=0.8
    ) where {T, MNames, FNames}

    car(
        data::IterateTimelineTable,
        rrs::AbstractVector{<:BasicReg};
        minobs=0.8
    )

Calculates the cumulative returns of a firm over a benchmark (through addition of each return).
If a regression is passed, then the benchmark is based on the coefficients from that regression and the performance of the benchmarks
in the regression. These are sometimes called Fama-French abnormal returns. Simple abnormal returns use a market index as the benchmark
(such as the S&P 500 or a value weighted return of all firms).

Similar to constructing the regression, passing an `IterateTimelineTable` will return a Vector and uses a more optimized method.
"""
function car(
    data::TimelineTable{Mssng, T, MNames, FNames},
    firm_col=FNames[1],
    mkt_col=MNames[1];
    minobs=0.8
) where {Mssng, T, MNames, FNames}
    simple_diff(data, TimelineColumn(firm_col), TimelineColumn(mkt_col), car; minobs)
end

function car(
    data::TimelineTable,
    rr::RegressionModel;
    minobs=0.8
)
    pred_diff(data, rr, car; minobs)
end

function car(
    data::IterateTimelineTable{T, MNames, FNames},
    firm_col=FNames[1],
    mkt_col=MNames[1];
    minobs=0.8
) where {T, MNames, FNames}
    vector_simple_diff(data, TimelineColumn(firm_col), TimelineColumn(mkt_col), car; minobs)
end

function car(
    data::IterateTimelineTable{T},
    rrs::AbstractVector{<:BasicReg};
    minobs=0.8
) where {T}
    vector_pred_diff(data, rrs, car; minobs)
end

function Statistics.var(
    data::TimelineTable{Mssng, T, MNames, FNames},
    firm_col=FNames[1],
    mkt_col=MNames[1];
    minobs=0.8
) where {Mssng, T, MNames, FNames}
    simple_diff(data, TimelineColumn(firm_col), TimelineColumn(mkt_col), var_diff; minobs)
end

function Statistics.std(
    data::TimelineTable{Mssng, T, MNames, FNames},
    firm_col=FNames[1],
    mkt_col=MNames[1];
    minobs=0.8
) where {Mssng, T, MNames, FNames}
    sqrt(var(data, col_firm, col_market; minobs))
end

function get_coefficient_pos(rr::RegressionModel, coefname::String...)
    for x in coefname
        if x ∈ coefnames(rr)
            return findfirst(x .== coefnames(rr))
        end
    end
    @error("None of $(coefname) is in the RegressionModel model.")
end

function get_coefficient_val(rr::RegressionModel, coefname::String...)
    ismissing(coefnames(rr)) && return missing
    coef(rr)[get_coefficient_pos(rr, coefname...)]
end

# as an optimization, if all are the same regression, then just find the coefname once
function get_coefficient_val(rrs::Vector{BasicReg{L, R, N}}, coefname::String...) where {L, R, N}
    out = Vector{Union{Missing, Float64}}(missing, length(rrs))
    pos = 0
    for (i, rr) in enumerate(rrs)
        if pos == 0 && !ismissing(coefnames(rr))
            pos = get_coefficient_pos(rr, coefname...)
        end
        if pos != 0
            out[i] = coef(rr)[pos]
        end
    end
    out
end
"""
    alpha(rr::RegressionModel, coefname::String...="intercept")

"alpha" in respect to the the CAPM model, i.e., the intercept in the model.
This is the alpha from the estimation period.

This function finds the position of the coefficient name provided, defaults to "intercept".
If the coefname is not in the regression, then this function returns an error.
"""
alpha(rr::RegressionModel, coefname::String...="(Intercept)") = get_coefficient_val(rr, coefname...)
alpha(rrs::Vector{BasicReg{L, R, N}}, coefname::String...="(Intercept)") where {L, R, N} = get_coefficient_val(rrs, coefname...)


"""
    beta(rr::RegressionModel, coefname::String...=["mkt", "mktrf", "vwretd", "ewretd"])

"beta" in respect to the CAPM model, i.e., the coefficient on the market return minus the risk free rate.
This is the beta from the estimation period.

This function finds the position of the coefficient name provided, defaults to several common market returns.
If the coefname is not in the regression, then this function returns an error.
"""
beta(rr::RegressionModel, coefname::String...=("mkt", "mktrf", "vwretd", "ewretd")...) = get_coefficient_val(rr, coefname...)
beta(rrs::Vector{BasicReg{L, R, N}}, coefname::String...=("mkt", "mktrf", "vwretd", "ewretd")...) where {L, R, N} = get_coefficient_val(rrs, coefname...)