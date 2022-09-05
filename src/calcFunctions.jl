
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

function bh_return(pred::FixedTable{N, T}, coef::SVector{N, T}) where {N, T}
    #@assert size(pred, 2) == length(coef) "Got Matrix of size $(size(pred)) and coefficients of $coef $pred"
    out = zero(T)
    @simd for i in 1:size(pred)[1]
        out *= (point_pred(pred, coef, i) + 1)
    end
    out - 1
end

bhar(resp, pred, coef) = bh_return(resp) - bh_return(pred, coef)

bhar(x::AbstractVector, y::AbstractVector) = bh_return(x) - bh_return(y)

function cumulative_return(pred::FixedTable{N, T}, coef::SVector{N, T}) where {N, T}
    #@assert size(pred, 2) == length(coef) "Got Matrix of size $(size(pred)) and coefficients of $coef $pred"
    out = zero(T)
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
function bh_return(data::FixedTable{1}; minobs=0.8)
    d = data[:, 1]
    if length(d) < adjust_minobs(minobs, data)
        missing
    else
        bh_return(data[:, 1])
    end
end

function bh_return(
    data::IterateFixedTable{1};
    minobs=0.8
)
    out = Vector{Union{Missing, Float64}}(missing, length(data))
    Threads.@threads for i in 1:length(data)
        out[i] = bh_return(data[i], minobs)
    end
    if any(ismissing.(out))
        out
    else
        disallowmissing(out)
    end
end

function sum_return(data::FixedTable{1}; minobs=0.8)
    d = data[:, 1]
    if length(d) < adjust_minobs(minobs, data)
        missing
    else
        sum(d)
    end
end

function simple_diff(
    data::FixedTable{2},
    fun;
    minobs=0.8
)
    x = data[:, 1]
    y = data[:, 2]
    if length(x) < adjust_minobs(minobs, data)
        missing
    else
        fun(x, y)
    end
end

function pred_diff(
    data::FixedTable{N1},
    rr::RegressionModel,
    fun;
    minobs::Float64=0.8
) where {N1}
    #@assert N1 == N2 + 1 "Dimensions are mismatched"
    if !isdefined(rr, :coef)
        return missing
    end
    if size(data, 1) < adjust_minobs(minobs, data)
        return missing
    else
        fun(data[:, 1], resp_matrix(data), coef(rr))
    end
end

function simple_diff(
    data::IterateFixedTable{T, 2},
    fun;
    minobs=0.8
) where {T}
    out = Vector{Union{Missing, Float64}}(missing, length(data))
    Threads.@threads for i in 1:length(data)
        out[i] = simple_diff(data[i], fun; minobs)
    end
    if any(ismissing.(out))
        out
    else
        disallowmissing(out)
    end
end

function pred_diff(
    data::IterateFixedTable{T, N1},
    rrs::AbstractVector{BasicReg{L, R, N2}},
    fun;
    minobs=0.8
) where {L, R, N1, N2, T}
    @assert N1 == N2 + 1 "Dimensions are mismatched"
    @assert length(data) == length(rrs) "Vectors are not the same length"
    out = Vector{Union{Missing, Float64}}(missing, length(data))
    Threads.@threads for i in 1:length(data)
        out[i] = pred_diff(data[i], rrs[i], fun; minobs)
    end
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
        data::IterateFixedTable{T, MNames, FNames},
        firm_col=FNames[1],
        mkt_col=MNames[1];
        minobs=0.8
    ) where {T, MNames, FNames}

    bhar(
        data::IterateFixedTable,
        rrs::AbstractVector{<:BasicReg};
        minobs=0.8
    )

Calculates the difference between buy and hold returns (also referred to as geometric returns) for a firm and a benchmark.
If a regression is passed, then the benchmark is based on the coefficients from that regression and the performance of the benchmarks
in the regression. These are sometimes called Fama-French abnormal returns. Simple abnormal returns use a market index as the benchmark
(such as the S&P 500 or a value weighted return of all firms).

Similar to constructing the regression, passing an `IterateFixedTable` will return a Vector and uses a more optimized method.
"""
bhar(data::Union{IterateFixedTable, FixedTable}; minobs=0.8) = simple_diff(data, bhar; minobs)
bhar(data::Union{IterateFixedTable, FixedTable}, rr; minobs=0.8) = pred_diff(data, rr, bhar; minobs)

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
        data::IterateFixedTable{T, MNames, FNames},
        firm_col=FNames[1],
        mkt_col=MNames[1];
        minobs=0.8
    ) where {T, MNames, FNames}

    car(
        data::IterateFixedTable,
        rrs::AbstractVector{<:BasicReg};
        minobs=0.8
    )

Calculates the cumulative returns of a firm over a benchmark (through addition of each return).
If a regression is passed, then the benchmark is based on the coefficients from that regression and the performance of the benchmarks
in the regression. These are sometimes called Fama-French abnormal returns. Simple abnormal returns use a market index as the benchmark
(such as the S&P 500 or a value weighted return of all firms).

Similar to constructing the regression, passing an `IterateFixedTable` will return a Vector and uses a more optimized method.
"""
car(data::Union{IterateFixedTable, FixedTable}; minobs=0.8) = simple_diff(data, car; minobs)
car(data::Union{IterateFixedTable, FixedTable}, rr; minobs=0.8) = pred_diff(data, rr, car; minobs)

Statistics.var(data::Union{IterateFixedTable, FixedTable}; minobs=0.8) = simple_diff(data, var_diff; minobs)
Statistics.std(data::Union{IterateFixedTable, FixedTable}; minobs=0.8) = std(var(data; minobs))


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