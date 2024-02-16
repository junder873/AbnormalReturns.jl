




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
    out = one(T)
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

Calculates the buy and hold returns (also called geometric return).

These functions treat missing returns in the period implicitly as a zero return.
"""
function bh_return(data::FixedTable{1}; minobs=0.8)
    d = data[:, 1]
    if length(d) < adjust_minobs(minobs, data)
        missing
    else
        bh_return(d)
    end
end

function bh_return(
    data::IterateFixedTable{T, 1};
    minobs=0.8
) where {T}
    out = Vector{Union{Missing, Float64}}(missing, length(data))
    Threads.@threads for i in 1:length(data)
        out[i] = bh_return(data[i]; minobs)
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
    minobs=0.8
) where {N1}
    #@assert N1 == N2 + 1 "Dimensions are mismatched"
    if !isdefined(rr, :coefnames)
        return missing
    end
    if size(data, 1) < adjust_minobs(minobs, data)
        return missing
    else
        fun(data[:, 1], pred_matrix(data), coef(rr))
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

function pred_diff(
    data::Tuple,
    rr,
    args...;
    vargs...
)
    f = if isa(rr, AbstractVector)
        first(rr).formula
    else
        rr.formula
    end
    pred_diff(data[1][data[2:end]..., f, check_intercept=false], rr, args...; vargs...)
end



"""
    bhar(
        data::FixedTable{T, 2};
        minobs=0.8
    ) where {T}

    bhar(
        data::FixedTable,
        rr::RegressionModel;
        minobs=0.8
    )

    bhar(
        data::IterateFixedTable{T, 2};
        minobs=0.8
    ) where {T, MNames, FNames}

    bhar(
        data::IterateFixedTable,
        rrs::AbstractVector{<:BasicReg};
        minobs=0.8
    )

Calculates the difference between buy and hold returns (also referred to as geometric returns) for a firm and a benchmark.
If a regression is passed, then the benchmark is based on the coefficients from that regression and the performance of the benchmarks
in the regression. These are sometimes called Fama-French abnormal returns. If no regression is passed,
abnormal returns are calculated as the difference between the first and second columns in
the FixedTable (second column is typically the benchmark such as the S&P 500 or a value weighted return of all firms).

Similar to constructing the regression, passing an `IterateFixedTable` will return a Vector and uses a more optimized method.
"""
bhar(data::Union{IterateFixedTable, FixedTable}; minobs=0.8) = simple_diff(data, bhar; minobs)
bhar(data::Union{IterateFixedTable, FixedTable, Tuple}, rr; minobs=0.8) = pred_diff(data, rr, bhar; minobs)

"""
    car(
        data::FixedTable{T, 2};
        minobs=0.8
    ) where {T}

    car(
        data::FixedTable,
        rr::RegressionModel;
        minobs=0.8
    )

    car(
        data::IterateFixedTable{T, 2};
        minobs=0.8
    ) where {T, MNames, FNames}

    car(
        data::IterateFixedTable,
        rrs::AbstractVector{<:BasicReg};
        minobs=0.8
    )

Calculates the cumulative returns of a firm over a benchmark (through addition of each return).
If a regression is passed, then the benchmark is based on the coefficients from that regression and the performance of the benchmarks
in the regression. These are sometimes called Fama-French abnormal returns. If no regression is passed,
abnormal returns are calculated as the difference between the first and second columns in
the FixedTable (second column is typically the benchmark such as the S&P 500 or a value weighted return of all firms).

Similar to constructing the regression, passing an `IterateFixedTable` will return a Vector and uses a more optimized method.
"""
car(data::Union{IterateFixedTable, FixedTable}; minobs=0.8) = simple_diff(data, car; minobs)
car(data::Union{IterateFixedTable, FixedTable, Tuple}, rr; minobs=0.8) = pred_diff(data, rr, car; minobs)


"""
    var[std](rr::Union{AbstractVector{<:RegressionModel}, RegressionModel})
    var[std](data::Union{IterateFixedTable, FixedTable}; minobs=0.8)

If a regression model is passed, then this calculates the variance (standard deviation)
based on the residual sum of squares divided by the degrees of freedom. A vector of
RegressionModel will return the same length of vector results.

If a FixedTable is passed (or an IterateFixedTable), and that contains only one column,
then the variance (standard deviation) is calculated for that column. If it has two
columns, then the calculation is based on the difference between the columns.
"""
Statistics.var(data::Union{IterateFixedTable{T, 2}, FixedTable{2}}; minobs=0.8) where {T} = simple_diff(data, var_diff; minobs)
function Statistics.var(data::FixedTable{1}; minobs=0.8)
    d = data[:, 1]
    if length(d) < adjust_minobs(minobs, data)
        missing
    else
        var(d)
    end
end
Statistics.var(data::IterateFixedTable{T, 1}; minobs=0.8) where {T} = var.(data; minobs)

"""
See [`var`](@ref).
"""
Statistics.std(data::FixedTable; minobs=0.8) = sqrt(var(data; minobs))
Statistics.std(data::IterateFixedTable; minobs=0.8) = sqrt.(var(data; minobs))

Statistics.var(rr::RegressionModel) = rss(rr) / dof_residual(rr)
Statistics.std(rr::RegressionModel) = sqrt(var(rr))
Statistics.var(rrs::AbstractVector{<:RegressionModel}) = var.(rrs)
Statistics.std(rrs::AbstractVector{<:RegressionModel}) = sqrt.(var(rrs))


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
function get_coefficient_val(rrs::Vector{<:RegressionModel}, coefname::String...)
    out = Vector{Union{Missing, Float64}}(missing, length(rrs))
    pos = 0
    for (i, rr) in enumerate(rrs)
        if pos == 0 && !ismissing(coefnames(rr))
            pos = get_coefficient_pos(rr, coefname...)
        end
        if pos != 0 && !ismissing(coefnames(rr))
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
alpha(rrs::Vector{<:RegressionModel}, coefname::String...="(Intercept)") = get_coefficient_val(rrs, coefname...)


"""
    beta(rr::RegressionModel, coefname::String...=["mkt", "mktrf", "vwretd", "ewretd"])

"beta" in respect to the CAPM model, i.e., the coefficient on the market return minus the risk free rate.
This is the beta from the estimation period.

This function finds the position of the coefficient name provided, defaults to several common market returns.
If the coefname is not in the regression, then this function returns an error.
"""
beta(rr::RegressionModel, coefname::String...=("mkt", "mktrf", "vwretd", "ewretd")...) = get_coefficient_val(rr, coefname...)
beta(rrs::Vector{<:RegressionModel}, coefname::String...=("mkt", "mktrf", "vwretd", "ewretd")...) = get_coefficient_val(rrs, coefname...)