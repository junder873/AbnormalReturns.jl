
struct BasicReg{L, R, N} <: RegressionModel
    nobs::Int
    formula::FormulaTerm{L,R}
    coef::SVector{N, Float64}
    coefnames::SVector{N, String}
    yname::String
    tss::Float64
    rss::Float64
    residuals::Union{Vector{Float64}, Nothing}
    function BasicReg(nobs, formula::FormulaTerm{L,R}, coef::SVector{N}, coefnames::SVector{N}, yname, tss, rss, residuals) where {L,R,N}
        # @assert rss >= 0 "Residual sum of squares must be greater than 0"
        # @assert tss >= 0 "Total sum of squares must be greater than 0"
        # @assert nobs >= 0 "Observations must be greater than 0"
        # @assert length(coef) == length(coefnames) "Number of coefficients must be same as number of coefficient names"
        new{L,R,N}(nobs, formula, coef, coefnames, yname, tss, rss, residuals)
    end
    BasicReg(x::Int, f::FormulaTerm{L,R}, sch) where {L, R} = new{L,R,width(sch.rhs)}(x, f)
end


"""
    function BasicReg(
        resp::AbstractVector{Float64},
        pred::AbstractMatrix{Float64},
        yname::String,
        xnames::Vector{String},
        f::FormulaTerm{L,R};
        save_residuals::Bool=false,
        minobs::Int=1
    )::BasicReg{L,R} where {L,R}

## Arguments
- resp::AbstractVector{Float64}: The "Y" or response in a linear regression
- pred::AbstractMatrix{Float64}: The "X" matrix in a linear regression
- yname::String: The name of the response variable
- xnames::Vector{String}: The names of the prediction variables
- f::FormulaTerm{L,R}: A StatsModels.jl formula, saved in the resulting struct
- save_residuals::Bool=false: Whether or not to save the vector of residuals from
    the regression. Note for large numbers of regressions this can significantly slow
    down the speed
- minobs::Int=1: The minimum length of the response vector for the regression to
    run. The regression will also not run if the length of the response vector is
    less than or equal to the number of columns in the prediction matrix.

BasicReg is an intentionally simplistic linear regression. It also attempts to produce
a minimum number of allocations if views of vectors are passed.
"""
function BasicReg(
    resp::AbstractVector,
    pred::AbstractMatrix,
    yname::String,
    xnames::SVector{N, String},
    f::FormulaTerm{L, R};
    save_residuals::Bool=false,
    minobs=1
) where {L, R, N}
    if length(resp) <= size(pred)[2] || length(resp) < minobs
        return BasicReg{N}(length(resp), f)
    end

    coef = cholesky(pred' * pred) \ (pred' * resp)

    BasicReg(
        length(resp),
        f,
        coef,
        xnames,
        yname,
        calc_tss(resp),
        calc_rss(resp, pred, coef),
        save_residuals ? resp - pred * coef : nothing
    )
end

function BasicReg(
    tab::FixedTable{N},
    f::FormulaTerm;
    vargs...
) where {N}
    resp = resp_matrix(tab)
    BasicReg(tab[:, 1], resp, tab.cols[1], resp.cols, f; vargs...)
end


function BasicReg(
    tab::FixedTable{N},
    args...;
    vargs...
) where {N}
    BasicReg(tab[:, 1], resp_matrix(tab), args...; vargs...)
end



"""
    quick_reg(
        data::TimelineTable{false},
        f::FormulaTerm;
        minobs::Real=0.8,
        save_residuals::Bool=false
    )

    quick_reg(
        data::IterateTimelineTable,
        f::FormulaTerm;
        minobs::Real=0.8,
        save_residuals::Bool=false
    )

Calculates a linear regression for the supplied data based on the formula (formula from StatsModels.jl).
Unless the formula explicitly excludes the intercept (i.e., `@formula(y ~ 0 + x)`), an intercept is added.

If `data` is of the type `IterateTimelineTable`, then the formula is applied to each `TimelineTable` in an
optimized way and returns a `Vector{BasicReg}`.

## Arguments
- `minobs::Real`: The minimum number of observations to return a completed regression. If less than 1,
    the value is used as a percentage relative to the total number of business days in the time period.
    Therefore, the default of 0.8 corresponds to at least 80% of the business days over the time period have values.
- `save_residuals::Bool=false`: Whether to save the residuals into `BasicReg`, This can have significant performance implications.


"""
function quick_reg(
    data,
    f::FormulaTerm{L, R};
    minobs=0.8,
    save_residuals::Bool=false
) where {L, R}
    final_data = data[1][data[2:end]..., f]
    N = length(final_data.col_names)
    out = Vector{BasicReg{L, R, N-1}}(undef, length(final_data))
    Threads.@threads for i in 1:length(final_data)
        x = final_data[i]
        out[i] = BasicReg(x, f; minobs=adjust_minobs(minobs, x), save_residuals)
    end
    out
end

function quick_reg(
    data::FixedTable,
    f::FormulaTerm;
    minobs=0.8,
    save_residuals::Bool=false
)
    BasicReg(data, f; minobs=adjust_minobs(minobs, data), save_residuals)
end


StatsBase.predict(mod::BasicReg{L, R, N}, x::FixedTable{N}) where {L, R, N} = x * coef(mod)
StatsBase.predict(mod::BasicReg, x) = x * coef(mod)

StatsBase.coef(x::BasicReg) = isdefined(x, :coef) ? x.coef : missing
StatsBase.coefnames(x::BasicReg) = isdefined(x, :coef) ? x.coefnames : missing
StatsBase.responsename(x::BasicReg) = isdefined(x, :coef) ? x.yname : missing
StatsBase.nobs(x::BasicReg) = x.nobs
StatsBase.dof_residual(x::BasicReg{L, R, N}) where {L, R, N} = isdefined(x, :coef) ? nobs(x) - N : missing
StatsBase.r2(x::BasicReg) = 1 - (rss(x) / deviance(x))
StatsBase.adjr2(x::BasicReg) = 1 - rss(x) / deviance(x) * (nobs(x) - 1) / dof_residual(x)
StatsBase.islinear(x::BasicReg) = true
StatsBase.deviance(x::BasicReg) = isdefined(x, :coef) ? x.tss : missing
StatsBase.rss(x::BasicReg) = isdefined(x, :coef) ? x.rss : missing
function StatsBase.residuals(x::BasicReg)
    if !isdefined(x, :coef)
        return missing
    end
    if x.residuals === nothing
        @error("To access residuals, run `quick_reg` with the option `save_residuals=true`")
    else
        x.residuals
    end
end

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