
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

# These are created to minimize the amount of allocations Julia does
# Julia typically allocates a vector for each loop, which when using
# so many loops, can create real garbage collection problems
# As it turns out, doing sum(abs2, resp - mean(resp)) also does
# an allocation, which could mean allocating a huge amount
# caluclating rss was even worse, so these functions are only
# meant to be used internally but do not allocate if passed a view
function calc_tss(resp)
    out = 0.0
    m = mean(resp)
    @simd for x in resp
        out += (x - m) ^ 2
    end
    out
end
function fast_pred(pred, coef, i)
    out = 0.0
    @simd for j in 1:length(coef)
        @inbounds out += pred[i, j] * coef[j]
    end
    out
end
function calc_rss(resp, pred, coef)
    out = 0.0
    @simd for i in 1:length(resp)
        @inbounds out += (resp[i] - fast_pred(pred, coef, i)) ^ 2
    end
    out
end

function mul_vec(x, y)
    out = 0.0
    @simd for i in 1:length(x)
        @inbounds out += x[i] * y[i]
    end
    out
end
function square_mult(x)
    cols = size(x, 2)
    out = zeros((cols, cols))
    for i in 1:cols
        for j in 1:i
            @views out[i, j] = mul_vec(
                x[:, i],
                x[:, j]
            )
        end
    end
    for i in 1:cols
        for j in 1:i-1
            out[j, i] = out[i, j]
        end
    end

    out
end

function second_mult(x, y)
    cols = size(x, 2)
    out = zeros(cols)
    for i in 1:cols
        @views out[i] = mul_vec(
            x[:, i],
            y
        )
    end
    out
end
function quick_cholesky(x; l=size(x, 1))
    #out = zeros((l, l))
    for i in 1:l
        for j in 1:i
            sum1 = 0.0
            @simd for k in 1:j-1
                sum1 += x[j, k] * x[i, k]
            end
            if j == i
                @inbounds x[j, j] = sqrt(x[j, j] - sum1)
            else
                if x[j, j] > 0
                    @inbounds x[i, j] = (x[i, j] - sum1) / x[j, j]
                end
            end
        end
    end
    LowerTriangular(x)
    #out
end
function my_ldiv!(L::LowerTriangular, y)
    for i in 1:length(y)
        sum1 = 0.0
        @simd for j in 1:i-1
            sum1 += L[i, j] * y[j]
        end
        y[i] = (y[i] - sum1) / L[i, i]
    end
    y
end

function my_ldiv!(L::UpperTriangular, y)
    for i in length(y):-1:1
        sum1 = 0.0
        @simd for j in length(y):-1:i+1
            sum1 += L[i, j] * y[j]
        end
        y[i] = (y[i] - sum1) / L[i, i]
    end
    y
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
    resp::AbstractVector{Float64},
    pred::AbstractMatrix{Float64},
    yname::String,
    xnames::Vector{String},
    f::FormulaTerm{L,R};
    save_residuals::Bool=false,
    minobs=1
)::BasicReg{L,R} where {L,R}
    if length(resp) <= size(pred, 2) || length(resp) < minobs
        return BasicReg(length(resp), f)
    end

    chols = quick_cholesky(square_mult(pred))
    coef = my_ldiv!(chols', my_ldiv!(chols, second_mult(pred, resp)))
    #coef = cholesky!(Symmetric(pred' * pred)) \ (pred' * resp)

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

function create_pred_matrix(data::IterateTimelineTable, sch)
    id = first(data)[1]
    create_pred_matrix(
        parent(data)[id, internal_termvars(sch.rhs), AllowMissing{true}],
        sch
    )
end

function create_pred_matrix(data::TimelineTable{true}, sch)
    # This allows missing to make sure all dates are in the data

    mat = modelcols(sch.rhs, data)
    DataMatrix(
        mat,
        data_dates(data),
        calendar(data)
    )
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
    data::TimelineTable{false},# only works if no missing data for multiplication
    f::FormulaTerm;
    minobs::V=0.8,
    save_residuals::Bool=false
) where {V<:Real}

    if !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
        f = FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs)
    end
    sch = apply_schema(f, schema(f, data))
    select!(data, internal_termvars(sch))

    
    resp, pred = modelcols(sch, data)
    yname, xnames = coefnames(sch)
    BasicReg(
        resp,
        pred,
        yname,
        xnames,
        f;
        save_residuals,
        minobs=adjust_minobs(minobs, calendar(data), norm_dates(data))
    )
end

function fill_vector_reg!(
    out_vector::Vector{BasicReg{L,R}},
    iter_data::IterateTimelineTable{T},
    cache::DataMatrix,
    sch::FormulaTerm,
    f::FormulaTerm{L,R};
    save_residuals::Bool=false,
    minobs=0.8
) where {T, L, R}
    @assert validate_iterator(iter_data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"
    cols = internal_termvars(sch)
    yname, xnames = coefnames(sch)
    
    Threads.@threads for (u_id, iter_indexes) in iter_data
        data = parent(iter_data)[u_id, cols, AllowMissing{true}]
        resp = datavector_modelcols(int_lhs(sch), data)
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
            @inbounds out_vector[iter_index(iter)] = BasicReg(
                x,
                y,
                yname,
                xnames,
                f;
                save_residuals,
                minobs=adjust_minobs(minobs, calendar(parent(data)), iter_dates(iter))
            )
        end
    end
    out_vector
end

function quick_reg(
    data::IterateTimelineTable,
    f::FormulaTerm;
    minobs::V=0.8,
    save_residuals::Bool=false
) where {V<:Real}
    if !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
        f = FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs)
    end
    sch = apply_schema(f, schema(f, parent(data)))
    cache = create_pred_matrix(data, sch)
    out = fill(BasicReg(0, f), total_length(data))
    fill_vector_reg!(
        out,
        data,
        cache,
        sch,
        f;
        save_residuals,
        minobs
    )
end

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
    if !isdefined(x, :coef)
        return missing
    end
    if x.residuals === nothing
        @error("To access residuals, run `cache_reg` with the option `save_residuals=true`")
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