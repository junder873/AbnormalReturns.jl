
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
    resp::AbstractVector{T},
    pred::FixedTable{N, T},
    yname::String,
    xnames::SVector{N, String},
    f::FormulaTerm{L, R};
    save_residuals::Bool=false,
    minobs=1
)::BasicReg{L, R, N} where {L, R, N, T}
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
    args...;
    vargs...
) where {N}
    BasicReg(tab[:, 1], resp_matrix(tab), args...; vargs...)
end



# """
#     quick_reg(
#         data::TimelineTable{false},
#         f::FormulaTerm;
#         minobs::Real=0.8,
#         save_residuals::Bool=false
#     )

#     quick_reg(
#         data::IterateTimelineTable,
#         f::FormulaTerm;
#         minobs::Real=0.8,
#         save_residuals::Bool=false
#     )

# Calculates a linear regression for the supplied data based on the formula (formula from StatsModels.jl).
# Unless the formula explicitly excludes the intercept (i.e., `@formula(y ~ 0 + x)`), an intercept is added.

# If `data` is of the type `IterateTimelineTable`, then the formula is applied to each `TimelineTable` in an
# optimized way and returns a `Vector{BasicReg}`.

# ## Arguments
# - `minobs::Real`: The minimum number of observations to return a completed regression. If less than 1,
#     the value is used as a percentage relative to the total number of business days in the time period.
#     Therefore, the default of 0.8 corresponds to at least 80% of the business days over the time period have values.
# - `save_residuals::Bool=false`: Whether to save the residuals into `BasicReg`, This can have significant performance implications.


# """
# function quick_reg(
#     data::TimelineTable{false},# only works if no missing data for multiplication
#     f::FormulaTerm;
#     minobs::V=0.8,
#     save_residuals::Bool=false
# ) where {V<:Real}

#     if !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
#         return quick_reg(
#             data,
#             FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs);
#             minobs,
#             save_residuals
#         )
#     end
#     sch = apply_schema(f, schema(f, data))
#     select!(data, internal_termvars(sch))

    
#     resp = datavector_modelcols(sch.lhs, data)
#     pred = datavector_modelcols(sch.rhs, data)
#     if nnz(data_missing_bdays(data)) == 0
#         y = view(resp, norm_dates(data))
#         x = FixedWidthMatrix(pred, norm_dates(data))
#     else
#         new_mssngs = get_missing_bdays(data, norm_dates(data))
#         y = view(resp, norm_dates(data), new_mssngs)
#         x = FixedWidthMatrix(pred, norm_dates(data), new_mssngs)
#     end
#     yname = coefnames(sch.lhs)
#     xnames2 = SVector{width(sch.rhs)}(coefnames(sch.rhs))
#     BasicReg(
#         y,
#         x,
#         yname,
#         xnames2,
#         f;
#         save_residuals,
#         minobs=adjust_minobs(minobs, calendar(data), norm_dates(data))
#     )
# end

# function fill_vector_reg!(
#     out_vector::Vector{BasicReg{L,R,N}},
#     iter_data::IterateTimelineTable{T},
#     cache::NTuple{N, DataVector},
#     sch::FormulaTerm,
#     f::FormulaTerm{L,R};
#     save_residuals::Bool=false,
#     minobs=0.8
# ) where {T, L, R, N}
#     @assert validate_iterator(iter_data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"
#     cols = internal_termvars(sch)
#     yname, xnames = coefnames(sch)
#     xnames2 = SVector{N}(xnames)
    
#     Threads.@threads for (u_id, iter_indexes) in iter_data
#         data = parent(iter_data)[u_id, cols, AllowMissing{true}]
#         resp = datavector_modelcols(int_lhs(sch), data)
#         for iter in iter_indexes
#             cur_dates = dates_min_max(data_dates(data), iter_dates(iter))
#             if nnz(data_missing_bdays(data)) == 0
#                 y = view(resp, cur_dates)
#                 x = view(cache, cur_dates)
#             else
#                 new_mssngs = get_missing_bdays(data, cur_dates)
#                 y = view(resp, cur_dates, new_mssngs)
#                 x = view(cache, cur_dates, new_mssngs)
#             end
#             @inbounds out_vector[iter_index(iter)] = BasicReg(
#                 y,
#                 x,
#                 yname,
#                 xnames2,
#                 f;
#                 save_residuals,
#                 minobs=adjust_minobs(minobs, calendar(parent(data)), iter_dates(iter))
#             )
#         end
#     end
#     out_vector
# end

# function general_fill_vector!(
#     out_vector::Vector,
#     iter_data::IterateTimelineTable,
#     cache,#NTuple{N,, DataVector}, DataVector, or nothing
#     cols::Vector{TimelineColumn},
#     col1,# TimelineColumn, Symbol, or abstractTerm
#     fun;
#     minobs,
#     kwargs...
# )
#     @assert validate_iterator(iter_data, out_vector) "Length of out_vector does not match the number of indexes in the iterator"
#     Threads.@threads for (u_id, iter_indexes) in iter_data
#         data = parent(iter_data)[u_id, cols, AllowMissing{true}]
#         resp = data[col1]
#         for iter in iter_indexes
#             cur_dates = dates_min_max(data_dates(data), iter_dates(iter))
#             if nnz(data_missing_bdays(data)) == 0
#                 single_col = view(resp, cur_dates)
#                 cache_cols = view(cache, cur_dates)
#             else
#                 new_mssngs = get_missing_bdays(data, cur_dates)
#                 single_col = view(resp, cur_dates, new_mssngs)
#                 cache_cols = view(cache, cur_dates, new_mssngs)
#             end

#             length(single_col) < adjust_minobs(minobs, calendar(parent(iter_data)), iter_dates(iter)) && continue

#             @inbounds out_vector[iter_index(iter)] = fun(
#                 single_col,
#                 cache_cols;
#                 kwargs...
#             )
#         end
#     end
#     out_vector
# end

# function quick_reg(
#     data::IterateTimelineTable,
#     f::FormulaTerm;
#     minobs::V=0.8,
#     save_residuals::Bool=false
# ) where {V<:Real}
#     if !StatsModels.omitsintercept(f) & !StatsModels.hasintercept(f)
#         f = FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs)
#     end
#     sch = apply_schema(f, schema(f, parent(data)))
#     cache = create_pred_matrix(data, sch)
#     out = fill(BasicReg(0, f, sch), total_length(data))
#     cols = internal_termvars(sch)
#     yname, xnames = coefnames(sch)
#     xnames2 = SVector{N}(xnames)
#     col1 = int_lhs(sch)
#     general_fill_vector!(
#         out,
#         data,
#         cache,
#         cols,
#         col1,
#         BasicReg;
#         save_residuals,
#         minobs,
#         yname=yname,
#         xnames=xnames2,
#         f=f
#     )
# end

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