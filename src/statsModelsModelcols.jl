"""
To make the regression as fast as possible, it would help to not
recreate a matrix that is very similar to the one that was previously
made. This package is designed around the idea that the LHS variables
will completely change but the RHS variables will often be very
similar since they are from the market data. While the dates
that underly this data might change, the basic structure usually
does not from one regression to the next.

This therefore implements some shortcuts to save/load a matrix for
the RHS variables. If this matrix is already saved and is based
on the same formula, then it loads the saved matrix and slices that
data, resulting in much faster execution times.

In addition, this creates a slightly different interpretation of the
lag/lead operators that StatsModels implements. Typically, a lag
operator would still allow data that is outside the originally
provided dates (i.e., if the requested dates were from 1/1/2021-1/31/2021,
a lag should return a value from before 1/1/2021 if the data exists).
"""

function shift(data::DataVector, shifts::Int, cal::MarketCalendar)
    if shifts == 0
        return data
    end
    # shifts = 2 # shifts = -2
    obj_end_to_cal_end = bdayscount(cal, data.dates.right, cal.dtmax) # 1
    obj_start_to_cal_start = bdayscount(cal, data.dates.left, cal.dtmin) # -1
    dt_min_change = max(shifts, obj_start_to_cal_start) # 2 # -1
    dt_max_change = min(shifts, obj_end_to_cal_end) # 1 # -2
    dt_min = advancebdays(cal, data.dates.left, dt_min_change) # = x.dates.left + 2 # = x.datesleft - 1 = cal.dtmin
    dt_max = advancebdays(cal, data.dates.right, dt_max_change) # = x.dates.right + 1 = cal.dtend # x.dates.right - 2

    new_data = raw_values(data)[1-(shifts - dt_min_change):end-(shifts - dt_max_change)]
    # [1 - (2 - 2):end - (2 - 1)] = [1:end-1]
    # [1 - (-2 - -1):end - (-2 - -2)] = [1 + 1:end]
    new_missings = if data_missing_bdays(data) === nothing
        nothing
    else
        filter(x -> x <= length(new_data), data_missing_bdays(data)) |> Set
    end
    DataVector(new_data, new_missings, dt_min .. dt_max)
end


StatsModels.modelcols(t::FormulaTerm, d::TimelineTableNoMissing) = (modelcols(t.lhs,d), modelcols(t.rhs, d))
StatsModels.modelcols(t::InterceptTerm{true}, d::TimelineTableNoMissing) = ones(length(d))
StatsModels.modelcols(t::InterceptTerm{false}, d::TimelineTableNoMissing) = Matrix{Float64}(undef, length(d), 0)

function StatsModels.modelcols(terms::MatrixTerm, data::TimelineTableNoMissing)
    # I did build another getindex function for this, but I need to combine the missing days from the new dataset
    if data.regression_cache !== nothing && terms == data.regression_cache.terms
        r = date_range(data.calendar, data.regression_cache.dates, data.dates)
        if data.missing_bdays !== nothing
            r = r[Not(data.missing_bdays)]
        end

        return data.regression_cache.data[r, :]
    end
    reduce(hcat, [modelcols(tt, data) for tt in terms.terms])
    
    #modelcols.(terms, Ref(data))
    #out = modelcols(terms, columntable(data))
    # if save_matrix
    #     data.parent.regression_cache = RegressionCache(
    #         terms,
    #         out
    #     )
    # end
    #return out
end

StatsModels.modelcols(t::ContinuousTerm, data::TimelineTableNoMissing) = data[:, t.sym]
function StatsModels.modelcols(ft::FunctionTerm{Fo, Fa, Names}, data::TimelineTableNoMissing) where {Fo,Fa,Names}
    ft.fanon.((data[:, n] for n in Names)...)
end

function StatsModels.modelcols(t::InteractionTerm, data::TimelineTableNoMissing)
    StatsModels.row_kron_insideout(*, (modelcols(term, data) for term in t.terms)...)
end

function StatsModels.modelcols(ll::StatsModels.LeadLagTerm{<:Any, F}, data::TimelineTableNoMissing) where F
    #println(F.instance)
    data[:, F.instance(ll.term.sym, ll.nsteps)]
end

"""
These are largely copied from StatsModels, I just am returning a TimelineColumn and have a special condition
for the lag/lead term
"""
internal_termvars(::AbstractTerm) = TimelineColumn[]
internal_termvars(t::Union{Term, CategoricalTerm, ContinuousTerm}) = [TimelineColumn(t.sym)]
internal_termvars(t::InteractionTerm) = mapreduce(internal_termvars, union, t.terms)
internal_termvars(t::StatsModels.TupleTerm) = mapreduce(internal_termvars, union, t, init=TimelineColumn[])
internal_termvars(t::MatrixTerm) = internal_termvars(t.terms)
internal_termvars(t::FormulaTerm) = union(internal_termvars(t.lhs), internal_termvars(t.rhs))
internal_termvars(t::FunctionTerm{Fo,Fa,names}) where {Fo,Fa,names} = collect(TimelineColumn.(names))
function internal_termvars(t::FunctionTerm{F}) where F<:Union{typeof(lead), typeof(lag)}
    f = nameof(F.instance) == :lead ? lead : lag
    if length(t.args_parsed) == 1
        [f(first(t.args_parsed).sym)]
    elseif length(t.args_parsed) == 2
        [f(first(t.args_parsed).sym, last(t.args_parsed).n)]
    else
        throw(ArgumentError("`$opname` terms require 1 or 2 arguments."))
    end
end

StatsModels.term(x::TimelineColumn) = StatsModels.term(x.name)