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

function shift_dates(data::DataVector, shifts::Int)
    if shifts == 0
        return data_dates(data)
    end

    # shifts = 2 # shifts = -2
    obj_end_to_cal_end = bdayscount(data.calendar, dt_max(data), data.calendar.dtmax) # 1
    obj_start_to_cal_start = bdayscount(data.calendar, dt_min(data), data.calendar.dtmin) # -1
    dt_min_change = max(shifts, obj_start_to_cal_start) # 2 # -1
    dt_max_change = min(shifts, obj_end_to_cal_end) # 1 # -2
    dt_min = advancebdays(data.calendar, dt_min(data), dt_min_change) # = x.dates.left + 2 # = x.datesleft - 1 = cal.dtmin
    dt_max = advancebdays(data.calendar, dt_max(data), dt_max_change) # = x.dates.right + 1 = cal.dtend # x.dates.right - 2

    dt_min .. dt_max
end


function shift(data::DataVector, shifts::Int)
    if shifts == 0
        return data
    end
    new_dates = shift_dates(data, shifts)
    len = bdayscount(data.calendar, dt_min(new_dates), dt_max(new_dates)) + 1
    r = if len == length(raw_values(data))
        1:length(raw_values(data))
    else
        if shifts < 0
            1:len
        else
            1+(length(raw_values(data)) - len):length(raw_values(data))
        end
    end
    # [1 - (2 - 2):end - (2 - 1)] = [1:end-1]
    # [1 - (-2 - -1):end - (-2 - -2)] = [1 + 1:end]
    DataVector(raw_values(data)[r], data_missing_bdays(data)[r], new_dates, data.calendar)
end


StatsModels.modelcols(t::FormulaTerm, d::TimelineTable) = (modelcols(t.lhs,d), modelcols(t.rhs, d))
StatsModels.modelcols(t::InterceptTerm{true}, d::TimelineTable) = ones(length(d))
StatsModels.modelcols(t::InterceptTerm{false}, d::TimelineTable) = Matrix{Float64}(undef, length(d), 0)

function StatsModels.modelcols(terms::MatrixTerm, data::TimelineTable)
    reduce(hcat, modelcols.(terms.terms, Ref(data)))
end

StatsModels.modelcols(t::ContinuousTerm, data::TimelineTable)::Vector{Float64} = data[:, t.sym]
function StatsModels.modelcols(ft::FunctionTerm{Fo, Fa, Names}, data::TimelineTable) where {Fo,Fa,Names}
    ft.fanon.((data[:, n] for n in Names)...)
end

function StatsModels.modelcols(t::InteractionTerm, data::TimelineTable)
    StatsModels.row_kron_insideout(*, (modelcols(term, data) for term in t.terms)...)
end

function StatsModels.modelcols(ll::StatsModels.LeadLagTerm{<:Any, F}, data::TimelineTable) where F
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
function internal_termvars(t::StatsModels.LeadLagTerm{T, F}) where {T, F<:Union{typeof(lead), typeof(lag)}}
    [F.instance(t.term.sym, t.nsteps)]
end

StatsModels.term(x::TimelineColumn) = StatsModels.term(x.name)

function StatsModels.schema(f::FormulaTerm, d::MarketData)
    StatsModels.Schema(
        Dict(
            term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0)
        )
    )
end

function StatsModels.schema(f::FormulaTerm, d::TimelineTable)
    StatsModels.Schema(
        Dict(
            term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0)
        )
    )
end

function int_lhs(f::FormulaTerm{L})::L where {L}
    f.lhs
end
function int_rhs(f::FormulaTerm{L,R})::R where {L,R}
    f.rhs
end