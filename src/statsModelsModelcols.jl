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
        data_dates(data)
    elseif shifts > 0 # lag
        if shifts > bdayscount(calendar(data), dt_max(data), cal_dt_max(data))
            advancebdays(calendar(data), dt_min(data), shifts) .. cal_dt_max(data)
        else
            advancebdays(calendar(data), dt_min(data), shifts) .. advancebdays(calendar(data), dt_max(data), shifts)
        end
    else # lead
        if shifts < bdayscount(calendar(data), dt_min(data), cal_dt_min(data))
            cal_dt_min(data) .. advancebdays(calendar(data), dt_max(data), shifts)
        else
            advancebdays(calendar(data), dt_max(data), shifts) .. advancebdays(calendar(data), dt_max(data), shifts)
        end
    end

    # # shifts = 2 # shifts = -2
    # obj_end_to_cal_end = bdayscount(data.calendar, dt_max(data), cal_dt_max(data)) # 1
    # obj_start_to_cal_start = bdayscount(data.calendar, dt_min(data), cal_dt_min(data)) # -1
    # dt_min_change = max(shifts, obj_start_to_cal_start) # 2 # -1
    # dt_max_change = min(shifts, obj_end_to_cal_end) # 1 # -2
    # dtmin = advancebdays(data.calendar, dt_min(data), dt_min_change) # = x.dates.left + 2 # = x.datesleft - 1 = cal.dtmin
    # dtmax = advancebdays(data.calendar, dt_max(data), dt_max_change) # = x.dates.right + 1 = cal.dtend # x.dates.right - 2

    # dtmin .. dtmax
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
        if shifts > 0
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

StatsModels.modelcols(t::ContinuousTerm, data::TimelineTable{Mssng}) where {Mssng} = data[:, t.sym]
function StatsModels.modelcols(ft::FunctionTerm{Fo, Fa, Names}, data::TimelineTable) where {Fo,Fa,Names}
    ft.fanon.((data[:, n] for n in Names)...)
end

function StatsModels.modelcols(t::InteractionTerm, data::TimelineTable)
    StatsModels.row_kron_insideout(*, (modelcols(term, data) for term in t.terms)...)
end

function StatsModels.modelcols(ll::StatsModels.LeadLagTerm{<:Any, F}, data::TimelineTable{Mssngs}) where {F, Mssngs}
    data[:, F.instance(ll.term.sym, ll.nsteps)]
end

datavector_modelcols(t::ContinuousTerm, data::TimelineTable) = data[t.sym]
function datavector_modlecols(t::StatsModels.LeadLagTerm{<:Any, F}, data::TimelineTable)  where {F, Mssngs}
    data[F.instance(ll.term.sy, ll.nsteps)]
end
function datavector_modelcols(t::AbstractTerm, data::TimelineTable)
    DataVector(
            modelcols(t, data),
            norm_dates(data),
            calendar(data)
    )
end
function datavector_modelcols(ts::MatrixTerm{Ts}, data::TimelineTable) where {Ts <: Tuple{InterceptTerm{true}, Vararg{ContinuousTerm{Float64}}}}
    datavector_modelcols.(ts.terms, Ref(data))
end
function datavector_modelcols(ts::MatrixTerm{Ts}, data::TimelineTable) where {Ts <: Tuple{InterceptTerm{false}, Vararg{ContinuousTerm{Float64}}}}
    datavector_modelcols.(ts.terms[2:end], Ref(data))
end

function Base.getindex(data::TimelineTable, col::AbstractTerm)
    datavector_modelcols(col, data)
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