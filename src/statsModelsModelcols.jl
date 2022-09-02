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

StatsModels.modelcols(t::FormulaTerm, d::FixedTable) = (modelcols(t.lhs,d), modelcols(t.rhs, d))
StatsModels.modelcols(t::InterceptTerm{true}, d::FixedTable) = ones(length(d))
StatsModels.modelcols(t::InterceptTerm{false}, d::FixedTable) = Matrix{Float64}(undef, length(d), 0)

function StatsModels.modelcols(terms::MatrixTerm, data::FixedTable)
    reduce(hcat, modelcols.(terms.terms, Ref(data)))
end

StatsModels.modelcols(t::ContinuousTerm, data::FixedTable) = data[:, t.sym]
function StatsModels.modelcols(ft::FunctionTerm{Fo, Fa, Names}, data::FixedTable) where {Fo,Fa,Names}
    ft.fanon.((data[:, n] for n in Names)...)
end

function StatsModels.modelcols(t::InteractionTerm, data::FixedTable)
    StatsModels.row_kron_insideout(*, (modelcols(term, data) for term in t.terms)...)
end

function StatsModels.modelcols(ll::StatsModels.LeadLagTerm{<:Any, F}, data::FixedTable) where {F}
    data[:, F.instance(ll.term.sym, ll.nsteps)]
end

"""
These are largely copied from StatsModels, I just am returning a TimelineColumn and have a special condition
for the lag/lead term
"""
combine_columns(cols) = mapreduce(data_column, union, cols, init=TimelineColumn[])
data_column(::AbstractTerm) = TimelineColumn[]
data_column(t::Union{Term, CategoricalTerm, ContinuousTerm}) = [TimelineColumn(t.sym)]
data_column(t::InteractionTerm) = mapreduce(data_column, union, t.terms)
data_column(t::StatsModels.TupleTerm) = mapreduce(data_column, union, t, init=TimelineColumn[])
data_column(t::MatrixTerm) = data_column(t.terms)
data_column(t::FormulaTerm) = union(data_column(t.lhs), data_column(t.rhs))
data_column(t::FunctionTerm{Fo,Fa,names}) where {Fo,Fa,names} = collect(TimelineColumn.(names))
data_column(t::Symbol) = [TimelineColumn(t)]
data_column(t::TimelineColumn) = [t]
function data_column(t::StatsModels.LeadLagTerm{T, F}) where {T, F<:Union{typeof(lead), typeof(lag)}}
    [F.instance(t.term.sym, t.nsteps)]
end

convert_cols(t::Symbol) = TimelineColumn(t)
convert_cols(t) = t
convert_cols(t::ContinuousTerm) = TimelineColumn(t.sym)
convert_cols(t::StatsModels.LeadLagTerm) = TimelineColumn(t)

StatsModels.term(x::TimelineColumn) = StatsModels.term(x.name)

function StatsModels.schema(f::FormulaTerm, d::MarketData)
    StatsModels.Schema(
        Dict(
            term.(StatsModels.termvars(f)) .=> ContinuousTerm.(StatsModels.termvars(f), 0.0, 0.0, 0.0, 0.0)
        )
    )
end

function StatsModels.schema(f::FormulaTerm, d::FixedTable)
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