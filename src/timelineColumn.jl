
struct TimelineColumn
    name::Symbol
    shifts::Int # negative for leads, positive for lags, 0 for normal
    TimelineColumn(n::Symbol) = new(n, 0)
    TimelineColumn(n::Symbol, s::Int) = new(n, s)
end

TimelineColumn(x::TimelineColumn) = x
TimelineColumn(x, s=0) = TimelineColumn(Symbol(x), Int(s))

Base.Symbol(x::TimelineColumn) = x.name
function Base.String(x::TimelineColumn)
    if x.shifts == 0
        x.name |> string
    elseif x.shifts < 0
        v = x.shifts < -1 ? ", $(-1 * x.shifts)" : ""
        "lead($(x.name)$v)"
    else
        v = x.shifts > 1 ? ", $(x.shifts)" : ""
        "lag($(x.name)$v)"
    end
end

function Base.string(x::TimelineColumn)
    String(x)
end

Base.isequal(x::TimelineColumn, y::TimelineColumn) = x.name == y.name && x.shifts == y.shifts

Base.:(==)(x::TimelineColumn, y::TimelineColumn) = isequal(x, y)

Base.isless(x::TimelineColumn, y::TimelineColumn) = isless(x.name, y.name) && isless(x.shifts, y.shifts)

Base.show(io::IO, x::TimelineColumn) = print(io, String(x))

ShiftedArrays.lead(x::Symbol, n::Int=1) = TimelineColumn(x, -n)
ShiftedArrays.lag(x::Symbol, n::Int=1) = TimelineColumn(x, n)

Base.convert(::Type{Symbol}, x::TimelineColumn) = Symbol(String(x))
shift_count(x::TimelineColumn) = x.shifts

TimelineColumn(t::StatsModels.LeadLagTerm{<:Any, F}) where {F<:Union{typeof(lead), typeof(lag)}} = F.instance(t.term.sym, t.nsteps)
TimelineColumn(t::Union{ContinuousTerm, InteractionTerm, FunctionTerm}) = TimelineColumn(coefnames(t))