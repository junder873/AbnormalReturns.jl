
struct TimelineColumn
    name::Symbol
    shifts::Int # negative for leads, positive for lags, 0 for normal
    TimelineColumn(n::Symbol) = new(n, 0)
    TimelineColumn(n::Symbol, s::Int) = new(n, s)
end

struct DictIndex
    cols::Vector{TimelineColumn}
    lookup::Dict{Int, TimelineColumn}# reverse of other itmes since data is just stored in vectors
end

DictIndex(cols::Vector{Symbol}) = DictIndex(cols, zeros(Int, length(cols)))

TimelineColumn(x::TimelineColumn) = x
TimelineColumn(x, s=0) = TimelineColumn(Symbol(x), Int(s))

function DictIndex(cols::Vector{Symbol}, shifts::Vector{Int})
    DictIndex(TimelineColumn.(cols, shifts))
end

function DictIndex(cols::Vector{TimelineColumn})
    DictIndex(
        cols,
        Dict(1:length(cols) .=> cols)
    )
end


function Base.getindex(x::DictIndex, i::Int)
    x.lookup[i]
end

function Base.push!(x::DictIndex, n::Symbol)
    n = TimelineColumn(n)
    push!(x.cols, n)
    x.lookup[length(x.cols)] = n
end

Base.length(x::DictIndex) = length(x.cols)
Base.names(x::DictIndex) = x.cols

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
Base.isequal(x::DictIndex, y::DictIndex) = all(isequal.(names(x), names(y)))

Base.:(==)(x::TimelineColumn, y::TimelineColumn) = isequal(x, y)
Base.:(==)(x::DictIndex, y::DictIndex) = isequal(x, y)

Base.isless(x::TimelineColumn, y::TimelineColumn) = isless(x.name, y.name) && isless(x.shifts, y.shifts)

Base.show(io::IO, x::DictIndex) = show(io, names(x))

Base.show(io::IO, x::TimelineColumn) = print(io, String(x))

ShiftedArrays.lead(x::Symbol, n::Int=1) = TimelineColumn(x, -n)
ShiftedArrays.lag(x::Symbol, n::Int=1) = TimelineColumn(x, n)

function Base.iterate(x::DictIndex, i=1)
    if i > length(x)
        nothing
    else
        (names(x)[i], i+1)
    end
end

Base.convert(::Type{Symbol}, x::TimelineColumn) = Symbol(String(x))
shift_count(x::TimelineColumn) = x.shifts