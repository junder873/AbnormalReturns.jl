# Iteration Util

struct IterateOutput
    i::Int
    dates::ClosedInterval{Date}
    minobs::Int
end

iter_index(x::IterateOutput) = x.i
iter_dates(x::IterateOutput) = x.dates
iter_minobs(x::IterateOutput) = x.minobs

function construct_id_dict(
    ids::Vector{T},
    date_starts::Vector{Date},
    date_ends::Vector{Date},
    cal::MarketCalendar;
    minobs::V=0.8
) where {T, V<:Real}
    @assert length(ids) == length(date_starts) == length(date_ends) "All Vectors must be of same length"
    out_minobs = if minobs < 1
        bdayscount(cal, date_starts, date_ends) .+ 1
    else
        fill(Int(floor(minobs)), length(ids))
    end
    out = Dict{T, Vector{IterateOutput}}()
    sizehint!(out, length(ids))
    for (i, id) in enumerate(ids)
        if !haskey(out, id)
            out[id] = IterateOutput[]
        end
        push!(out[id], IterateOutput(i, date_starts[i] .. date_ends[i], out_minobs[i]))
    end
    out
end

function validate_iterator(
    data::Dict,
    out_vector::Vector
)
    all_ids = iter_index.(vcat(values(data)...)) |> sort
    all_ids == 1:length(out_vector)
end

function dates_min_max(dates::Vector{ClosedInterval{Date}})::ClosedInterval{Date}
    maximum(dt_min.(dates)) .. minimum(dt_max.(dates))
end
function dates_min_max(date1::ClosedInterval{Date}, date2::ClosedInterval{Date})::ClosedInterval{Date}
    max(dt_min(date1), dt_min(date2)) .. min(dt_max(date1), dt_max(date2))
end

function check_col(x::Symbol, g1, g2)
    x ∈ g1 || x ∈ g2
end

function check_col(x::Vector{Symbol}, g1, g2)
    all(check_col.(x, Ref(g1), Ref(g2)))
end

check_col(x::Symbol, g1) = x ∈ g1

function check_col(x::Vector{Symbol}, g1)
    any(check_col.(x, Ref(g1)))
end