
function dates_min_max(dates::Vector{ClosedInterval{Date}})::ClosedInterval{Date}
    maximum(dt_min.(dates)) .. minimum(dt_max.(dates))
end
function dates_min_max(date1::ClosedInterval{Date}, date2::ClosedInterval{Date})::ClosedInterval{Date}
    max(dt_min(date1), dt_min(date2)) .. min(dt_max(date1), dt_max(date2))
end

function range_available(ranges::UnitRange{Int}...)
    max(first.(ranges)) .. min(last.(ranges))
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

adjust_minobs(x::Integer, ::MarketCalendar, ::ClosedInterval{Date}) = x
function adjust_minobs(x::Real, cal::MarketCalendar, dates::ClosedInterval{Date})
    if x < 1
        (bdayscount(cal, dt_min(dates), dt_max(dates)) .+ isbday(cal, dt_max(dates))) * x
    else
        x
    end
end