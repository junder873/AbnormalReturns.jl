# this is nearly an exact copy of `GenericHolidayCalendar` that BusinessDays.jl
# creates, except that it allows business dyas to be on the weekend and
# instead of a list of holidays it is a list of business days

"""
    MarketCalendar
* `bdays`: a vector of business days
* `dtmin`: minimum date allowed to check for bdays in bdays set. Defaults to `min(bdays...)`.
* `dtmax`: maximum date allowed to check for bdays in bdays set. Defaults to `max(bdays...)`.
* `cache`: instance of HolidayCalendarCache.
"""
struct MarketCalendar <: BusinessDays.HolidayCalendar
    bdays::Vector{Date}
    dtmin::Date
    dtmax::Date
    isbday_array::Vector{Bool}
    bdayscounter_array::Vector{UInt32}
end

Base.:(==)(g1::MarketCalendar, g2::MarketCalendar) = g1.bdays == g2.bdays && g1.dtmin == g2.dtmin && g1.dtmax == g2.dtmax
Base.hash(g::MarketCalendar) = hash(g.bdays) + hash(g.dtmin) + hash(g.dtmax)

"""
    MarketCalendar(bdays, [dtmin], [dtmax], [_initcache_])
* `bdays`: a vector of dates
* `dtmin`: minimum date allowed to check for bdays in bdays set. Defaults to `min(bdays...)`.
* `dtmax`: maximum date allowed to check for bdays in bdays set. Defaults to `max(bdays...)`.
* `_initcache_`: initializes the cache for this calendar. Defaults to `true`.
"""
function MarketCalendar(bdays::Vector{Date}, dtmin::Date=min(bdays...), dtmax::Date=max(bdays...))
    bdays = sort(bdays)
    isbday_array = zeros(Bool, Dates.value(dtmax - dtmin)+1)
    bdayscounter_array = zeros(UInt32, length(isbday_array))

    for (i, d) in enumerate(dtmin:Day(1):dtmax)
        isbday_array[i] = d ∈ bdays
        if i > 1
            bdayscounter_array[i] = bdayscounter_array[i-1] + isbday_array[i]
        end
    end
    market_calendar = MarketCalendar(bdays, dtmin, dtmax, isbday_array, bdayscounter_array)
    return market_calendar
end

function BusinessDays.checkbounds(cal::MarketCalendar, dt::Date)
    if dt < cal_dt_min(cal)  || dt > cal_dt_max(cal)
        throw(AssertionError("Date out of calendar bounds: $dt. Allowed dates interval is from $(cal.dtmin) to $(cal.dtmax)."))
    end
end

@inline BusinessDays._linenumber(cal::MarketCalendar, dt::Date) = Dates.days(dt) - Dates.days(cal.dtmin) + 1

function BusinessDays.isholiday(cal::MarketCalendar, dt::Date)
    !isbday(cal, dt)
end

function BusinessDays.isbday(hc::MarketCalendar, dt::Date)::Bool
    BusinessDays.checkbounds(hc, dt)
    hc.isbday_array[BusinessDays._linenumber(hc, dt)]
end

function BusinessDays.bdayscount(hc::MarketCalendar, dt0::Date, dt1::Date)::Int
    dt0 = tobday(hc, dt0)
    dt1 = tobday(hc, dt1)
    Int(hc.bdayscounter_array[BusinessDays._linenumber(hc, dt1)]) - Int(hc.bdayscounter_array[BusinessDays._linenumber(hc, dt0)])
end

function BusinessDays.listbdays(hc::MarketCalendar, dt0::Date, dt1::Date)
    BusinessDays.checkbounds(hc, dt0)
    BusinessDays.checkbounds(hc, dt1)
    hc.bdays[dt0 .<= hc.bdays .<= dt1]
end

function date_range(hc::MarketCalendar, dates::ClosedInterval{Date})
    date_pos(hc, dates.left):date_pos(hc, dates.right, true)
end

function date_range(hc::MarketCalendar, dates::ClosedInterval{Vector{Date}})
    @assert length(dates.left) == length(dates.right) "Date sets are not equal length"
    BusinessDays.checkbounds(hc, minimum(dates.left))
    BusinessDays.checkbounds(hc, maximum(dates.left))
    BusinessDays.checkbounds(hc, minimum(dates.right))
    BusinessDays.checkbounds(hc, maximum(dates.right))
    out = Vector{UnitRange{Int}}(undef, length(dates.left))
    @inbounds for i in eachindex(dates.left, dates.right)
        out[i] = date_pos(hc, dates.left[i]; perform_check=false):date_pos(hc, dates.right[i], true; perform_check=false)
    end
    out
end

function date_pos(hc::MarketCalendar, date::Date, sub_end=false; perform_check=true)
    perform_check && BusinessDays.checkbounds(hc, date)
    v = hc.bdayscounter_array[BusinessDays._linenumber(hc, date)]
    if sub_end && !isbday(hc, date)
        v - 1
    else
        v
    end
end


function date_range(hc::MarketCalendar, dt_min::Date, dt_max::Date)
    date_range(hc, dt_min .. dt_max)
end

function Base.show(io::IO, cal::MarketCalendar)
    print(io, "MarketCalendar: $(cal.dtmin) .. $(cal.dtmax) with $(sum(cal.isbday_array)) business days")
end

cal_dt_min(x::MarketCalendar) = x.dtmin
cal_dt_max(x::MarketCalendar) = x.dtmax