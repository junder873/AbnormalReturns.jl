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
mutable struct MarketCalendar <: BusinessDays.HolidayCalendar
    bdays::Vector{Date}
    dtmin::Date
    dtmax::Date
    cache::BusinessDays.HolidayCalendarCache
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
function MarketCalendar(bdays::Vector{Date}, dtmin::Date=min(bdays...), dtmax::Date=max(bdays...), _initcache_::Bool=true)
    bdays = sort(bdays)
    market_calendar = MarketCalendar(bdays, dtmin, dtmax, BusinessDays.HolidayCalendarCache())
    market_calendar.cache.hc = market_calendar

    if _initcache_
        BusinessDays.initcache!(market_calendar.cache, market_calendar, dtmin, dtmax)
    end
    return market_calendar
end

@inline BusinessDays.checkbounds(cal::MarketCalendar, dt::Date) = @assert cal.dtmin <= dt && dt <= cal.dtmax "Date out of calendar bounds: $dt. Allowed dates interval is from $(cal.dtmin) to $(cal.dtmax)."

function BusinessDays.isholiday(cal::MarketCalendar, dt::Date)
    BusinessDays.checkbounds(cal, dt)
    return dt ∉ cal.bdays
end

@inline BusinessDays._getcachestate(hc::MarketCalendar) = BusinessDays._getcachestate(hc.cache)
@inline BusinessDays._getholidaycalendarcache(hc::MarketCalendar) = hc.cache
@inline BusinessDays.cleancache(cal::MarketCalendar) = BusinessDays.cleancache!(cal.cache)
@inline BusinessDays.needs_cache_update(hc::MarketCalendar, d0::Date, d1::Date) = BusinessDays._getcachestate(hc) && hc.cache.dtmin == d0 && hc.cache.dtmax == d1

function initcache(hc::MarketCalendar, d0::Date=hc.dtmin, d1::Date=hc.dtmax)
    BusinessDays.checkbounds(hc, d0)
    BusinessDays.checkbounds(hc, d1)
    BusinessDays.initcache!(hc.cache, hc, d0, d1)
end

function BusinessDays.isbday(hc::MarketCalendar, dt::Date)::Bool
    if BusinessDays._getcachestate(hc)
        return isbday(BusinessDays._getholidaycalendarcache(hc), dt)
    else
        return !isholiday(hc, dt)
    end
end