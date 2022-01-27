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


# function StatsModels.modelcols(ll::StatsModels.LagLeadTerm{<:Any, F}, data::TimelineTableNoMissing) where F
#     x = data[ll.term.sym]
#     if F == lag
#         l_change = ll.nsteps
#         r_change = min(ll.nsteps, bdayscount(data.calendar, x.dates.right, data.celendar.dtmax))
#         lost_obs = r_change - l_change
#     else
#         l_change = -1 * min(ll.nsteps, bdayscount(data.calendar, data.calendar.dtmin, x.dates.left))
#         r_change = -1 * ll.nsteps
#         lost_obs = l_change - r_change
#     end
#     dt_min = advancebdays(data.calendar, x.dates.left, l_change)
#     dt_max = advancebdays(data.calendar, x.dates.right, r_change)
#     # need to deal with case that there are not enough extra days in calendar, in which case
#     # remove some from the end
#     new_data = if F == lag # if lost_obs == 0, then just have the same vector
#         raw_values(x)[1:end-lost_obs]
#     else
#         raw_values(x)[1+lost_obs:end]
#     end
#     new_missings = if x.missing_bdays === nothing
#         nothing
#     else
#         filter(x -> x <= length(new_data), x.missing_bdays) |> Set
#     end
#     DataVector(new_data, new_missings, dt_min .. dt_max)
# end

# function StatsModels.modelcols(t::ContinuousTerm, d::TimelineTableNoMissing)
#     return d[t.sym]
# end

# function StatsModels.modelcols(t::InteractionTerm, d::TimelineTableNoMissing)
#     return *([modelcols(x, d) for x in t.terms]...)
# end

# function StatsModels.modelcols(t::InterceptTerm{true}, d::TimelineTableNoMissing)
#     dates = d.calendar.dtmin .. d.calendar.dtmax
#     l = bdayscount(d.calendar, dates.left, dates.right) + 1
#     return DataVector(ones(l), nothing, dates, d.calendar)
# end

function StatsModels.modelcols(terms::MatrixTerm, data::TimelineTableNoMissing; save_matrix::Bool=true)
    if data.regression_cache !== nothing && terms == data.regression_cache.terms

        return data.regression_cache.rhs_mat
    end
    out = hcat([modelcols(tt, data) for tt in terms.terms]...)
    #println("here2")
    if save_matrix
        data.parent.regression_cache = RegressionCache(
            terms,
            out
        )
    end
    return out
end

# function Base.hcat(data::DataVector...)
#     dates = dates_min_max(data_dates.(data)...)
#     data = [y[dates] for y in data]
#     # println(dates)
#     # println(data_dates.(data))
#     # println(length(data))
#     # println(length.(data))
#     mat = hcat(raw_values.(data)...)
#     # display(mat)
#     # println(typeof(mat))
#     missing_bdays = combine_missing_bdays(data_missing_bdays.(data)...)
#     #display(missing_bdays)
#     return DataMatrix(mat, missing_bdays, dates, data[1].calendar)
# end



# function Base.:(*)(x::DataVector...)
#     dates = dates_min_max(data_dates.(x)...)
#     x = (i[dates] for i in x)
#     vec = (*).(raw_values.(x)...)
#     missing_bdays = combine_missing_bdays(data_missing_bdays.(x)...)
#     DataVector(vec, missing_bdays, dates, x[1].calendar)
# end

# function Base.:(+)(x::DataVector...)
#     dates = dates_min_max(data_dates.(x)...)
#     x = (i[dates] for i in x)
#     vec = (+).(raw_values.(x)...)
#     missing_bdays = combine_missing_bdays(data_missing_bdays.(x)...)
#     DataVector(vec, missing_bdays, dates, x[1].calendar)
# end

# function Base.:(-)(x::DataVector...)
#     dates = dates_min_max(data_dates.(x)...)
#     x = (i[dates] for i in x)
#     vec = (-).(raw_values.(x)...)
#     missing_bdays = combine_missing_bdays(data_missing_bdays.(x)...)
#     DataVector(vec, missing_bdays, dates, x[1].calendar)
# end

# function Base.:(/)(x::DataVector...)
#     dates = dates_min_max(data_dates.(x)...)
#     x = (i[dates] for i in x)
#     vec = (/).(raw_values.(x)...)
#     missing_bdays = combine_missing_bdays(data_missing_bdays.(x)...)
#     DataVector(vec, missing_bdays, dates, x[1].calendar)
# end
