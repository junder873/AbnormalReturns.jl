using DataFrames, CSV, DataFramesMeta, Dates
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##

@time data = MarketData(df_mkt, df_firm; id_col=:firm_id)
# First run: 38.816633 seconds (20.05 M allocations: 12.740 GiB, 44.66% gc time, 20.50% compilation time)
# Second run: 18.013269 seconds (629.12 k allocations: 11.719 GiB, 16.72% gc time)

##

@time @chain df_events begin
    @transform(
        :reg_mkt = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt)),
        :reg_ffm = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt + smb + hml + umd)),
    )
    @transform(
        :bhar_mkt = AbnormalReturns.bhar(data, :firm_id, :event_window_start, :event_window_end, :reg_mkt),
        :bhar_ffm = AbnormalReturns.bhar(data, :firm_id, :event_window_start, :event_window_end, :reg_ffm),
    )
    @rtransform(
        :var_mkt = var(:reg_mkt),
        :var_ffm = var(:reg_ffm),
        :alpha = alpha(:reg_mkt),
        :beta = beta(:reg_mkt),
    )
end
# 601.901257 seconds (412.12 M allocations: 52.929 GiB, 95.80% gc time, 0.24% compilation time)

##

@time AbnormalReturns.vector_reg(data, df_events.firm_id, df_events.est_window_start, df_events.est_window_end, @formula(ret ~ mkt))

##

function update_timeline_table(data::AbnormalReturns.TimelineTable{Mssng, T}, id::T, dates::ClosedInterval{Date}) where {Mssng, T}
    data.firmdata = data.parent.firmdata[id]
    data.dates = dates
    data
end

function get_all_items(parent_data::MarketData{T}, ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    data = parent_data[ids[1], date_starts[1] .. date_ends[1], :]
    out = zeros(length(ids))
    for i in 1:length(ids)
        update_timeline_table(data, ids[i], date_starts[i] .. date_ends[i])
        x = data[:, :ret]
        out[i] = sum(x)
    end
    out
end



@time get_all_items(data,  df_events.firm_id, df_events.est_window_start, df_events.est_window_end)

##

function get_all_items2(parent_data::MarketData{T}, ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    out = zeros(length(ids))
    for i in 1:length(ids)
        data = parent_data[ids[i], date_starts[i] .. date_ends[i], [:ret]]
        x = data[:, :ret]
        out[i] = sum(x)
    end
    out
end

#@code_warntype get_all_items2(data,  df_events.firm_id[1:10], df_events.est_window_start[1:10], df_events.est_window_end[1:10])
@time get_all_items2(data,  df_events.firm_id, df_events.est_window_start, df_events.est_window_end)

##

function get_ind_column(parent_data::MarketData{T}, ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    out = zeros(length(ids))
    for i in 1:length(ids)
        x = parent_data.firmdata[ids[i]].ret
        r = AbnormalReturns.date_range(parent_data.calendar, x.dates, date_starts[i] .. date_ends[i])
        out[i] = sum(@view x.data[r])
    end
    out
end

@time get_ind_column(data,  df_events.firm_id, df_events.est_window_start, df_events.est_window_end)

##
@code_warntype get_ind_column(data,  df_events.firm_id[1:10], df_events.est_window_start[1:10], df_events.est_window_end[1:10])

##

@code_warntype data[500, Date(2015) .. Date(2016), [TimelineColumn(:ret)]]