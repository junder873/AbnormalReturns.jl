using CSV, Dates, BenchmarkTools, InMemoryDatasets
using Revise
using AbnormalReturns

##
ds_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> Dataset
ds_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> Dataset
ds_events = CSV.File(joinpath("data", "event_dates.csv")) |> Dataset |> unique

##
@time ds_all = innerjoin(
    ds_firm,
    ds_mkt,
    on=[:date]
)
# 8.563892 seconds (4.58 k allocations: 7.422 GiB, 19.65% gc time, 0.05% compilation time)
##

@time ds_event_joined = innerjoin(
    ds_all,
    ds_events[:, [:firm_id, :event_date, :est_window_start, :est_window_end]],
    on=[:firm_id => :firm_id, :date => (:est_window_start, :est_window_end)]
)
# 52.283242 seconds (11.54 k allocations: 23.178 GiB, 7.14% gc time)
##
@time groupby!(ds_event_joined, [:firm_id, :event_date])
# 45.565012 seconds (1.01 M allocations: 26.786 GiB, 16.40% gc time)
##
function simple_reg(xs...)
    pred = disallowmissing(hcat(xs[2:end]...))
    resp = disallowmissing(xs[1])
    [cholesky!(Symmetric(pred' * pred)) \ (pred' * resp)]
end

@time InMemoryDatasets.combine(ds_event_joined, (:ret, :mkt, :hml, :umd, :smb) => simple_reg)
# 15.090089 seconds (35.68 M allocations: 19.679 GiB, 56.77% gc time, 8.81% compilation time)

##

@time ds_event_joined = innerjoin(
    ds_all,
    ds_events[:, [:firm_id, :event_date, :event_window_start, :event_window_end]],
    on=[:firm_id => :firm_id, :date => (:event_window_start, :event_window_end)]
)

##

@time groupby!(ds_event_joined, [:firm_id, :event_date])