using DataFrames, CSV, DataFramesMeta, Dates, BenchmarkTools, Cthulhu
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##

@time data = MarketData(df_mkt, df_firm; id_col=:firm_id)
# First run R5 3600: 29.455642 seconds (32.05 M allocations: 13.431 GiB, 7.98% gc time, 38.13% compilation time)
# Second run R5 3600: 17.942238 seconds (669.24 k allocations: 11.774 GiB, 11.31% gc time)
##
df_firm = nothing
df_mkt = nothing
GC.gc()
##
@time df_temp = @chain df_events begin
    @transform(:reg_mkt = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt)),)
    @transform(:reg_ffm = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt + smb + hml + umd)),)
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
# First run R5 3600: 19.867720 seconds (63.84 M allocations: 5.861 GiB, 13.96% gc time, 52.57% compilation time)
# Second run R5 3600: 8.928954 seconds (37.73 M allocations: 4.539 GiB, 17.06% gc time, 0.94% compilation time)

##

cols = TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])
@time @chain df_events[1:100000, :] begin
    @rtransform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end, cols], @formula(ret ~ mkt + smb + hml + umd)),)
end

##
quick_reg(data[1, Date(2017) .. Date(2018), cols], @formula(ret ~ mkt + smb + hml + umd))
##


@benchmark $data[i, TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])] setup=(i=rand(1:10000))
##


i = 151465
ids = df_events.firm_id[1:100]
date_mins = df_events.est_window_start[1:100]
date_maxs = df_events.est_window_end[1:100]
f = @formula(ret ~ mkt + smb + hml + umd)
cols = AbnormalReturns.internal_termvars(f)
data_dict = AbnormalReturns.construct_id_dict(ids)
sch = apply_schema(f, schema(f, data))
cache = AbnormalReturns.create_pred_matrix(data[1], sch)
out = fill(BasicReg(0, f), length(ids))

int_data = data[1]



b = @benchmarkable AbnormalReturns.vector_reg(
    $data,
    $ids,
    $date_mins,
    $date_maxs,
    $f,
)
b.params.evals=100
b.params.seconds=100
run(b)
##
