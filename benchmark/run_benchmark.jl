using DataFrames, CSV, DataFramesMeta, Dates, BenchmarkTools, Cthulhu, SparseArrays, StaticArrays, LinearAlgebra
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##
@time data = MarketData(df_mkt, df_firm; id_col=:firm_id, valuecols_firms=[:ret])
# First run R5 3600: 17.334942 seconds (29.90 M allocations: 12.970 GiB, 12.24% gc time, 62.19% compilation time)
# Second run R5 3600: 7.513322 seconds (687.97 k allocations: 11.436 GiB, 22.34% gc time)
# First run i7 6700: 37.776454 seconds (32.06 M allocations: 13.430 GiB, 30.23% gc time, 32.16% compilation time)
# Second run i7 6700: 18.069534 seconds (669.95 k allocations: 11.774 GiB, 15.21% gc time)

##
df_firm = nothing
df_mkt = nothing
GC.gc()

##
@time df_temp = @chain df_events begin
    @transform(:reg_mkt = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt)),)
    @transform(:reg_ffm = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
    @transform(
        :bhar_mkt = AbnormalReturns.bhar(data[:firm_id, :est_window_start .. :est_window_end], :reg_mkt),
        :bhar_ffm = AbnormalReturns.bhar(data[:firm_id, :est_window_start .. :est_window_end], :reg_ffm),
        # :bhar_simple = bhar(data[:firm_id, :est_window_start .. :est_window_end], ^(:ret), ^(:mkt)),
        # :car_simple = car(data[:firm_id, :est_window_start .. :est_window_end], ^(:ret), ^(:mkt))
    )
    @rtransform(
        :var_mkt = var(:reg_mkt),
        :var_ffm = var(:reg_ffm),
        :alpha = alpha(:reg_mkt),
        :beta = beta(:reg_mkt),
    )
end
# First run R5 3600: 12.886998 seconds (35.27 M allocations: 4.308 GiB, 13.43% gc time, 70.63% compilation time)
# Second run R5 3600:  3.537313 seconds (12.73 M allocations: 3.092 GiB, 24.15% gc time, 2.50% compilation time)
# First run i7 6700: 33.919134 seconds (64.13 M allocations: 5.936 GiB, 45.39% gc time, 65.25% compilation time)
# Second run i7 6700: 19.343996 seconds (37.73 M allocations: 4.602 GiB, 62.63% gc time, 0.44% compilation time)

##


cols = TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])
@time @chain df_events[1:1000000, :] begin
    @rtransform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end, cols], @formula(ret ~ 1 + mkt + smb + hml + umd)),)
end
# Run R5 3600: 55.141968 seconds (589.18 M allocations: 59.681 GiB, 12.77% gc time, 0.19% compilation time)

##

@time @chain df_events[1:1000000, :] begin
    @transform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
end
# Run R5 3600: 0.655350 seconds (1.53 M allocations: 850.312 MiB, 5.46% compilation time)
##
@time @chain df_events[1:1000000, :] begin
    @transform(:bhar = bhar(data[:firm_id, :est_window_start .. :est_window_end]),)
end
##
quick_reg(data[1, Date(2017) .. Date(2018), cols], @formula(ret ~ mkt + smb + hml + umd))
##


@benchmark $data[i, TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])] setup=(i=rand(1:10000))
##


# i = 151465
# ids = df_events.firm_id[1:100]
# date_mins = df_events.est_window_start[1:100]
# date_maxs = df_events.est_window_end[1:100]
# f = @formula(ret ~ mkt + smb + hml + umd)
# cols = AbnormalReturns.internal_termvars(f)
# data_dict = AbnormalReturns.construct_id_dict(ids)
# sch = apply_schema(f, schema(f, data))
# cache = AbnormalReturns.create_pred_matrix(data[1], sch)
# out = fill(BasicReg(0, f), length(ids))

# int_data = data[1]



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
