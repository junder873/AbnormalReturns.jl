using DataFrames, CSV, DataFramesMeta, Dates, BenchmarkTools, SparseArrays, StaticArrays, LinearAlgebra
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##
@time data = MarketData(df_mkt, df_firm; id_col=:firm_id, valuecols_firms=[:ret])
# First run R7 5700X: 8.356466 seconds (14.89 M allocations: 11.201 GiB, 3.74% gc time, 57.41% compilation time: <1% of which was recompilation)
# Second run R7 5700X: 4.575627 seconds (527.51 k allocations: 10.303 GiB, 19.86% gc time)
# First run i7 6700: 16.871077 seconds (14.90 M allocations: 11.202 GiB, 11.98% gc time, 50.89% compilation time: <1% of which was recompilation)
# Second run i7 6700: 8.927993 seconds (527.49 k allocations: 10.303 GiB, 21.14% gc time, 0.05% compilation time)

##
df_firm = nothing
df_mkt = nothing
GC.gc()

##
@time df_temp = @chain df_events begin
    @transform(:reg_mkt = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt)),)
    @transform(:reg_ffm = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
    @transform(
        :bhar_mkt = bhar(data[:firm_id, :est_window_start .. :est_window_end], :reg_mkt),
        :bhar_ffm = bhar(data[:firm_id, :est_window_start .. :est_window_end], :reg_ffm),
    )
    @rtransform(
        :var_mkt = var(:reg_mkt),
        :var_ffm = var(:reg_ffm),
    )
    @transform(
        :alpha = alpha(:reg_mkt),
        :beta = beta(:reg_mkt),
    )
end
# First run R7 5700X: 9.017057 seconds (30.39 M allocations: 3.184 GiB, 10.33% gc time, 85.86% compilation time)
# Second run R7 5700X: 1.453299 seconds (7.18 M allocations: 1.988 GiB, 13.25% gc time, 5.04% compilation time)
# First run i7 6700: 28.656179 seconds (33.49 M allocations: 3.360 GiB, 53.20% gc time, 90.52% compilation time)
# Second run i7 6700: 3.095861 seconds (7.13 M allocations: 1.985 GiB, 2.87% compilation time)

##


@time @chain df_events begin # parse formula every time and perform all checks each time
    @rtransform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
end
# Run R7 5700X: 26.472138 seconds (360.50 M allocations: 25.194 GiB, 11.69% gc time, 2.24% compilation time)
# i7 6700: 231.401211 seconds (363.87 M allocations: 25.336 GiB, 83.97% gc time, 0.04% compilation time)

##

@time @chain df_events begin # full multithread
    @transform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
end
# Run R7 5700X: 0.437344 seconds (2.07 M allocations: 637.053 MiB, 5.50% compilation time)

##

@time @chain df_events begin # not multi threaded
    @transform(:reg = quick_reg.(data[:firm_id, :est_window_start .. :est_window_end, @formula(ret ~ mkt + smb + hml + umd)], Ref(@formula(ret ~ mkt + smb + hml + umd))),)
end
# Run R7 5700X: 2.890529 seconds (4.91 M allocations: 795.016 MiB, 5.99% gc time, 6.97% compilation time)
# i7 6700: 4.215225 seconds (4.08 M allocations: 752.040 MiB, 1.12% compilation time)
