using DataFrames, CSV, DataFramesMeta, Dates, BenchmarkTools, Cthulhu, SparseArrays, StaticArrays, LinearAlgebra
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##
@time data = MarketData(df_mkt, df_firm; id_col=:firm_id, valuecols_firms=[:ret])
# First run R5 3600: 20.229416 seconds (34.98 M allocations: 13.211 GiB, 12.81% gc time, 64.56% compilation time)
# Second run R5 3600: 8.131775 seconds (681.52 k allocations: 11.431 GiB, 26.59% gc time, 0.06% compilation time)
# First run R7 5700X: 16.141040 seconds (35.39 M allocations: 13.239 GiB, 14.37% gc time, 56.56% compilation time)
# Second run R7 5700X: 5.998991 seconds (547.61 k allocations: 11.431 GiB, 29.92% gc time)
# First run i7 6700: 30.367487 seconds (29.85 M allocations: 12.966 GiB, 40.87% gc time, 38.86% compilation time)
# Second run i7 6700: 10.302273 seconds (688.03 k allocations: 11.436 GiB, 22.65% gc time)

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
# First run R5 3600: 12.886998 seconds (35.27 M allocations: 4.308 GiB, 13.43% gc time, 70.63% compilation time)
# Second run R5 3600: 3.537313 seconds (12.73 M allocations: 3.092 GiB, 24.15% gc time, 2.50% compilation time)
# First run R7 5700X: 9.017057 seconds (30.39 M allocations: 3.184 GiB, 10.33% gc time, 85.86% compilation time)
# Second run R7 5700X: 1.453299 seconds (7.18 M allocations: 1.988 GiB, 13.25% gc time, 5.04% compilation time)
# First run i7 6700: 17.520028 seconds (26.97 M allocations: 4.145 GiB, 28.89% gc time, 77.12% compilation time)
# Second run i7 6700: 5.772432 seconds (4.73 M allocations: 2.944 GiB, 44.27% gc time, 1.59% compilation time)

##


@time @chain df_events[1:1000000, :] begin
    @rtransform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
end
# Run R5 3600: 55.141968 seconds (589.18 M allocations: 59.681 GiB, 12.77% gc time, 0.19% compilation time)
# Run R7 5700X: 26.472138 seconds (360.50 M allocations: 25.194 GiB, 11.69% gc time, 2.24% compilation time)

##

@time @chain df_events[1:1000000, :] begin
    @transform(:reg = quick_reg.(data[:firm_id, :est_window_start .. :est_window_end, @formula(ret ~ mkt + smb + hml + umd)], Ref(@formula(ret ~ mkt + smb + hml + umd))),)
end
# Run R7 5700X: 2.890529 seconds (4.91 M allocations: 795.016 MiB, 5.99% gc time, 6.97% compilation time)
