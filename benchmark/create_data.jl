using Distributions, DataFrames, DataFramesMeta, Statistics, BusinessDays, Dates, CSV
using Revise
using AbnormalReturns

##
BusinessDays.initcache(:USNYSE)

##

dates = listbdays(:USNYSE, Date(1980), Date(2020))

df_mkt_benchmark = DataFrame(
        rand(MvNormal(mkt_mat), length(dates))' |> Matrix, [:mkt, :smb, :hml, :umd])
df_mkt_benchmark[!, :date] = dates
select!(df_mkt_benchmark, [:date, :mkt, :smb, :hml, :umd])

##
count = 10000

df_firm_resp = DataFrame(hcat(
    1:count |> collect,
    rand(MvNormal(coef_means, coef_cov), count)' |> Matrix,
    rand(d, count)
), [:firm_id, :int, :mkt, :smb, :hml, :umd, :var]
)
@transform!(df_firm_resp, :firm_id = Int.(:firm_id))
##

main_mkt_mat = df_mkt_benchmark[:, [:mkt, :smb, :hml, :umd]] |> Matrix
firm_mats = collect.(Tuple.(Tables.rowtable(df_firm_resp[:, [:mkt, :smb, :hml, :umd]])))
firm_errors = rand.(Normal.(df_firm_resp.var), nrow(df_mkt_benchmark))


##
df_rand = flatten(DataFrame(
    firm_id = 1:count,
    ret = Ref(main_mkt_mat) .* firm_mats .+ firm_errors
), :ret)
##
df_rand = @chain df_rand begin
    groupby(:firm_id)
    @transform(:date = df_mkt_benchmark.date)
    select([:firm_id, :date, :ret])
end
##

CSV.write(joinpath("ab_benchmark_data", "mkt_ret.csv"), df_mkt_benchmark)
CSV.write(joinpath("ab_benchmark_data", "firm_ret.csv"), df_rand)

##

df_events = DataFrame(
    firm_id = rand(1:count, 1_000_000),
    event_date = rand(Date(1982):Day(1):Date(2019, 12), 1_000_000)
)
df_events = @chain df_events begin
    @rtransform(
        :event_window_start = advancebdays("USNYSE", :event_date, -10),
        :event_window_end = advancebdays("USNYSE", :event_date, 10),
        :est_window_start = :event_date - Year(1),
        :est_window_end = advancebdays("USNYSE", :event_date, -11),
    )
end

CSV.write(joinpath("ab_benchmark_data", "event_dates.csv"), df_events)