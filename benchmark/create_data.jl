using Distributions, DataFrames, DataFramesMeta, Statistics, BusinessDays, Dates, CSV, Arrow

##
BusinessDays.initcache(:USNYSE)

##

dates = listbdays(:USNYSE, Date(1980), Date(2020))

mkt_mat = [
  0.000116248  -1.03659e-5   1.1843e-5   -1.04678e-5;
 -1.03659e-5    3.49562e-5  -2.37018e-6   1.9103e-6;
  1.1843e-5    -2.37018e-6   3.77871e-5  -1.05367e-5;
 -1.04678e-5    1.9103e-6   -1.05367e-5   5.98879e-5;
]

df_mkt_benchmark = DataFrame(
        rand(MvNormal(mkt_mat), length(dates))' |> Matrix, [:mkt, :smb, :hml, :umd])
df_mkt_benchmark[!, :date] = dates
select!(df_mkt_benchmark, [:date, :mkt, :smb, :hml, :umd])

##
count = 10000

coef_means = [
        0.0003416652538979692,
        0.7557301603393168,
        0.6035365140401592,
        0.13841355976862427,
        -0.0917780093857371,
    ]

coef_cov = [
    2.52326e-6  -4.71181e-5   3.31182e-5  -7.05013e-5    1.32415e-5;
    -4.71181e-5   0.27552      0.173033     0.0742338    -0.0448335;
     3.31182e-5   0.173033     0.389159     0.0460528    -0.0280572;
    -7.05013e-5   0.0742338    0.0460528    0.40226      -0.000262376;
     1.32415e-5  -0.0448335   -0.0280572   -0.000262376   0.183439;
    ]
d = Gamma(0.5792772418365709, 0.00448646291687775)

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
firm_errors = rand.(Normal.(0, df_firm_resp.var), nrow(df_mkt_benchmark))


##
df_rand = flatten(DataFrame(
    firm_id = 1:count,
    ret = Ref(main_mkt_mat) .* firm_mats .+ firm_errors
), :ret)

df_rand = @chain df_rand begin
    groupby(:firm_id)
    @transform(:date = df_mkt_benchmark.date)
    select([:firm_id, :date, :ret])
end
##

CSV.write(joinpath("data", "mkt_ret.csv"), df_mkt_benchmark)
CSV.write(joinpath("data", "firm_ret.csv"), df_rand)

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

CSV.write(joinpath("data", "event_dates.csv"), df_events)