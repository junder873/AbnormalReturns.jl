using CSV, DataFrames, Dates, Test, BusinessDays, StaticArrays
using AbnormalReturns

##


df_firm = CSV.File(joinpath("data", "daily_ret.csv")) |> DataFrame
# one day (Permno 10104, 2019-01-17) is deleted to provide a test for missing data
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_res = CSV.File(joinpath("data", "car_results.csv")) |> DataFrame

##

data = MarketData(
    df_mkt,
    df_firm
)

DataFrames.transform!(data, @formula(mkt ~ mktrf + rf))
##
@test AbnormalReturns.date_range(data.calendar, Date(2019, 3, 31) .. Date(2019, 4, 5)) == 62:66
@test AbnormalReturns.date_range(data.calendar, Date(2019, 3, 31) .. Date(2019, 4, 6)) == 62:66
@test AbnormalReturns.date_range(data.calendar, Date(2019, 4) .. Date(2019, 4, 5)) == 62:66
@test AbnormalReturns.date_range(data.calendar, Date(2019, 4) .. Date(2019, 4, 6)) == 62:66

##
@test isapprox(std(data[18428, Date(2019, 4) .. Date(2019, 10), ["ret"]]), 0.0224085; atol=.00001)
@test isapprox(var(data[18428, Date(2019, 4) .. Date(2019, 10), ["ret"]]), 0.0005021; atol=.00001)
@test isapprox.(std(data[[18428, 18428], [Date(2019, 4), Date(2019, 4)] .. [Date(2019, 10), Date(2019, 10)], ["ret"]]), 0.0224085; atol=.00001) |> all
@test isapprox.(var(data[[18428, 18428], [Date(2019, 4), Date(2019, 4)] .. [Date(2019, 10), Date(2019, 10)], ["ret"]]), 0.0005021; atol=.00001) |> all

@test isapprox(std(data[18428, Date(2019, 4) .. Date(2019, 10), ["ret", "mktrf"]]), 0.0189731; atol=.00001)
@test isapprox(var(data[18428, Date(2019, 4) .. Date(2019, 10), ["ret", "mktrf"]]), 0.00036; atol=.00001)
@test isapprox.(std(data[[18428, 18428], [Date(2019, 4), Date(2019, 4)] .. [Date(2019, 10), Date(2019, 10)], ["ret", "mktrf"]]), 0.0189731; atol=.00001) |> all
@test isapprox.(var(data[[18428, 18428], [Date(2019, 4), Date(2019, 4)] .. [Date(2019, 10), Date(2019, 10)], ["ret", "mktrf"]]), 0.00036; atol=.00001) |> all
##

rr = quick_reg(data[18428, Date(2019, 4) .. Date(2019, 10)], @formula(ret ~ mktrf + hml))
@test coefnames(rr) == SVector{3}("(Intercept)", "mktrf", "hml")
@test responsename(rr) == "ret"
@test nobs(rr) == 126
@test all(isapprox.(coef(rr), [-.00125105, 1.40602071, 1.19924984]; atol=.00001))
@test isapprox(r2(rr), .42340085667)
@test isapprox(adjr2(rr), .41402526)
@test dof_residual(rr) == 123
@test islinear(rr)
@test alpha(rr) == rr.coef[1]
@test beta(rr) == rr.coef[2]

##
rr = quick_reg(data[10104, Date(2019, 1, 2) .. Date(2019, 6)], @formula(ret ~ mktrf + hml))
@test nobs(rr) == 103 # one less due to missing date
@test all(isapprox.(coef(rr), [0.00036009, 0.886629, -0.0247203]; atol=.0001))
@test isapprox(r2(rr), 0.568637; atol=.0001)
##
# test that function term works and interactions work
# lag and lead terms with boundary
rr = quick_reg(data[11762, Date(2019, 1, 2) .. Date(2019, 6)], @formula(ret ~ log1p(mktrf) * lead(mktrf) + hml))
@test nobs(rr) == 104
@test isapprox(adjr2(rr), 0.5872816; atol=.000001)
@test all(isapprox.(coef(rr), [0.00010658, 1.233606, 0.0290669, 0.2685282, -1.80086095]; atol=.000001))
##
rr = quick_reg(data[11762, Date(2019, 1, 2) .. Date(2019, 6)], @formula(ret ~ log1p(mktrf) * lag(mktrf) + hml))
@test nobs(rr) == 103
@test isapprox(adjr2(rr), 0.61329402; atol=.000001)
@test all(isapprox.(coef(rr), [0.000478025, 1.262271804, -0.0727015, 0.3222624, 20.7364917]; atol=.000001))
##
# and again with a missing data point
rr = quick_reg(data[10104, Date(2019, 1, 10) .. Date(2019, 6)], @formula(ret ~ log1p(mktrf) * lead(mktrf) + hml))
@test nobs(rr) == 97
@test all(isapprox.(coef(rr), [0.000156196, 0.85712201, -0.13628144, -0.03920588, -12.1167821]; atol=.000001))
@test isapprox(adjr2(rr), 0.4994264; atol=.000001)
##
# lag and lead terms with boundary
rr = quick_reg(data[10104, Date(2019, 1, 10) .. Date(2019, 6)], @formula(ret ~ log1p(mktrf) * lag(mktrf) + hml))
@test nobs(rr) == 97
@test isapprox(adjr2(rr), 0.48525097; atol=.000001)

##

@test_throws KeyError quick_reg(data[1, Date(2020) .. Date(2021)], @formula(ret ~ mktrf))

rr = quick_reg(data[18428, Date(2019, 4) .. Date(2019, 10)], @formula(ret ~ mktrf + hml); minobs=1000)# not enough data case

@test ismissing(coefnames(rr))
@test ismissing(coef(rr))
@test nobs(rr) == 126
@test ismissing(r2(rr))
@test ismissing(adjr2(rr))
@test alpha(rr) |> ismissing
@test beta(rr) |> ismissing

##

# the SAS code that I tested this against appears to round results to 3 significant digits

# Since these are over specific periods, I specify those as functions to make this easier

event_start(x; data=data) = advancebdays(data.calendar, tobday(data.calendar, x), -10)
event_end(x; data=data) = advancebdays(data.calendar, tobday(data.calendar, x), 10)
est_end(x; data=data) = advancebdays(data.calendar, event_start(x), -16)
est_start(x; data=data) = advancebdays(data.calendar, est_end(x), -149)

rr_market = quick_reg(
    data[df_res.permno, est_start.(df_res.event_date) .. est_end.(df_res.event_date)],
    @formula(ret ~ mktrf)
)

@test isapprox(round.(alpha.(rr_market), digits=5), df_res.alpha_market_model_)
@test isapprox(round.(beta.(rr_market), digits=3), df_res.beta_market_model)
cars = car(data[df_res.permno, event_start.(df_res.event_date) .. event_end.(df_res.event_date)], rr_market)
@test isapprox(round.(cars, sigdigits=3), df_res.car_mm)
bhars = bhar(data[df_res.permno, event_start.(df_res.event_date) .. event_end.(df_res.event_date)], rr_market)
@test isapprox(round.(bhars, sigdigits=3), df_res.bhar_mm)
stds = std(rr_market)
vars = var(rr_market)
@test isapprox(round.(vars, digits=10), df_res.estimation_period_variance_market_model_)

##

rr_ff = quick_reg(
    data[df_res.permno, est_start.(df_res.event_date) .. est_end.(df_res.event_date)],
    @formula(ret ~ mktrf + smb + hml)
)

cars = car(data[df_res.permno, event_start.(df_res.event_date) .. event_end.(df_res.event_date)], rr_ff)
@test isapprox(round.(cars, sigdigits=3), df_res.car_ff)
bhars = bhar(data[df_res.permno, event_start.(df_res.event_date) .. event_end.(df_res.event_date)], rr_ff)
@test isapprox(round.(bhars, sigdigits=3), df_res.bhar_ff)
stds = std.(rr_ff)
vars = var.(rr_ff)
@test isapprox(round.(vars, digits=10), df_res.estimation_period_variance_ff_model_)

##

rr_ffm = quick_reg(
    data[df_res.permno, est_start.(df_res.event_date) .. est_end.(df_res.event_date)],
    @formula(ret ~ mktrf + smb + hml + umd)
)

cars = car(data[df_res.permno, event_start.(df_res.event_date) .. event_end.(df_res.event_date)], rr_ffm)
@test isapprox(round.(cars, sigdigits=3), df_res.car_ffm)
bhars = bhar(data[df_res.permno, event_start.(df_res.event_date) .. event_end.(df_res.event_date)], rr_ffm)
@test isapprox(round.(bhars, sigdigits=3), df_res.bhar_ffm)
stds = std.(rr_ffm)
vars = var.(rr_ffm)
@test isapprox(round.(vars, digits=10), df_res.estimation_period_variance_carhart_model_)

##

