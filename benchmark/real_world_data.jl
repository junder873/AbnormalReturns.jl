# This file runs the benchmark using real world data (which includes missing data)
# the goal is partially to test that the functions work properly

# the data used in this test is proprietary, but it is data downloaded from the WRDS database
# and is the CRSP daily file and the Fama French factors file, the CRSP daily file goes back
# well before the test starts (test starts in 2015, CRSP daily file goes to 1986 in this version)
# and Fama French file goes back to 1926

using CSV, DataFramesMeta, Dates
using Revise
using AbnormalReturns

df_crsp = CSV.File(joinpath("data", "crsp_entry.csv")) |> DataFrame
df_mkt = @chain CSV.File(joinpath("data", "ff_entry.csv")) begin
    DataFrame
    @rtransform(:mkt = :mktrf + :rf)
    select([:date, :mkt, :mktrf, :smb, :hml, :umd, :rf])
end
df_res = CSV.File(joinpath("data", "results.csv")) |> DataFrame
mkt_data = MarketData(df_mkt, df_crsp)

##

df_test = @chain df_crsp begin
    @rsubset(Date(2015) <= :date <= Date(2019))
    @by(
        :permno,
        :date = minimum(:date):Day(1):maximum(:date) # redo to every date to include weekends
    )
end

##

@time df_test = @chain df_test begin
    @transform(:reg = quick_reg(mkt_data[:permno, :date-Year(1) .. :date-Day(1)], @formula(ret ~ mkt)))
    @transform(
        :bh_ret = bh_return(mkt_data[:permno, :date .. :date + Week(2), ["ret"]]),
        :car_raw = car(mkt_data[:permno, :date .. :date + Week(2), ["ret", "mkt"]]),
        :bhar_raw = bhar(mkt_data[:permno, :date .. :date + Week(2), ["ret", "mkt"]]),
        :car_mm = car(mkt_data[:permno, :date .. :date + Week(2)], :reg),
        :bhar_mm = bhar(mkt_data[:permno, :date .. :date + Week(2)], :reg),
        :r2_mm = r2.(:reg),
        :std_mm = std.(:reg)
    )
    select(Not(:reg))
    @transform(:reg = quick_reg(mkt_data[:permno, :date-Year(1) .. :date-Day(1)], @formula(ret ~ mkt + smb + umd + hml)))
    @transform(
        :car_ffm = car(mkt_data[:permno, :date .. :date + Week(2)], :reg),
        :bhar_ffm = bhar(mkt_data[:permno, :date .. :date + Week(2)], :reg),
        :r2_ffm = r2.(:reg),
        :std_ffm = std.(:reg)
    )
    select(Not(:reg))
    @transform(:reg = quick_reg(mkt_data[:permno, :date-Year(1) .. :date-Day(1)], @formula(ret ~ 0 + mkt + smb + umd + hml)))
    @transform(
        :car_ffm2 = car(mkt_data[:permno, :date .. :date + Week(2)], :reg),
        :bhar_ffm2 = bhar(mkt_data[:permno, :date .. :date + Week(2)], :reg),
        :r2_ffm2 = r2.(:reg),
        :std_ffm2 = std.(:reg)
    )
    select(Not(:reg))
end
# A little over 10 million observations
# First Run R7 5700X: 38.191368 seconds (200.79 M allocations: 55.747 GiB, 24.94% gc time, 215.70% compilation time: <1% of which was recompilation)
# Second Run: 35.790831 seconds (182.35 M allocations: 60.448 GiB, 37.53% gc time, 1.01% compilation time)
# changing sparsevec to vector in iteratetimelinetable reduced time by about 10 seconds
# dealing with missing values at the start or end of a range in the setup saves another few seconds
# changing the r[Not(missing_bdays)] to a custom function that creates the vector first saves around 30 seconds

# the main cost are the cases where there is a single or short set of missing values, getting rid of those saves over 30 seconds
# this is largely due to setting up the data where it creates a view that is noncontinuous, which means creating the entire vector
# that way

##

approx_or_missing(x::Missing, y::Missing) = true
approx_or_missing(x::Real, y::Real) = isapprox(x, y)
approx_or_missing(x, y) = false

df_res2 = @chain df_res begin
    outerjoin(
        _,
        df_test,
        on=[:permno, :date],
        validate=(true, true),
        makeunique=true
    )
end


df_res2[!, :all_equal] .= true
for col in names(df_res)
    println(col)
    col ∈ ("permno", "date") && continue
    col1 = col * "_1"
    col1 ∉ names(df_res2)
    df_res2[!, :all_equal] = df_res2.all_equal .* approx_or_missing.(df_res2[:, col], df_res2[:, col1])
end

@chain df_res2 begin
    @rsubset(!:all_equal)
end
