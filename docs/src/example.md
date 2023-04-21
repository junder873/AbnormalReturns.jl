# [AbnormalReturns Example](@id Example)

As a quick example:
```@example
data_dir = joinpath("..", "..", "test", "data") # hide
using CSV, DataFramesMeta, Dates, AbnormalReturns

df_firm = CSV.File(joinpath(data_dir, "daily_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath(data_dir, "mkt_ret.csv")) |> DataFrame
df_mkt[!, :mkt] = df_mkt.mktrf .+ df_mkt.rf
df_events = CSV.File(joinpath(data_dir, "firm_earnings_announcements.csv")) |> DataFrame
mkt_data = MarketData(
    df_mkt,
    df_firm
)
df_events = @chain df_events begin
    @rtransform(
        :est_start = advancebdays(mkt_data.calendar, :ea, -120),
        :est_end = advancebdays(mkt_data.calendar, :ea, -2),
        :event_start = advancebdays(mkt_data.calendar, :ea, -1),
        :event_end = advancebdays(mkt_data.calendar, :ea, 1),
    )
    @transform(:reg = quick_reg(mkt_data[:permno, :est_start .. :est_end], @formula(ret ~ mkt + smb + hml)))
    @transform(
        :bhar_reg = bhar(mkt_data[:permno, :event_start .. :event_end], :reg),
        :bhar_simple = bhar(mkt_data[:permno, :event_start .. :event_end, ["ret", "mkt"]]),
        :car_reg = car(mkt_data[:permno, :event_start .. :event_end], :reg),
        :car_simple = car(mkt_data[:permno, :event_start .. :event_end, ["ret", "mkt"]]),
        :total_ret = bh_return(mkt_data[:permno, :event_start .. :event_end, ["ret"]]),
        :total_mkt_ret = bh_return(mkt_data[:permno, :event_start .. :event_end, ["mkt"]]),
    )
    @rtransform(
        :std = std(:reg),
        :var = var(:reg),

    )
    select(Not([:est_start, :est_end, :event_start, :event_end, :reg]))
    # columns eliminated to save space:
    select(Not([:car_reg, :car_simple, :var, :total_mkt_ret]))
end
show(df_events) # hide
```

## Data

For the basic data, this uses the files in the test folder of this package ("test\data"). The "daily\_ret.csv" file is a selection of firm returns, while "mkt\_ret.csv" includes the average market return along with some Fama-French factor returns, you can download similar Fama-French data from [here](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) or from [FamaFrenchData.jl](https://github.com/tbeason/FamaFrenchData.jl) and stock market data from [AlphaVantage.jl](https://github.com/ellisvalentiner/AlphaVantage.jl) or [WRDSMerger.jl](https://github.com/junder873/WRDSMerger.jl) (requires access to the WRDS database).

The firm data uses "Permno" to identify a stock. This package will work with other identifiers, as long as the identifier-date pair is unique. However, Integers and Symbols will often be fastest (as opposed to String identifiers).

```@setup main_run
data_dir = joinpath("..", "..", "test", "data") # hide
using CSV, DataFramesMeta, Dates, AbnormalReturns
```

Load the firm data:
```@example main_run
df_firm = CSV.File(joinpath(data_dir, "daily_ret.csv")) |> DataFrame
show(df_firm) # hide
```

and the market data:
```@example main_run
df_mkt = CSV.File(joinpath(data_dir, "mkt_ret.csv")) |> DataFrame
show(df_mkt) # hide
```

## Arranging and Accessing the Data

Next, load the data into a `MarketData` object:
```@example main_run
mkt_data = MarketData(
    df_mkt,
    df_firm;
    id_col=:permno,# default
    date_col_firms=:date,# default
    date_col_market=:date,# default
    add_intercept_col=true,# default
    valuecols_firms=[:ret],# defaults to nothing, in which case
    # all columns other than id and date are used
    valuecols_market=[:mktrf, :rf, :smb, :hml, :umd]# defaults to
    # nothing, in which case all columns other than date are used
)
show(mkt_data) # hide
```
!!! note
    For performance, especially when loading large datasets of firm data, it is best to make sure the firm dataframe is presorted by ID then Date.

This object rearranges the data so it can be quickly accessed later. The `mkt_data` now contains 3 things:
1. A [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl) calendar that exactly matches the days loaded in the market data.
2. Each column of the `df_mkt` stored
3. Each column of the `df_firm` stored in a `Dict` for each firm.

Data is accessed on a by firm basis, for a given date range and specific columns. For example, say you wanted to get the data for Oracle (ORCL) ("Permno" of 10104), for a specific date (using [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl)) and set of columns:
```@repl main_run
orcl_data = mkt_data[10104, Date(2020) .. Date(2020, 6, 30), [:ret, :mktrf, :smb]]
```

Sometimes it is helpful to add a new column (either for convenience or performance reasons, discussed later). To do so, this package borrows the `transform!` function from [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), using a formula where the left side is the column that is created:
```@repl main_run
transform!(mkt_data, @formula(mkt ~ mktrf + rf));
```

It is also easy to specify the columns as a formula from [StatsModels.jl](https://github.com/JuliaStats/StatsModels.jl). This allows for arbitrary functions, interactions and lags/leads:
```@repl main_run
orcl_data = mkt_data[10104, Date(2020) .. Date(2020, 6, 30), @formula(ret ~ mkt + lag(mkt) + log1p(smb) * hml)]
```
!!! note
    While interactions and arbitrary functions are supported, they can significantly slow down performance since a new vector is allocated in each call. Therefore, it is generally recommended to create a new pice of data by calling `transform!` on the dataset to create the new columns. This advice does not apply to lag/lead terms since those do not need to allocate a new column.

The data returned by accessing `mkt_data` is a `FixedTable`, which is essentially a matrix with a fixed width (useful for multiplication and returning a `StaticMatrix` from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)). Access into this data is done either by a slice as you would any other matrix:
```@repl main_run
orcl_data[:, 1]
```
Or via the names used to access it in the first place:
```@repl main_run
orcl_data[:, :mkt]
```

## Estimating Regressions

The main goal of this package is quickly running regressions for firm events. The example used here is a firm's earnings announcement. Starting with one example, Oracle announced its Q3 2020 earnings on 2020-9-10. Calculating abnormal returns typically follows three steps:

1. Estimate how the firm typically responds to market factors during a control (or estimation) window
2. Use the coefficients from that regression to estimate how the firm should do during the event window
3. Subtract the estimated return from the actual firm return during the event window. Depending on how this difference is aggregated, these are typically buy and hold abnormal returns (bhar) or cumulative abnormla returns (CAR)

First, to create the table for the estimation window, define an estimation window and an event window:
```@repl main_run
event_date = Date("2020-09-10")
est_start = advancebdays(mkt_data.calendar, event_date, -120)
est_end = advancebdays(mkt_data.calendar, event_date, -2)
event_start = advancebdays(mkt_data.calendar, event_date, -1)
event_end = advancebdays(mkt_data.calendar, event_date, 1)
```

Next, run the estimation regression (the regression automatically selects the correct columns from the data, so it is not necessary to do that beforehand):
```@example main_run
orcl_data = mkt_data[10104, est_start .. est_end]
rr = quick_reg(orcl_data, @formula(ret ~ mkt + smb + hml))
```

Then get the data for the event window:
```@repl main_run
orcl_data = mkt_data[10104, event_start .. event_end];
```

Now it is easy to run some statistics for the event window:
```@repl main_run
bhar(orcl_data, rr) # BHAR based on regression

car(orcl_data, rr) # CAR based on regression
```

It is also easy to calculate some statistics for the estimation window:
```@repl main_run
var(rr) # Variance of firm returns (similar equation for standard deviation)

beta(rr) # Firm's market beta

alpha(rr) # Firm's market alpha
```

## More Data Using DataFramesMeta

While the above works well, abnormal returns are often calculated on thousands or more firm-events. Here, I use earnings announcements for about 100 firms from March to November 2020:
```@repl main_run
df_events = CSV.File(joinpath(data_dir, "firm_earnings_announcements.csv")) |> DataFrame
```

Using [DataFramesMeta.jl](https://github.com/JuliaData/DataFramesMeta.jl) and the `@chain` macro from [Chain.jl](https://github.com/jkrumbiegel/Chain.jl), the above steps become:

```@example main_run
df_events = @chain df_events begin
    @rtransform(
        :est_start = advancebdays(mkt_data.calendar, :ea, -120),
        :est_end = advancebdays(mkt_data.calendar, :ea, -2),
        :event_start = advancebdays(mkt_data.calendar, :ea, -1),
        :event_end = advancebdays(mkt_data.calendar, :ea, 1),
    )
    @rtransform(:reg = quick_reg(mkt_data[:permno, :est_start .. :est_end], @formula(ret ~ mkt + smb + hml)))
    @rtransform(
        :bhar_reg = bhar(mkt_data[:permno, :event_start .. :event_end], :reg),
        :bhar_simple = bhar(mkt_data[:permno, :event_start .. :event_end, ["ret", "mkt"]]),
        :std = std(:reg),
        :total_ret = bh_return(mkt_data[:permno, :event_start .. :event_end, ["ret"]]),
    )
    select(Not([:est_start, :est_end, :event_start, :event_end, :reg]))
end
show(df_events) # hide
```

## Vectorizing the Data

While the above works, and is reasonably fast (Doing a test on 1 million regressions takes about 26 seconds on a Ryzen 7 5700X), faster is better.

In particular, a significant reason the above is slow method is that the formula is parsed for each iteration. If the formula is the same for all of the cases, it is better if it is simply parsed once. Therefore, it is optimal to do as much as possible using vectors.

To make this possible, this package provides a type `IterateFixedTable` which will return a `FixedTable` based on a supplied set of ids, dates and columns (or formula as above):
```@repl main_run
est_starts = advancebdays.(mkt_data.calendar, df_events.ea, -120)
est_ends = advancebdays.(mkt_data.calendar, df_events.ea, -2)
vec_data = mkt_data[df_events.permno, est_starts .. est_ends, [:ret, :mkt, :smb]]
```

Each element of `vec_data` is then easily accessible by an integer or can be looped over in a for loop:
```@example main_run
# for x in vec_data
#     x
# end
# or
vec_data[10]
```

This object can be similarly passed to the above functions, just like a firm level table. The function will iterate through the data and return a vector of results.

However, the above is rather ugly. A more practical way to use this is to continue using the `@chain` macro:
```@example main_run
df_events = @chain df_events begin
    @rtransform(
        :est_start = advancebdays(mkt_data.calendar, :ea, -120),
        :est_end = advancebdays(mkt_data.calendar, :ea, -2),
        :event_start = advancebdays(mkt_data.calendar, :ea, -1),
        :event_end = advancebdays(mkt_data.calendar, :ea, 1),
    )
    @transform(:reg = quick_reg(mkt_data[:permno, :est_start .. :est_end], @formula(ret ~ mkt + smb + hml)))
    @transform(
        :bhar_reg = bhar(mkt_data[:permno, :event_start .. :event_end], :reg),
        :bhar_simple = bhar(mkt_data[:permno, :event_start .. :event_end, ["ret", "mkt"]]),
    )
    @transform(
        :std = std.(:reg),
        :total_ret = bh_return(mkt_data[:permno, :event_start .. :event_end, ["ret"]]),
    )
    select(Not([:est_start, :est_end, :event_start, :event_end, :reg]))
end
show(df_events) # hide
```
Notice that the only difference between these two `@chain` macros is that this one uses `@transform` instead of `@rtransform`. This sends the entire column vector to the function, and allows for much faster overall results. Those same 1 million regressions now takes just 0.44 seconds on the same computer.
