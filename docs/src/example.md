# Example

```julia
using CSV, DataFramesMeta, Dates, AbnormalReturns
```




As a quick example:
```julia
df_firm = CSV.File(joinpath("data", "daily_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_mkt[!, :mkt] = df_mkt.mktrf .+ df_mkt.rf
df_events = CSV.File(joinpath("data", "firm_earnings_announcements.csv")) |> DataFrame
mkt_data = MarketData(
    df_mkt,
    df_firm
)
@chain df_events begin
    @rtransform(
        :est_start = advancebdays(mkt_data.calendar, :ea, -120),
        :est_end = advancebdays(mkt_data.calendar, :ea, -2),
        :event_start = advancebdays(mkt_data.calendar, :ea, -1),
        :event_end = advancebdays(mkt_data.calendar, :ea, 1),
    )
    @transform(:reg = quick_reg(mkt_data[:permno, :est_start .. :est_end], @formula(ret ~ mkt + smb + hml)))
    @transform(
        :bhar_reg = bhar(mkt_data[:permno, :event_start .. :event_end], :reg),
        :bhar_simple = bhar(mkt_data[:permno, :event_start .. :event_end], "ret", "mkt"),
        :car_reg = car(mkt_data[:permno, :event_start .. :event_end], :reg),
        :car_simple = car(mkt_data[:permno, :event_start .. :event_end], "ret", "mkt"),
    )
    @rtransform(
        :std = std(:reg),
        :var = var(:reg),
        :total_ret = bh_return(mkt_data[:permno, :event_start .. :event_end], "ret"),
        :total_mkt_ret = bh_return(mkt_data[:permno, :event_start .. :event_end], "mkt"),
    )
    select(Not([:est_start, :est_end, :event_start, :event_end, :reg]))
    # columns eliminated to save space:
    select(Not([:car_reg, :car_simple, :var, :total_mkt_ret]))
end
```

```
283×6 DataFrame
 Row │ permno  ea          bhar_reg     bhar_simple  std         total_ret
     │ Int64   Date        Float64?     Float64      Float64?    Float64
─────┼─────────────────────────────────────────────────────────────────────
───
   1 │  49373  2020-03-05  -0.0114753   -0.0368966   0.0120205   -0.0494716
   2 │  23660  2020-03-19  -0.0761877   -0.0799651   0.0100327   -0.16273
   3 │  17144  2020-03-18   0.00584866  -0.00749583  0.0122784    0.0056200
7
   4 │  52708  2020-03-19   0.0449673    0.057314    0.023234    -0.0254504
   5 │  57665  2020-03-24   0.105811     0.0931054   0.0105352    0.171386
   6 │  61621  2020-03-25   0.125142     0.129935    0.0132682    0.303036
   7 │  10104  2020-03-12   0.0295405    0.0515017   0.00903142  -0.01338
   8 │  75154  2020-03-19   0.15375      0.0269029   0.0312084   -0.0558615
  ⋮  │   ⋮         ⋮            ⋮            ⋮           ⋮            ⋮
 277 │  90373  2020-10-29  -0.00714584  -0.00486512  0.016151    -0.0422143
 278 │  90386  2020-11-02  -0.120063    -0.0769287   0.0256068   -0.0606557
 279 │  90993  2020-10-29  -0.00087777   0.00495886  0.0109232   -0.0323903
 280 │  90880  2020-10-28   0.00943539  -0.00174027  0.0109558   -0.0271723
 281 │  91668  2020-10-30   0.0290468    0.0200758   0.017888     0.0283726
 282 │  12449  2020-11-05  -0.0402344   -0.0572696   0.0183221   -0.0128859
 283 │  61399  2020-11-18  -0.0761987   -0.0706244   0.0134081   -0.0759727
                                                              268 rows omit
ted
```

The below sections step through this example more in depth.



## Data

For the basic data, this uses the files in the test folder of this package ("test\data"). The "daily_ret.csv" file is a selection of firm returns, while "mkt_ret.csv" includes the average market return along with some Fama-French factor returns, you can download similar Fama-French data from [here](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) and stock market data from [AlphaVantage.jl](https://github.com/ellisvalentiner/AlphaVantage.jl) or [WRDSMerger.jl](https://github.com/junder873/WRDSMerger.jl) (requires access to the WRDS database).

The firm data uses "Permno" to identify a stock. This package will work with other identifiers, as long as the identifier-date pair is unique.

Load the firm data:
```julia
df_firm = CSV.File(joinpath("data", "daily_ret.csv")) |> DataFrame
```

```
48659×4 DataFrame
   Row │ permno  date        ret           vol
       │ Int64   Date        Float64?      Float64
───────┼──────────────────────────────────────────────────
     1 │  10104  2019-01-02   0.00155038        1.43204e7
     2 │  10104  2019-01-03  -0.00973026        1.98687e7
     3 │  10104  2019-01-04   0.0430996         2.0984e7
     4 │  10104  2019-01-07   0.0158425         1.79679e7
     5 │  10104  2019-01-08   0.00906218        1.62557e7
     6 │  10104  2019-01-09  -0.0020886         1.91002e7
     7 │  10104  2019-01-10   0.00083719        1.66527e7
     8 │  10104  2019-01-11   0.00982855        1.63975e7
   ⋮   │   ⋮         ⋮            ⋮              ⋮
 48653 │  91668  2020-12-22   0.0112328    332197.0
 48654 │  91668  2020-12-23   0.00639975   182199.0
 48655 │  91668  2020-12-24   0.00367913    39932.0
 48656 │  91668  2020-12-28   0.0184641    162785.0
 48657 │  91668  2020-12-29  -0.0155966    140853.0
 48658 │  91668  2020-12-30   0.0124131    119724.0
 48659 │  91668  2020-12-31  -0.00222926   196481.0
                                        48644 rows omitted
```





and the market data:
```julia
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_mkt[!, :mkt] = df_mkt.mktrf .+ df_mkt.rf
df_mkt
```

```
672×7 DataFrame
 Row │ date        mktrf    smb      hml      rf       umd      mkt
     │ Date        Float64  Float64  Float64  Float64  Float64  Float64
─────┼──────────────────────────────────────────────────────────────────
   1 │ 2019-01-02   0.0023   0.006    0.0113   0.0001  -0.023    0.0024
   2 │ 2019-01-03  -0.0245   0.0036   0.012    0.0001  -0.0078  -0.0244
   3 │ 2019-01-04   0.0355   0.0041  -0.007    0.0001  -0.0097   0.0356
   4 │ 2019-01-07   0.0094   0.0101  -0.0074   0.0001  -0.0077   0.0095
   5 │ 2019-01-08   0.0101   0.0054  -0.0063   0.0001   0.0011   0.0102
   6 │ 2019-01-09   0.0056   0.0046   0.001    0.0001  -0.0083   0.0057
   7 │ 2019-01-10   0.0042   0.0003  -0.0046   0.0001  -0.0037   0.0043
   8 │ 2019-01-11  -0.0001   0.0012   0.0022   0.0001  -0.002    0.0
  ⋮  │     ⋮          ⋮        ⋮        ⋮        ⋮        ⋮        ⋮
 666 │ 2021-08-23   0.0108   0.0092  -0.0083   0.0      0.0076   0.0108
 667 │ 2021-08-24   0.0041   0.0034   0.0038   0.0      0.0094   0.0041
 668 │ 2021-08-25   0.0026  -0.0004   0.0034   0.0      0.0072   0.0026
 669 │ 2021-08-26  -0.0066  -0.004   -0.0031   0.0     -0.0031  -0.0066
 670 │ 2021-08-27   0.0108   0.0162   0.003    0.0      0.0094   0.0108
 671 │ 2021-08-30   0.0035  -0.0035  -0.0138   0.0     -0.0117   0.0035
 672 │ 2021-08-31  -0.0013   0.0043  -0.0005   0.0     -0.002   -0.0013
                                                        657 rows omitted
```





## Arranging and Accessing the Data

Next, load the data into a `MarketData` object:
```julia
mkt_data = MarketData(
    df_mkt,
    df_firm;
    id_col=:permno,# default
    date_col_firms=:date,# default
    date_col_market=:date,# default
    valuecols_firms=[:ret],# defaults to nothing, in which case
    # all columns other than id and date are used
    valuecols_market=[:mkt, :smb, :hml, :umd]# defaults to
    # nothing, in which case all columns other than date are used
)
```

```
MarketData with ID type Int64 with 99 unique firms
MarketCalendar: 2019-01-02 .. 2021-08-31 with 672 business days
Market Columns: mkt, smb, hml, umd
Firm Columns: ret
```




!!! note
    For performance, especially when loading large datasets of firm data, it is best to make sure the firm dataframe is presorted by ID then Date.

This object rearranges the data so it can be quickly accessed later. The `mkt_data` now contains 3 things:
1. A [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl) calendar that exactly matches the days loaded in the market data.
2. Each column of the `df_mkt` stored
3. Each column of the `df_firm` stored in a `Dict` for each firm.

Data is accessed on a by firm basis. For example, the "Permno" for Oracle (ORCL) is 10104:
```julia
orcl_data = mkt_data[10104]
```

```
TimelineTable for 10104 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-01-02..2020-12-31
Currently selected dates: 2019-01-02..2020-12-31
============= ============== ========= ========= ========= ==========
        Date            ret       mkt       smb       hml       umd
============= ============== ========= ========= ========= ==========
  2019-01-02     0.00155038    0.0024     0.006    0.0113    -0.023
  2019-01-03    -0.00973026   -0.0244    0.0036     0.012   -0.0078
  2019-01-04      0.0430996    0.0356    0.0041    -0.007   -0.0097
  2019-01-07      0.0158425    0.0095    0.0101   -0.0074   -0.0077
  2019-01-08     0.00906218    0.0102    0.0054   -0.0063    0.0011
  2019-01-09     -0.0020886    0.0057    0.0046     0.001   -0.0083
  2019-01-10     0.00083719    0.0043    0.0003   -0.0046   -0.0037
  2019-01-11     0.00982855       0.0    0.0012    0.0022    -0.002
  2019-01-14    -0.00227792   -0.0059    -0.006    0.0094   -0.0059
  2019-01-15     0.00809466    0.0107       0.0   -0.0088    0.0114
  2019-01-16     -0.0066143    0.0029    0.0003    0.0092   -0.0076
  2019-01-17      0.0108198    0.0076    0.0011   -0.0024    -0.002
  2019-01-18      0.0142033     0.013   -0.0037    0.0011   -0.0085
  2019-01-22    -0.00669782   -0.0152   -0.0041    0.0032    0.0094
  2019-01-23     0.00613002    0.0016   -0.0039   -0.0014    0.0099
  2019-01-24   -0.000812366    0.0024    0.0044       0.0   -0.0102
      ⋮             ⋮            ⋮         ⋮         ⋮         ⋮
============= ============== ========= ========= ========= ==========
                                                     489 rows omitted
```





You can also request a specific range of dates and columns using [EllipsisNotation.jl](https://github.com/ChrisRackauckas/EllipsisNotation.jl):
```julia
orcl_data = mkt_data[10104, Date(2020) .. Date(2020, 6, 30), [:ret, :mkt, lag(:mkt)]]
```

```
TimelineTable for 10104 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-01-03..2020-12-31
Currently selected dates: 2020-01-01..2020-06-30
============= ============= ========== ===========
        Date           ret        mkt   lag(mkt)
============= ============= ========== ===========
  2020-01-02     0.0183088    0.00866    0.00287
  2020-01-03   -0.00352182   -0.00664    0.00866
  2020-01-06    0.00520838    0.00366   -0.00664
  2020-01-07    0.00222056   -0.00184    0.00366
  2020-01-08    0.00387742    0.00476   -0.00184
  2020-01-09    0.00461851    0.00656    0.00476
  2020-01-10    0.00128723   -0.00334    0.00656
  2020-01-13    0.00238753    0.00736   -0.00334
  2020-01-14     0.0054965   -0.00054    0.00736
  2020-01-15   -0.00218664    0.00166   -0.00054
  2020-01-16     0.0122352    0.00886    0.00166
  2020-01-17   -0.00541222    0.00286    0.00886
  2020-01-21    0.00163251   -0.00314    0.00286
  2020-01-22   -0.00905469    0.00086   -0.00314
  2020-01-23    0.00475143    0.00086    0.00086
  2020-01-24    -0.0165515   -0.00964    0.00086
      ⋮             ⋮           ⋮          ⋮
============= ============= ========== ===========
                                  109 rows omitted
```





The dates requested only matter for what is output, internally all data is still stored. In this way, the parameters are largely a "view" into the `mkt_data` object, allowing for quick updating. For example, changing the dates is just:
```julia
AbnormalReturns.update_dates!(orcl_data, Date(2020, 7, 1) .. Date(2020, 7, 5))
```

```
TimelineTable for 10104 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-01-03..2020-12-31
Currently selected dates: 2020-07-01..2020-07-05
============= ============ ======== ===========
        Date          ret      mkt   lag(mkt)
============= ============ ======== ===========
  2020-07-01   0.00398048   0.0041     0.0158
  2020-07-02   0.00810951    0.005     0.0041
============= ============ ======== ===========
```





There are similar functions for `update_id!` (changing the firm ID) and `select!` (changing the columns).

Finally, by default, returned data is not missing. This is to make the regressions faster/easier later, since missing data does not work for those. You can change between allowing or not allowing missing data with `allowmissing` and `dropmissing`:
```julia
no_missings = mkt_data[18428, Date(2019, 3, 28) .. Date(2019, 4, 5)]
```

```
TimelineTable for 18428 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-04-03..2020-12-31
Currently selected dates: 2019-03-28..2019-04-05
============= ============ ======== ======== ========= ==========
        Date          ret      mkt      smb       hml       umd
============= ============ ======== ======== ========= ==========
  2019-04-03       0.0112   0.0028   0.0022   -0.0039    -0.001
  2019-04-04    0.0497538    0.002    0.002    0.0093   -0.0132
  2019-04-05   -0.0413666   0.0053   0.0053   -0.0008   -0.0066
============= ============ ======== ======== ========= ==========
```



```julia
with_missings = allowmissing(no_missings)
```

```
TimelineTable for 18428 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-04-03..2020-12-31
Currently selected dates: 2019-03-28..2019-04-05
============= ============ ========= ========= ========= ==========
        Date          ret       mkt       smb       hml       umd
============= ============ ========= ========= ========= ==========
  2019-03-28      missing   missing   missing   missing   missing
  2019-03-29      missing   missing   missing   missing   missing
  2019-04-01      missing   missing   missing   missing   missing
  2019-04-02      missing   missing   missing   missing   missing
  2019-04-03       0.0112    0.0028    0.0022   -0.0039    -0.001
  2019-04-04    0.0497538     0.002     0.002    0.0093   -0.0132
  2019-04-05   -0.0413666    0.0053    0.0053   -0.0008   -0.0066
============= ============ ========= ========= ========= ==========
```





## Estimating Regressions

The main goal of this package is quickly running regressions for firm events. The example used here is a firm's earnings announcement. Starting with one example, Oracle announced its Q3 2020 earnings on 2020-9-10. Calculating abnormal returns typically follows three steps:
1. Estimate how the firm typically responds to market factors during a control (or estimation) window
2. Use the coefficients from that regression to estimate how the firm should do during the event window
3. Subtract the estimated return from the actual firm return during the event window. Depending on how this difference is aggregated, these are typically buy and hold abnormal returns (bhar) or cumulative abnormla returns (CAR)

First, to create the table for the estimation window, define an estimation window and an event window:
```julia
event_date = Date("2020-09-10")
est_start = advancebdays(mkt_data.calendar, event_date, -120)
est_end = advancebdays(mkt_data.calendar, event_date, -2)
est_start .. est_end
```

```
2020-03-20..2020-09-08
```



```julia
event_start = advancebdays(mkt_data.calendar, event_date, -1)
event_end = advancebdays(mkt_data.calendar, event_date, 1)
event_start .. event_end
```

```
2020-09-09..2020-09-11
```





Next, run the estimation regression (the regression automatically selects the correct columns from the data, so it is not necessary to do that beforehand):
```julia
orcl_data = mkt_data[10104, est_start .. est_end]
rr = quick_reg(orcl_data, @formula(ret ~ mkt + smb + hml))
```

```
Obs: 119, ret ~ -0.0 + 0.77*mkt + -0.302*smb + -0.012*hml, AdjR2: 52.744%
```





Then change the table to the event window:
```julia
AbnormalReturns.update_dates!(orcl_data, event_start .. event_end)
```

```
TimelineTable for 10104 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-01-02..2020-12-31
Currently selected dates: 2020-09-09..2020-09-11
============= ============= ========= ========= ==========
        Date           ret       mkt       smb       hml
============= ============= ========= ========= ==========
  2020-09-09      0.029465    0.0207   -0.0017   -0.0199
  2020-09-10    0.00667254   -0.0164    0.0029   -0.0007
  2020-09-11   -0.00575618   -0.0006   -0.0083     0.009
============= ============= ========= ========= ==========
```





Now it is easy to run some statistics for the event window:
- BHAR based on regression: `bhar(orcl_data, rr)` 0.026330396109940146
- CAR based on regression: `car(orcl_data, rr)` 0.02612077404623775
- BHAR relative to market return: `bhar(orcl_data, "ret", "mkt")` 0.027010625593106186
- Total firm return during event window: `bh_return(orcl_data, "ret")` 0.0303687692811061

It is also easy to calculate some statistics for the estimation window:
- Variance of firm returns (similar equation for standard deviation): `var(rr)` 0.00020283685456709143
- Firm's market beta: `beta(rr)` 0.7698040453871026
- Firm's market alpha: `alpha(rr)` -0.00029218876560305055

## More Data Using DataFramesMeta

While the above works well, abnormal returns are often calculated on thousands or more firm-events. Here, I earnings announcements for about 100 firms from March to November 2020:
```julia
df_events = CSV.File(joinpath("data", "firm_earnings_announcements.csv")) |> DataFrame
```

```
283×2 DataFrame
 Row │ permno  ea
     │ Int64   Date
─────┼────────────────────
   1 │  49373  2020-03-05
   2 │  23660  2020-03-19
   3 │  17144  2020-03-18
   4 │  52708  2020-03-19
   5 │  57665  2020-03-24
   6 │  61621  2020-03-25
   7 │  10104  2020-03-12
   8 │  75154  2020-03-19
  ⋮  │   ⋮         ⋮
 277 │  90373  2020-10-29
 278 │  90386  2020-11-02
 279 │  90993  2020-10-29
 280 │  90880  2020-10-28
 281 │  91668  2020-10-30
 282 │  12449  2020-11-05
 283 │  61399  2020-11-18
          268 rows omitted
```





Using [DataFramesMeta.jl](https://github.com/JuliaData/DataFramesMeta.jl) and the `@chain` macro from [Chain.jl](https://github.com/jkrumbiegel/Chain.jl), the above steps become:

```julia
@chain df_events begin
    @rtransform(
        :est_start = advancebdays(mkt_data.calendar, :ea, -120),
        :est_end = advancebdays(mkt_data.calendar, :ea, -2),
        :event_start = advancebdays(mkt_data.calendar, :ea, -1),
        :event_end = advancebdays(mkt_data.calendar, :ea, 1),
    )
    @rtransform(:reg = quick_reg(mkt_data[:permno, :est_start .. :est_end], @formula(ret ~ mkt + smb + hml)))
    @rtransform(
        :bhar_reg = bhar(mkt_data[:permno, :event_start .. :event_end], :reg),
        :bhar_simple = bhar(mkt_data[:permno, :event_start .. :event_end], "ret", "mkt"),
        :std = std(:reg),
        :total_ret = bh_return(mkt_data[:permno, :event_start .. :event_end], "ret"),
    )
    select(Not([:est_start, :est_end, :event_start, :event_end, :reg]))
end
```

```
283×6 DataFrame
 Row │ permno  ea          bhar_reg     bhar_simple  std         total_ret
     │ Int64   Date        Float64?     Float64      Float64?    Float64
─────┼─────────────────────────────────────────────────────────────────────
───
   1 │  49373  2020-03-05  -0.0114753   -0.0368966   0.0120205   -0.0494716
   2 │  23660  2020-03-19  -0.0761877   -0.0799651   0.0100327   -0.16273
   3 │  17144  2020-03-18   0.00584866  -0.00749583  0.0122784    0.0056200
7
   4 │  52708  2020-03-19   0.0449673    0.057314    0.023234    -0.0254504
   5 │  57665  2020-03-24   0.105811     0.0931054   0.0105352    0.171386
   6 │  61621  2020-03-25   0.125142     0.129935    0.0132682    0.303036
   7 │  10104  2020-03-12   0.0295405    0.0515017   0.00903142  -0.01338
   8 │  75154  2020-03-19   0.15375      0.0269029   0.0312084   -0.0558615
  ⋮  │   ⋮         ⋮            ⋮            ⋮           ⋮            ⋮
 277 │  90373  2020-10-29  -0.00714584  -0.00486512  0.016151    -0.0422143
 278 │  90386  2020-11-02  -0.120063    -0.0769287   0.0256068   -0.0606557
 279 │  90993  2020-10-29  -0.00087777   0.00495886  0.0109232   -0.0323903
 280 │  90880  2020-10-28   0.00943539  -0.00174027  0.0109558   -0.0271723
 281 │  91668  2020-10-30   0.0290468    0.0200758   0.017888     0.0283726
 282 │  12449  2020-11-05  -0.0402344   -0.0572696   0.0183221   -0.0128859
 283 │  61399  2020-11-18  -0.0761987   -0.0706244   0.0134081   -0.0759727
                                                              268 rows omit
ted
```





## Vectorizing the Data

While the above works, and is reasonably fast (Doing a test on 1 million regressions takes about 65 seconds on a Ryzen 5 3600), faster is better.

In particular, this process can be very fast if all the predictor variables in a regression are based on the market data (as they are above). This allows a single RHS matrix to be built and to just select the necessary rows from that matrix. Further, while building the table for each firm is generally fast, if it is only necessary to build it once that is preferred.

To make this possible, request a vector of firm IDs and date starts and ends:
```julia
est_starts = advancebdays.(mkt_data.calendar, df_events.ea, -120)
est_ends = advancebdays.(mkt_data.calendar, df_events.ea, -2)
vec_data = mkt_data[df_events.permno, est_starts .. est_ends]
```

```
Iterable set of TimelineTable with 96 unique firms
and a total number of iterations of 283 and a parent of
MarketData with ID type Int64 with 99 unique firms
MarketCalendar: 2019-01-02 .. 2021-08-31 with 672 business days
Market Columns: mkt, smb, hml, umd
Firm Columns: ret
```





This object can be similarly passed to the above functions, just like a firm level table. The function will iterate through the data and return a vector of results.

However, the above is rather ugly and is far less flexible (no IDs to update, etc.). A more practical way to use this is to continue using the `@chain` macro:
```julia
@chain df_events begin
    @rtransform(
        :est_start = advancebdays(mkt_data.calendar, :ea, -120),
        :est_end = advancebdays(mkt_data.calendar, :ea, -2),
        :event_start = advancebdays(mkt_data.calendar, :ea, -1),
        :event_end = advancebdays(mkt_data.calendar, :ea, 1),
    )
    @transform(:reg = quick_reg(mkt_data[:permno, :est_start .. :est_end], @formula(ret ~ mkt + smb + hml)))
    @transform(
        :bhar_reg = bhar(mkt_data[:permno, :event_start .. :event_end], :reg),
        :bhar_simple = bhar(mkt_data[:permno, :event_start .. :event_end], "ret", "mkt"),
    )
    @rtransform(
        :std = std(:reg),
        :total_ret = bh_return(mkt_data[:permno, :event_start .. :event_end], "ret"),
    )
    select(Not([:est_start, :est_end, :event_start, :event_end, :reg]))
end
```

```
283×6 DataFrame
 Row │ permno  ea          bhar_reg     bhar_simple  std         total_ret
     │ Int64   Date        Float64?     Float64      Float64?    Float64
─────┼─────────────────────────────────────────────────────────────────────
───
   1 │  49373  2020-03-05  -0.0114753   -0.0368966   0.0120205   -0.0494716
   2 │  23660  2020-03-19  -0.0761877   -0.0799651   0.0100327   -0.16273
   3 │  17144  2020-03-18   0.00584866  -0.00749583  0.0122784    0.0056200
7
   4 │  52708  2020-03-19   0.0449673    0.057314    0.023234    -0.0254504
   5 │  57665  2020-03-24   0.105811     0.0931054   0.0105352    0.171386
   6 │  61621  2020-03-25   0.125142     0.129935    0.0132682    0.303036
   7 │  10104  2020-03-12   0.0295405    0.0515017   0.00903142  -0.01338
   8 │  75154  2020-03-19   0.15375      0.0269029   0.0312084   -0.0558615
  ⋮  │   ⋮         ⋮            ⋮            ⋮           ⋮            ⋮
 277 │  90373  2020-10-29  -0.00714584  -0.00486512  0.016151    -0.0422143
 278 │  90386  2020-11-02  -0.120063    -0.0769287   0.0256068   -0.0606557
 279 │  90993  2020-10-29  -0.00087777   0.00495886  0.0109232   -0.0323903
 280 │  90880  2020-10-28   0.00943539  -0.00174027  0.0109558   -0.0271723
 281 │  91668  2020-10-30   0.0290468    0.0200758   0.017888     0.0283726
 282 │  12449  2020-11-05  -0.0402344   -0.0572696   0.0183221   -0.0128859
 283 │  61399  2020-11-18  -0.0761987   -0.0706244   0.0134081   -0.0759727
                                                              268 rows omit
ted
```




Notice that the only difference between these two `@chain` macros is that this one uses `@transform` instead of `@rtransform`. This sends the entire column vector to the function, and allows for much faster overall results. Those same 1 million regressions now takes just 3 seconds on the same computer.

## Lag and Lead Operators

Sometimes, you might want to include a lag or lead variable in your regression. This package is designed to handle those cases. For example:

```julia
quick_reg(mkt_data[10104], @formula(ret ~ lag(mkt, 2) + lag(mkt) + mkt + lead(mkt)))
```

```
Obs: 503, ret ~ 0.0 + -0.166*mkt_lag2 + -0.178*mkt_lag1 + 0.894*mkt + -0.08
5*mkt_lead1, AdjR2: 61.008%
```





Note that these lag and lead operators will get data from outside the requested period to prevent the sample size from shrinking. For example:
```julia
mkt_data[10104, Date(2020) .. Date(2020, 1, 10), [:ret, :mkt, lag(:mkt)]]
```

```
TimelineTable for 10104 and available columns mkt,smb,hml,umd,ret
Maximum dates given current column selection: 2019-01-03..2020-12-31
Currently selected dates: 2020-01-01..2020-01-10
============= ============= ========== ===========
        Date           ret        mkt   lag(mkt)
============= ============= ========== ===========
  2020-01-02     0.0183088    0.00866    0.00287
  2020-01-03   -0.00352182   -0.00664    0.00866
  2020-01-06    0.00520838    0.00366   -0.00664
  2020-01-07    0.00222056   -0.00184    0.00366
  2020-01-08    0.00387742    0.00476   -0.00184
  2020-01-09    0.00461851    0.00656    0.00476
  2020-01-10    0.00128723   -0.00334    0.00656
============= ============= ========== ===========
```




The first value under "lag(mkt)" is from 2019-12-31, but this makes sure that a regression or abnormal return has all the data necessary.