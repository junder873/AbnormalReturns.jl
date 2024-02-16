# AbnormalReturns Data Structure

The key to the performance in this package is the underlying data structure. These rely on a combination of [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl) and [Tables.jl](https://github.com/JuliaData/Tables.jl) to provide fast access to slices of data based on dates.

## DataVector

```@docs
AbnormalReturns.DataVector
```

These structures provide strongly typed data that is easy to slice based on a range of dates. The data is always stored as `Float64`, even though it accepts elements of type `Missing`. In storing the data, `Missing` values are converted to `0.0`, and the `missing_bdays` is a `SparseVector` that is `true` when that value is missing. `dates` are the minimum and maximum dates for the data.

Data in a `DataVector` is stored in an `OffsetVector` from [OffsetArrays.jl](https://github.com/JuliaArrays/OffsetArrays.jl), but this data is rarely accessed directly.

## MarketData

```julia
    struct MarketData{T}
        calendar::MarketCalendar
        marketdata::Dict{Symbol, DataVector} # column names as symbols
        firmdata::Dict{T, Dict{Symbol, DataVector}} # data stored by firm id and then by column name as symbol
    end
```

This struct is made up of a set of `DataVector`. The main purpose of this is to provide efficient storage of the underlying data along with the calendar so that dates are easily translated to access the underlying data.

## FixedTable

```@docs
FixedTable
```

## Market Calendar

```@docs
AbnormalReturns.MarketCalendar
```