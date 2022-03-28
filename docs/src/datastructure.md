# AbnormalReturns Data Structure

The key to the performance in this package is the underlying data structure. These rely on a combination of [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl) and [Tables.jl](https://github.com/JuliaData/Tables.jl) to provide fast access to slices of data based on dates.

## DataVector (and DataMatrix)

```@docs
DataVector
DataMatrix
```

These structures provide strongly typed data that is easy to slice based on a range of dates. The data is always stored as `Float64`, even though it accepts elements of type `Missing`. In storing the data, `Missing` values are converted to `0.0`, and the `missing_bdays` is a `SparseVector` that is `true` when that value is missing. `dates` are the minimum and maximum dates for the data.

Access to a `DataVector` or Matrix is through providing a `ClosedInterval` from [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl):

```julia
data[Date(2018) .. Date(2019)]
```

This will return a Vector or Matrix, depending on type.

!!! Note
    The data returned is always a Vector (or Matrix) of `Float64`, even if there was missing data in the original Vector. To truly get back to the original, you must combine the `missing_bdays` with the Vector.

## MarketData

```julia
    struct MarketData{T, MNames, FNames, N1, N2}
        calendar::MarketCalendar
        marketdata::NamedTuple{MNames, NTuple{N1, DataVector}} # column names as symbols
        firmdata::Dict{T, NamedTuple{FNames, NTuple{N2, DataVector}}} # data stored by firm id and then by column name as symbol
    end
```

This struct is made up of a set of `DataVector`. The main purpose of this is to provide efficient storage of the underlying data.

## TimelineTable

```@docs
TimelineTable
```
