
# AbnormalReturns API

This package reexports [StatsModels.jl](https://raw.githubusercontent.com/JuliaStats/StatsModels.jl) and [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl). See those packages and their respective documentation for methods that they export.

Notable methods from [StatsModels.jl](https://raw.githubusercontent.com/JuliaStats/StatsModels.jl):
- `@formula`

Notable methods from [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl):
- `advancebdays`
- `bdayscount`
- `tobday`

## Setting Up Data

```@docs
MarketData
AbnormalReturns.all_unique_obs
```

## Regression Related Methods

```@docs
quick_reg
BasicReg
alpha
beta
var
std
```

## Calculation Functions

```@docs
bhar
car
bh_return
```