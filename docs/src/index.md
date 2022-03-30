
# AbnormalReturns Documentation

This package provides functionality for getting slices of firm and market return data and running firm-specific regressions commonly used to calculate abnormal stock returns (actual stock return minus a benchmark). These are common in event studies in finance and economics and often require running a large number of regressions on different slices of firm and market data.

Most of the documentation is currently in the [example](@ref Example).

## Motivation

When estimating abnormal returns, it is common to estimate how the firm's return typically responds during an estimation window and use those predicted results in an event window:

![](images/event_timeline.PNG)

The exact length of the estimation and event windows varies, but are typically about 150 and 3-5, respectively. The estimation is typically is a linear regression of firm specific return on market-wide factors.

### The Problem

Estimating abnormal returns requires getting two separate slices of data (for the estimation window and event window) for each firm-event. This is relatively trivial for small datasets, but abnormal returns are often calculated for a large number of events. For example, there are over 600,000 firm earnings announcements since 1990.

Generally, creating the dataset is done through a range join (e.g., gather all firm data between the start and end of the estimation window), which is often time consuming and/or creates huge datasets.

### This Package

This package uses a custom data structure to avoid repeating the data. The data structure is built on [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl), making it easy to get a slice of data between two dates. It also implements threaded solutions to make the regressions and aggregation as fast as possible.

In a benchmark on 1 million firm events, it runs all the regressions in under 3 seconds. In a larger benchmark with two different models (so 2 million regressions) and calculating abnormal returns for the events, along with other basic statistics, it takes less than 9 seconds on a Ryzen 5 3600.

## Acknowledgements

This package would not be possible without [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl), which is used for all of the date operations in this package and [StatsModels.jl](https://github.com/JuliaStats/StatsModels.jl), which provides an incredible `@formula` macro and the functionality that comes with that.