# AbnormalReturns.jl

This package is designed to calculate abnormal returns on large datasets. In finance and economics research, abnormal returns are common for event studies related to firms to interpret how the stock market perceives the event. For example, if a firm makes an announcement, did the market see that as good news? To what degree (i.e., how big are the returns)?

Calculating abnormal returns typically requires running regressions on a slice of the data (during an estimation window) and using those coefficients to predict what a firm's returns would be during an event window. The difference between the the actual returns and the expected returns is used as a measure of abnormal returns.

## Problem with Other Methods

Calculating these returns over large datasets can be time consuming since there is one regression (and two slices of the data) for each firm event. The only method I am aware of to do this first merges the market data and firm data into one table then uses a group by to calculate the regressions. However, creating the initial table can be costly (if there are 100,000 events and 200 days in the estimation window, that table has 20 million rows), so the most common software is SAS which allows for such size while using advanced statistical methods. This package aims to solve this issue.

## Pakcage Overview

This package provides two main data structures, `MarketData` stores the firm specific and market average data and tracks which days are missing. This means that when requesting data for a specific firm over a number of days, a set of data is quickly returned. After requesting data from `MarketData`, a `TimelineTable` is returned, which is [Tables.jl](https://github.com/JuliaData/Tables.jl) compatible. These slices are especially quick since they are partially lazy and are built on [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl). This package also implements [StatsModels.jl](https://github.com/JuliaStats/StatsModels.jl).

## Example

