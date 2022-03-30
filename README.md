[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://junder873.github.io/AbnormalReturns.jl/dev/)
[![Build status](https://github.com/junder873/AbnormalReturns.jl/workflows/CI/badge.svg)](https://github.com/junder873/AbnormalReturns.jl/actions)

# AbnormalReturns.jl

This package is designed to quickly calculate abnormal returns on large datasets by running regressions on slices of dates. In finance and economics research, abnormal returns are common for event studies related to firms to interpret how the stock market perceives the event. For example, if a firm makes an announcement, did the market see that as good news? To what degree (i.e., how big are the returns)?

Calculating abnormal returns typically requires running regressions on a slice of the data (during an estimation window) and using those coefficients to predict what a firm's returns would be during an event window. The difference between the the actual returns and the expected returns is used as a measure of abnormal returns.

## Performance

This package is capable of calculating abnormal returns very quickly. On a Ryzen 5 3600, it can calculate 1 million different regressions on different slices of data in under 3 seconds. In a larger benchmark using simulated data for 1 million firm-events, this package can calculate abnormal returns for all events using two methods (so 2 million total regressions, 2 million estimations and some other statistics) in 9 seconds. See the benchmark folder for more details.

## Acknowledgements

This package would not be possible without the extensive work done in [BusinessDays.jl](https://github.com/JuliaFinance/BusinessDays.jl) and [StatsModels.jl](https://github.com/JuliaStats/StatsModels.jl).

## Disclaimer

While this package does have testing, it is in beta. Methods might change and there could be errors.
