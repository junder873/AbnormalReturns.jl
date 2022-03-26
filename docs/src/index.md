
# AbnormalReturns Documentation

This package provides functionality for running firm-specific regressions commonly used to calculate abnormal stock returns (actual stock return minus a benchmark). These are common in event studies in finance and economics and often require running a large number of regressions.

Most of the documentation is currently in the [example](@ref Example).

## Motivation

A common firm event is an earnings announcement (when a firm announces bottom line earnings for a quarter). Since 1990, There are over 600,000 of these announcements for firms in the United States. In some academic papers, it is further expected that multiple methods are tested per event (to ensure robustness), therefore, a paper might run several million regressions. Under most common methods, this can take a long time. For example, a common SAS macro to run these regressions took over 30 minutes to complete, which does not include downloading the initial data. This package, on its first run (so including compilation), was able to download the necessary data and run the regression in under 5 total minutes (4 minutes for downloading the data, 30 seconds for arranging the data, 20 for running the regressions). A second run took under 3 minutes (still downloading the data).

Overall, this package can run 1 million regressions in under 3 seconds, which, as far as I am aware, make this package easily the fastest method to calculate such regressions.