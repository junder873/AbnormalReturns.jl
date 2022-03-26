# Benchmarking AbnormalReturns.jl

Since the goal of this package is to calculate abnormal returns on a large number of events, a large dataset is needed. Since I cannot publish actual market data due to copyright, "create_approximate_real.jl" downloads real market data and runs estimates to simulate what some real data would look like. It is not necessary to simulate actual data. The actual data is downloaded from the Wharton Research Database (WRDS) and uses CRSP and data from [Fama-French](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html).

The file "create_data.jl" uses the real correlations and mean coefficients to build a sample of market level returns, firm specific returns for 10,000 returns, and 1 million firm events.

Finally, "run_benchmark.jl" loads the generated data and runs timings. I included the timings on two computers I have access to in that file.