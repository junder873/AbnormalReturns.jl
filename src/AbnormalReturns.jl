module AbnormalReturns

using Tables
using LinearAlgebra
using StatsBase
using Reexport
using StatsModels
using Statistics
using Dates
using DataFrames
using DataFrames: Index
@reexport using BusinessDays

##############################################################################
##
## Exported methods and types
##
##############################################################################
    
# types and functions for fast CAR calculations
export TimelineData, FirmData, car, alpha, beta,
    MarketData, get_firm_data, get_market_data,
    get_firm_market_data, BasicReg, cache_reg,
    bh_return, bhar, clear_firm_cached_data!,
    firm_in_cache, CrspMarketCalendar

export getindex, MatrixTable

# From Statistics
export var, std

# From StatsBase
export coef, coefnames, responsename, nobs, dof_residual,
    r2, adjr2, islinear, deviance, rss, predict

##############################################################################
##
## Load files
##
##############################################################################

include("marketCalendar.jl")
include("timelineDataCache.jl")
include("fastRegression.jl")

end