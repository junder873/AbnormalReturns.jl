module AbnormalReturns

using LinearAlgebra
using StatsBase
using Reexport
using Statistics
using Dates
using DataFrames
using Tables
using IntervalSets: ClosedInterval, Ellipsis, (..)
using SparseArrays
@reexport using BusinessDays
@reexport using StatsModels
using StaticArrays
using OffsetArrays

##############################################################################
##
## Exported methods and types
##
##############################################################################
    
# types and functions for fast CAR calculations
export MarketData, FixedTable, car, alpha, beta,
    BasicReg, quick_reg, IterateFixedTable,
    bh_return, bhar, MarketCalendar

export getindex, names

# From Statistics
export var, std

# From StatsBase
export coef, coefnames, responsename, nobs, dof_residual,
    r2, adjr2, islinear, deviance, rss, predict

export ClosedInterval, ..

##############################################################################
##
## Load files
##
##############################################################################

include("marketCalendar.jl")
include("timelineData.jl")
include("calcUtils.jl")
include("iterateTimelineTable.jl")
include("fastRegression.jl")
include("calcFunctions.jl")

end