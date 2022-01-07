module AbnormalReturns

using Tables
using LinearAlgebra
using StatsBase
using Reexport
using Statistics
using Dates
using DataFrames
using DataFrames: Index
using IntervalSets
@reexport using BusinessDays
@reexport using StatsModels

##############################################################################
##
## Exported methods and types
##
##############################################################################
    
# types and functions for fast CAR calculations
export MarketData, TimelineTable, car, alpha, beta,
    BasicReg, quick_reg, DataVector,
    bh_return, bhar, MarketCalendar

export getindex, values, names, istable, columnaccess, columns,
    getcolumn, columnnames

# From Statistics
export var, std

# From StatsBase
export coef, coefnames, responsename, nobs, dof_residual,
    r2, adjr2, islinear, deviance, rss, predict

export getindex

export dropmissing!, select!

##############################################################################
##
## Load files
##
##############################################################################

include("marketCalendar.jl")
include("timelineData.jl")
include("fastRegression.jl")
include("calcFunctions.jl")

end