module AbnormalReturns

using Tables
using LinearAlgebra
using StatsBase
using Reexport
using Statistics
using Dates
using DataFrames
using IntervalSets: ClosedInterval, Ellipsis, (..)
using SparseArrays
@reexport using BusinessDays
@reexport using StatsModels
using ShiftedArrays
using PrettyTables
using StaticArrays

##############################################################################
##
## Exported methods and types
##
##############################################################################
    
# types and functions for fast CAR calculations
export MarketData, TimelineTable, car, alpha, beta,
    BasicReg, quick_reg, DataVector,
    bh_return, bhar, MarketCalendar, TimelineColumn

export getindex, values, names, istable, columnaccess, columns,
    getcolumn, columnnames

# From Statistics
export var, std

# From StatsBase
export coef, coefnames, responsename, nobs, dof_residual,
    r2, adjr2, islinear, deviance, rss, predict

export getindex

export dropmissing, select!

export ClosedInterval, ..

export lag, lead

##############################################################################
##
## Load files
##
##############################################################################

include("marketCalendar.jl")
include("dictIndex.jl")
include("timelineData.jl")
include("utils.jl")
include("iterateTimelineTable.jl")
include("fastRegression.jl")
include("calcFunctions.jl")
include("statsModelsModelcols.jl")

end