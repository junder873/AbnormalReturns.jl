using DataFrames, CSV, DataFramesMeta, Dates, BenchmarkTools, Cthulhu, SparseArrays, StaticArrays, LinearAlgebra
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##
@time data = MarketData(df_mkt, df_firm; id_col=:firm_id, valuecols_firms=[:ret])
# First run R5 3600: 20.229416 seconds (34.98 M allocations: 13.211 GiB, 12.81% gc time, 64.56% compilation time)
# Second run R5 3600: 8.131775 seconds (681.52 k allocations: 11.431 GiB, 26.59% gc time, 0.06% compilation time)
# First run R7 5700X: 13.649963 seconds (35.34 M allocations: 13.235 GiB, 14.57% gc time, 64.03% compilation time)
# Second run R7 5700X: 5.804785 seconds (587.93 k allocations: 11.432 GiB, 29.49% gc time)
# First run i7 6700: 30.367487 seconds (29.85 M allocations: 12.966 GiB, 40.87% gc time, 38.86% compilation time)
# Second run i7 6700: 10.302273 seconds (688.03 k allocations: 11.436 GiB, 22.65% gc time)

##
df_firm = nothing
df_mkt = nothing
GC.gc()

##
@time df_temp = @chain df_events begin
    @transform(:reg_mkt = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ smb + mkt)),)
    @transform(:reg_ffm = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
    @transform(
        :bhar_mkt = bhar(data[:firm_id, :est_window_start .. :est_window_end, @formula(ret ~ smb + mkt)], :reg_mkt),
        :bhar_ffm = bhar(data[:firm_id, :est_window_start .. :est_window_end, @formula(ret ~ mkt + smb + hml + umd)], :reg_ffm),
        # :bhar_simple = bhar(data[:firm_id, :est_window_start .. :est_window_end], ^(:ret), ^(:mkt)),
        # :car_simple = car(data[:firm_id, :est_window_start .. :est_window_end], ^(:ret), ^(:mkt))
    )
    @rtransform(
        :var_mkt = var(:reg_mkt),
        :var_ffm = var(:reg_ffm),
    )
    @transform(
        #:alpha = alpha(:reg_mkt),
        :beta = beta(:reg_mkt),
    )
end
# First run R5 3600: 12.886998 seconds (35.27 M allocations: 4.308 GiB, 13.43% gc time, 70.63% compilation time)
# Second run R5 3600: 3.537313 seconds (12.73 M allocations: 3.092 GiB, 24.15% gc time, 2.50% compilation time)
# First run R7 5700X: 7.630244 seconds (31.75 M allocations: 3.091 GiB, 14.36% gc time, 79.25% compilation time)
# Second run R7 5700X: 1.492441 seconds (12.18 M allocations: 2.084 GiB, 23.40% gc time, 4.99% compilation time)
# First run i7 6700: 17.520028 seconds (26.97 M allocations: 4.145 GiB, 28.89% gc time, 77.12% compilation time)
# Second run i7 6700: 5.772432 seconds (4.73 M allocations: 2.944 GiB, 44.27% gc time, 1.59% compilation time)

##


@time @chain df_events[1:1000000, :] begin
    @rtransform(:reg = BasicReg(data[:firm_id, :est_window_start .. :est_window_end, @formula(ret ~ mkt + smb + hml + umd)], @formula(ret ~ mkt + smb + hml + umd)),)
end
# Run R5 3600: 55.141968 seconds (589.18 M allocations: 59.681 GiB, 12.77% gc time, 0.19% compilation time)

##

@time @chain df_events[1:1000000, :] begin
    @transform(:reg = BasicReg.(data[:firm_id, :est_window_start .. :est_window_end, @formula(ret ~ mkt + smb + hml + umd)], Ref(@formula(ret ~ mkt + smb + hml + umd))),)
end
##
x = NTuple{2, Vector{Int}}(([1, 2, 3], [4, 5, 6]))
##
const to2 = TimerOutput()

##
function test_reg(
    data::TimelineTable{false},# only works if no missing data for multiplication
    f::FormulaTerm{L, R};
    minobs::V=0.8,
    save_residuals::Bool=false,
    sch=sch
) where {L, R, V<:Real}

    @timeit to2 "start" if !isa(f.rhs[1], ConstantTerm)#!StatsModels.omitsintercept(f) && !StatsModels.hasintercept(f)
        f = FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs)
        # return test_reg(
        #     data,
        #     FormulaTerm(f.lhs, InterceptTerm{true}() + f.rhs);
        #     minobs,
        #     save_residuals
        # )
    end
    @timeit to2 "schema" sch = apply_schema(f, sch)

    @timeit to2 "termvars" x = AbnormalReturns.internal_termvars(sch)
    @timeit to2 "select" select!(data, x)
    
    
    @timeit to2 "modelcols response" resp = AbnormalReturns.datavector_modelcols(sch.lhs, data)
    @timeit to2 "modelcols pred" pred = AbnormalReturns.datavector_modelcols(sch.rhs, data)

    
    @timeit to2 "views" if nnz(AbnormalReturns.data_missing_bdays(data)) == 0
        y = view(resp, AbnormalReturns.norm_dates(data))
        x = AbnormalReturns.FixedWidthMatrix(pred, AbnormalReturns.norm_dates(data))
    else
        new_mssngs = AbnormalReturns.get_missing_bdays(data, AbnormalReturns.norm_dates(data))
        y = view(resp, AbnormalReturns.norm_dates(data), new_mssngs)
        x = AbnormalReturns.FixedWidthMatrix(pred, AbnormalReturns.norm_dates(data), new_mssngs)
    end
    
    yname = coefnames(sch.lhs)
    @timeit to2 "names" xnames2 = SVector{width(sch.rhs)}(coefnames(sch.rhs))
    
    @timeit to2 "reg" BasicReg(
        y,
        x,
        yname,
        xnames2,
        f;
        save_residuals,
        minobs=AbnormalReturns.adjust_minobs(minobs, AbnormalReturns.calendar(data), AbnormalReturns.norm_dates(data))
    )
    
end
temp_f = @formula(ret ~  mkt + smb + 1 + hml + umd)
test_reg(data[4400, Date(2003, 10, 7) .. Date(2004, 9, 22), cols], temp_f)

##
temp_f = term(:f) ~ term(:smb) + ConstantTerm(1)
apply_schema(temp_f, schema(temp_f, data))
##

@time @chain df_events[1:1000000, :] begin
    @rtransform(:reg = test_reg(data[:firm_id, :est_window_start .. :est_window_end, cols], @formula(ret ~ mkt + smb + hml + umd)),)
end

##
sch = schema(temp_f, data)
##
temp_f = @formula(ret ~ 0 + mkt + smb)
sch = apply_schema(temp_f, schema(temp_f, data))
isa(sch.rhs.terms, Tuple{InterceptTerm{false}, Vararg{ContinuousTerm{Float64}}})
##

@time @chain df_events[1:1000000, :] begin
    @transform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end], @formula(ret ~ mkt + smb + hml + umd)),)
end
# Run R5 3600: 0.655350 seconds (1.53 M allocations: 850.312 MiB, 5.46% compilation time)
# Run i7 6700: 1.060585 seconds (1.53 M allocations: 850.309 MiB, 3.80% compilation time)
##
@time @chain df_events[1:1000000, :] begin
    @transform(:bhar = bhar(data[:firm_id, :est_window_start .. :est_window_end]),)
end
##
quick_reg(data[1, Date(2017) .. Date(2018), cols], @formula(ret ~ mkt + smb + hml + umd))
##


@benchmark $data[i, TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])] setup=(i=rand(1:10000))
##


# i = 151465
# ids = df_events.firm_id[1:100]
# date_mins = df_events.est_window_start[1:100]
# date_maxs = df_events.est_window_end[1:100]
# f = @formula(ret ~ mkt + smb + hml + umd)
# cols = AbnormalReturns.internal_termvars(f)
# data_dict = AbnormalReturns.construct_id_dict(ids)
# sch = apply_schema(f, schema(f, data))
# cache = AbnormalReturns.create_pred_matrix(data[1], sch)
# out = fill(BasicReg(0, f), length(ids))

# int_data = data[1]



b = @benchmarkable AbnormalReturns.vector_reg(
    $data,
    $ids,
    $date_mins,
    $date_maxs,
    $f,
)
b.params.evals=100
b.params.seconds=100
run(b)
##
