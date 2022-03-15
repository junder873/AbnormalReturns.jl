using DataFrames, CSV, DataFramesMeta, Dates, BenchmarkTools, Cthulhu, Traceur
using Revise
using AbnormalReturns

##

df_firm = CSV.File(joinpath("data", "firm_ret.csv")) |> DataFrame
df_mkt = CSV.File(joinpath("data", "mkt_ret.csv")) |> DataFrame
df_events = CSV.File(joinpath("data", "event_dates.csv")) |> DataFrame

##

@time data = MarketData(df_mkt, df_firm; id_col=:firm_id)
# First run: 38.816633 seconds (20.05 M allocations: 12.740 GiB, 44.66% gc time, 20.50% compilation time)
# Second run: 18.013269 seconds (629.12 k allocations: 11.719 GiB, 16.72% gc time)

##
@time df_test = @chain df_events begin
    @transform(
        :reg_mkt = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt)),
        :reg_ffm = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt + smb + hml + umd)),
    )
    @transform(
        :bhar_mkt = AbnormalReturns.bhar(data, :firm_id, :event_window_start, :event_window_end, :reg_mkt),
        :bhar_ffm = AbnormalReturns.bhar(data, :firm_id, :event_window_start, :event_window_end, :reg_ffm),
    )
    @rtransform(
        :var_mkt = var(:reg_mkt),
        :var_ffm = var(:reg_ffm),
        :alpha = alpha(:reg_mkt),
        :beta = beta(:reg_mkt),
    )
end
# 601.901257 seconds (412.12 M allocations: 52.929 GiB, 95.80% gc time, 0.24% compilation time)
# 65.207394 seconds (220.15 M allocations: 42.688 GiB, 59.49% gc time, 0.12% compilation time)

##

@time @chain df_events[1:100000, :] begin
    @transform(:reg = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ 1 + mkt + smb + hml + umd)),)
end

##

x = temp_data[TimelineColumn(:ret)]
@code_warntype x[Date(2018) .. Date(2019)]

##
@benchmark $x[Date(2018) .. Date(2019)]

##
@benchmark view(x, Date(2018) .. Date(2019))
##
x = randn(10000)
@benchmark view($x, y:y+250) setup=(y=rand(1:9000))
##

@code_warntype temp_data[:, TimelineColumn(:ret)]

##
temp_data = data[1]
AbnormalReturns.update_dates!(temp_data, Date(2018) .. Date(2019))
@benchmark $temp_data[:, TimelineColumn(:ret)]
##

@time @chain df_events begin
    @transform(:reg = AbnormalReturns.vector_reg(data, :firm_id, :est_window_start, :est_window_end, @formula(ret ~ mkt + smb)),)
end
##
cols = TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])
@time @chain df_events[1:1000000, :] begin
    @rtransform(:reg = quick_reg(data[:firm_id, :est_window_start .. :est_window_end, cols], @formula(ret ~ mkt + smb + hml + umd)),)
end

##
quick_reg(data[1, Date(2017) .. Date(2018), cols], @formula(ret ~ mkt + smb + hml + umd))
##


##

@benchmark $temp_data[:, :ret]
##
temp_data[:, TimelineColumn(:ret)]
##

@benchmark $data[i, TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])] setup=(i=rand(1:10000))
##

@time AbnormalReturns.vector_reg(data, df_events.firm_id, df_events.est_window_start, df_events.est_window_end, @formula(ret ~ mkt + smb + hml + umd); minobs=100)

##
@code_warntype data[50, TimelineColumn.([:ret, :mkt, :smb, :hml, :umd])]

##

@code_warntype temp_data[Date(2018) .. Date(2019)]

##
function f_test(parent_data, ids, date_starts, date_ends, f)
    
    # out = AbnormalReturns.TimelineTable[]
    # sizehint!(out, length(ids))
    cols = AbnormalReturns.internal_termvars(f)
    data = parent_data[ids[1], date_starts[1] .. date_ends[1], cols]
    for i in 1:length(ids)
        AbnormalReturns.update_data!(data; new_id=ids[i], new_dates=date_starts[i] .. date_ends[i])
        #push!(out, parent_data[ids[i], date_starts[i] .. date_ends[i], cols])
    end 
    #out
    data
end

function f_test2(parent_data, ids, date_starts, date_ends, f)
    ids = unique(ids)
    cols = AbnormalReturns.internal_termvars(f)

    Dict(ids .=> [parent_data[id, Date(2018) .. Date(2019), cols] for id in ids])
end

@time f_test2(data, df_events.firm_id, df_events.est_window_start, df_events.est_window_end, @formula(ret ~ mkt + smb + hml + umd))


##

data

##

temp_data
##

@code_warntype AbnormalReturns.update_data!(temp_data, new_id=50, new_dates = Date(2018) .. Date(2019))

##
@code_warntype AbnormalReturns.dates_min_max(Date(2018) .. Date(2019), [AbnormalReturns.data_dates(data[50, col]) for col in names(temp_data)]...)
##

function f_test(pred, resp; f=f, save_residuals=false, yname="h", xnames=coefnames(sch.rhs))
    coef = cholesky!(Symmetric(pred' * pred)) \ (pred' * resp)
    # resid = resp .- pred * coef
    # rss = sum(abs2, resid)
    # tss = sum(abs2, (resp .- mean(resp)))
    # #yname, xnames = coefnames(sch)
    # BasicReg(
    #     length(resp),
    #     f,
    #     coef,
    #     xnames,
    #     yname,
    #     tss,
    #     rss,
    #     save_residuals ? resid : nothing
    # )
end
i = 151465
ids = df_events.firm_id
date_mins = df_events.est_window_start
date_maxs = df_events.est_window_end
f = @formula(ret ~ mkt + smb + hml + umd)
cols = AbnormalReturns.internal_termvars(f)
data_dict = AbnormalReturns.construct_id_dict(ids)
sch = apply_schema(f, schema(f, data))
cache = AbnormalReturns.create_pred_matrix(data[1], sch)
out = fill(BasicReg(0, f), length(ids))

int_data = data[1]

@time AbnormalReturns.fill_vector_reg(
    data,
    date_mins,
    date_maxs,
    out,
    cache,
    sch,
    data_dict,
    cols,
    f,
    minobs=0.8,
    save_residuals=false
)
##


# @code_warntype BasicReg(
#             data_dict[ids[i]][date_mins[i] .. date_maxs[i]],
#             cache,
#             sch,
#             f;
#             minobs=0.8
#         )

#@benchmark f_test($pred, $resp)
#@benchmark BasicReg($data, $cache, $sch, $f)
#@benchmark modelcols($sch.lhs, $data)
#@benchmark $cache[$data.dates, nothing]
#@benchmark coefnames($sch)
#@benchmark cholesky!(Symmetric($pred' * $pred)) \ ($pred' * $resp)
##
function f_test2(sch::FormulaTerm{L2, R2}) where {L2, R2}
    # println(L2)
    # println(typeof(L2))
    println(sch.lhs)
    println(typeof(sch.lhs))
end
@code_warntype f_test2(sch)

##
function f_test3(sch::FormulaTerm)
    # println(L2)
    # println(typeof(L2))
    println(sch.lhs)
    println(typeof(sch.lhs))
end
@code_warntype f_test3(sch)
##
x = data[:, [:mkt, :smb, :hml, :umd]]
x_h = x'
y = data[:, :ret]
##

@benchmark quick_reg($data[i], $f) setup=(i=rand(1:10000))
##
@code_warntype data.id
##
f(x) = x.id
@code_warntype f(data)
##
@benchmark AbnormalReturns.update_request!($data; id) setup=(id=rand(1:10000))
##

function rand_data_test(data::MarketData{T}, id::T) where {T}
    AbnormalReturns.update_request!(data; id, dates=Date(2011)..Date(2012))
    data[:, :mkt] |> sum
end

@benchmark rand_data_test($data, id) setup=(id=rand(1:10000))

##
@benchmark AbnormalReturns.dates_min_max(AbnormalReturns.all_dates($data))

##
@code_warntype AbnormalReturns.update_request!(data; id=100, dates = Date(2016) .. Date(2018))

##
function get_all_items2(parent_data::MarketData{T}, ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    # int_data = MarketData(
    #     getfield(parent_data, :calendar),
    #     getfield(parent_data, :marketdata),
    #     getfield(parent_data, :firmdata),
    #     getfield(parent_data, :request)
    # )
    out = zeros(length(ids))
    for (i, id) in enumerate(ids)
        #AbnormalReturns.update_request!(parent_data; id, dates= date_starts[i] .. date_ends[i])
        parent_data.request.id = id
        parent_data.request.req_dates = date_starts[i] .. date_ends[i]
        parent_data.request.dates = AbnormalReturns.dates_min_max(date_starts[i] .. date_ends[i], AbnormalReturns.all_dates(parent_data)...)
        AbnormalReturns.calculate_missing_bdays!(data)

        #data = parent_data[ids[i], date_starts[i] .. date_ends[i], [:ret]]
        #x = parent_data[:, :mkt]
        #out[i] = sum(x)
    end
    out
end

#@code_warntype get_all_items2(data,  df_events.firm_id[1:10], df_events.est_window_start[1:10], df_events.est_window_end[1:10])
@time get_all_items2(data,  df_events.firm_id, df_events.est_window_start, df_events.est_window_end)

##

function get_ind_column(parent_data::MarketData{T}, ids::Vector{T}, date_starts::Vector{Date}, date_ends::Vector{Date}) where {T}
    out = zeros(length(ids))
    for i in 1:length(ids)
        x = getfield(parent_data, :firmdata)[ids[i]].ret
        r = AbnormalReturns.date_range(parent_data.calendar, x.dates, date_starts[i] .. date_ends[i])
        out[i] = sum(@view x.data[r])
    end
    out
end

@time get_ind_column(data,  df_events.firm_id, df_events.est_window_start, df_events.est_window_end)

##
@code_warntype get_ind_column(data,  df_events.firm_id[1:10], df_events.est_window_start[1:10], df_events.est_window_end[1:10])

##

@code_warntype data[500, Date(2015) .. Date(2016), [TimelineColumn(:ret)]]