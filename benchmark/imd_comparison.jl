using Dates, BenchmarkTools, InMemoryDatasets
using Revise
using AbnormalReturns
using DLMReader
using LinearAlgebra, Random




ds_firm = filereader(joinpath("data", "firm_ret.csv"), types = Dict(2=>Date)) 
ds_mkt = filereader(joinpath("data", "mkt_ret.csv"), types = Dict(1=>Date));
ds_events=filereader(joinpath("data","event_dates.csv"),types = Dict(2:6 .=>Date)) |> unique;


@time IMD.leftjoin!(ds_firm, ds_mkt[!, [:date]], on = :date, obs_id = [false, true])
# 9.554551 seconds (18.55 M allocations: 3.019 GiB, 5.20% gc time, 59.21% compilation time)
# 1 threads second run: 2.472221 seconds (81.86 k allocations: 1.508 GiB, 9.02% gc time, 1.41% compilation time)
# 4 threads:  7.670734 seconds (24.06 M allocations: 3.323 GiB, 7.50% gc time, 83.31% compilation time)
# 4 threads second run: 1.172847 seconds (81.88 k allocations: 1.508 GiB, 20.61% gc time, 3.67% compilation time)
# 6 threads: 7.470740 seconds (24.06 M allocations: 3.323 GiB, 7.71% gc time, 85.09% compilation time)
# 6 threads second run: 0.978876 seconds (81.89 k allocations: 1.508 GiB, 17.10% gc time, 4.15% compilation time)




m = disallowmissing(ds_mkt[!, 2:5]|>Matrix)



@time ds_event_joined=IMD.innerjoin(
    ds_firm,
    view(ds_events,:,[:firm_id,:event_date,:est_window_start, :est_window_end]), 
    on=[:firm_id => :firm_id, :date => (:est_window_start, :est_window_end)],
    );
# 115.654700 seconds (14.67 M allocations: 17.010 GiB, 6.81% gc time, 6.50% compilation time)
# run with 4 threads:  58.524648 seconds (7.51 M allocations: 16.613 GiB, 12.74% gc time, 6.83% compilation time)
# run with 6 threads: 32.617801 seconds (14.67 M allocations: 17.010 GiB, 6.06% gc time, 17.52% compilation time)
# second run with 6 threads: 35.282279 seconds (11.14 k allocations: 16.205 GiB, 11.54% gc time)

@time groupby!(ds_event_joined, [:firm_id, :event_date],stable=false)
#  79.080849 seconds (2.73 M allocations: 19.911 GiB, 4.84% gc time, 1.64% compilation time)
# run with 6 threads: 35.726548 seconds (1.64 M allocations: 19.862 GiB, 9.33% gc time, 5.06% compilation time)
# second run with 6 threads: 0.000083 seconds (16 allocations: 896 bytes)

function simple_reg(obs_id,y; m = m
    )::Vector{Float64}

    pred = view(m,obs_id,:)
    resp = y
    cholesky!(Symmetric(pred' * pred)) \ (pred' * resp)
end;

@time IMD.combine(ds_event_joined, ( :obs_id_right,:ret) => simple_reg => :re_output )
# 17.010435 seconds (30.98 M allocations: 3.499 GiB, 9.47% gc time, 22.09% compilation time)
# run with 4 threads:  10.071312 seconds (14.77 M allocations: 965.885 MiB, 6.66% gc time, 15.58% compilation time)
# run with 6 threads: 8.359918 seconds (19.96 M allocations: 1.222 GiB, 11.31% gc time, 47.56% compilation time)
# second run with 6 threads: 3.187782 seconds (9.96 M allocations: 741.199 MiB, 10.31% gc time)




@time ds_event_joined=IMD.innerjoin(
    ds_firm,
    view(ds_events,:,[:firm_id,:event_date,:est_window_start, :est_window_end]), 
    on=[:firm_id => :firm_id, :date => (:est_window_start, :est_window_end)],
    );


@time groupby!(ds_event_joined, [:firm_id, :event_date],stable=false)


