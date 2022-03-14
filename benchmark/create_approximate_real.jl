using Distributions, DataFrames, DataFramesMeta, Arrow, Statistics, BusinessDays, Dates, CSV
using Revise
using AbnormalReturns

##

df_ff_data = @chain Arrow.Table(abspath("G:", "My Drive", "python", "tests", "werdsmerger tests", "temp_data", "ff_data.feather")) begin
    DataFrame
    copy
    @rtransform(:mkt = :mktrf + :rf)
end

ff_data_max = maximum(df_ff_data.date)

@time df_crsp_raw = @chain Arrow.Table(abspath("G:", "My Drive", "python", "Share Pledging Analysis", "data", "general csv files", "crsp_daily.feather")) begin
    DataFrame
    select([:permno, :date, :ret, :retx, :shrout_adj])
    @rsubset(:date <= ff_data_max)
    copy
    sort([:permno, :date])
end

data = MarketData(
    df_ff_data,
    df_crsp_raw
)

##

mkt_mat = @chain df_ff_data begin
    select([:mktrf, :smb, :hml, :umd])
    dropmissing
    Matrix
    cov
end
# 4×4 Matrix{Float64}:
#   0.000116248  -1.03659e-5   1.1843e-5   -1.04678e-5
#  -1.03659e-5    3.49562e-5  -2.37018e-6   1.9103e-6
#   1.1843e-5    -2.37018e-6   3.77871e-5  -1.05367e-5
#  -1.04678e-5    1.9103e-6   -1.05367e-5   5.98879e-5


##


@time df_all_regs = @chain df_crsp_raw begin
    select([:permno, :date])
    dropmissing
    @by(
        :permno,
        :date_start = minimum(:date),
        :date_end = maximum(:date)
    )
    @transform(:reg = AbnormalReturns.vector_reg(data, :permno, :date_start, :date_end, @formula(ret ~ mktrf + smb + hml + umd); minobs=100, save_residuals=true))
end


##

function winsorize(vect::AbstractArray, low_per::AbstractFloat=0.01, high_per::AbstractFloat=0.99)
    lower = quantile(skipmissing(vect), low_per)
    upper = quantile(skipmissing(vect), high_per)
    returnVect = vect[:]
    for (i, v) in enumerate(vect)
        if typeof(v) == Missing
            continue
        end
        if v > upper
            returnVect[i] = upper
        elseif v < lower
            returnVect[i] = lower
        end
    end
    return returnVect
end

df_coefs = @chain df_all_regs begin
    @rtransform(:r2 = r2(:reg))
    dropmissing
    @rtransform(:coefs = NamedTuple{(Symbol.(coefnames(:reg))..., ^(:var))}((coef(:reg)..., var(:reg))))
     _.coefs
    DataFrame
    rename("(Intercept)" => "int")
    @transform(
        :int = winsorize(:int),
        :mktrf = winsorize(:mktrf),
        :smb = winsorize(:smb),
        :hml = winsorize(:hml),
        :umd = winsorize(:umd),
        :var = winsorize(:var)
    )
end
##
coef_means = mean(df_coefs[:, 1:5] |> Matrix, dims=1)[1, :]
# [
#     0.0003416652538979692,
#     0.7557301603393168,
#     0.6035365140401592,
#     0.13841355976862427,
#     -0.0917780093857371,
# ]
coef_cov = cov(df_coefs[:, 1:5] |> Matrix)
# 2.52326e-6  -4.71181e-5   3.31182e-5  -7.05013e-5    1.32415e-5
# -4.71181e-5   0.27552      0.173033     0.0742338    -0.0448335
#  3.31182e-5   0.173033     0.389159     0.0460528    -0.0280572
# -7.05013e-5   0.0742338    0.0460528    0.40226      -0.000262376
#  1.32415e-5  -0.0448335   -0.0280572   -0.000262376   0.183439

fit_mle(Gamma, @rsubset(df_coefs, :var > 0).var)
# Gamma{Float64}(α=0.5792772418365709, θ=0.00448646291687775)