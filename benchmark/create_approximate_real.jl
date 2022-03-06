using Distributions, DataFrames, DataFramesMeta, Arrow, Statistics, FixedEffectModels, BusinessDays, Dates, CSV
using Revise
using AbnormalReturns

##

df_ff_data = @chain Arrow.Table(joinpath("temp_data", "ff_data.feather")) begin
    DataFrame
    copy
    @rtransform(:mkt = :mktrf + :rf)
end

ff_data_max = maximum(df_ff_data.date)

@time df_crsp_raw = @chain Arrow.Table(joinpath("..", "..", "Share Pledging Analysis", "data", "general csv files", "crsp_daily.feather")) begin
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
coef_cov = cov(df_coefs[:, 1:5] |> Matrix)
fit_mle(Gamma, @rsubset(df_coefs, :var > 0).var)