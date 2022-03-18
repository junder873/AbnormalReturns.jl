data _null_; 
    rc=dlgcdir("Documents/GitHub/AbnormalReturns.jl");
    put rc=;
run;
proc import file="benchmark/data/event_dates.csv"
	out = event_dates
	dbms = csv
	replace;
run;

proc import file="benchmark/data/firm_ret.csv"
	out = firm_data
	dbms = csv
	replace;
run;

proc import file="benchmark/data/mkt_ret.csv"
	out = mkt_data
	dbms = csv
	replace;
run;

proc sql; create table unique_evts as 
	select distinct firm_id, est_window_start, est_window_end, event_window_start, event_window_end from event_dates;
quit;

proc sql;
	create view evt_rets as
	select * from
	unique_evts a
	left join (
		select b.firm_id, b.ret, c.* from firm_data b
		inner join mkt_data c
		on b.date = c.date
	) d
	on a.firm_id = d.firm_id
	and a.est_window_start <= d.date <= a.est_window_end
	order by a.firm_id, a.est_window_start, a.est_window_end, d.date;
quit;

proc reg data=evt_rets edf outest=params  noprint;
   by firm_id est_window_start est_window_end event_window_start event_window_end;
   eq1: model ret=mkt;
   eq2: model ret=mkt smb hml umd;
run;

data abrets; merge
  evt_rets(where=(event_window_start<=date<=event_window_end) in=a)
  params (where=(_model_='eq1')
     keep=firm_id event_window_start event_window_end _model_ _rmse_ intercept mkt
     rename=(_rmse_=std1 intercept=alpha1 mkt=beta1))
  params (where=(_model_='eq2')
     keep=firm_id event_window_start event_window_end _model_ _rmse_ intercept mkt smb hml umd
     rename=(_rmse_=std2 intercept=alpha2 mkt=beta2 smb=sminb2 hml=hminl2 umd=umind2));
  by firm_id event_window_start event_window_end; 
  ar1=std1**2;var2=std2**2;
  expret1=alpha1+beta1*mkt; abret1=ret-expret1;
  expret2=alpha2+beta2*mkt+sminb2*smb+hminl2*hml+umind2*umd; abret2=ret-expret2;
run;

proc expand data=abrets out=car method=none;
    by firm_id event_window_start event_window_end;
    convert expret1 =bhar1 /transformout=(+1 cupord - 1);
    convert expret2 =bhar2 /transformout=(+1 cupord - 1);
run;