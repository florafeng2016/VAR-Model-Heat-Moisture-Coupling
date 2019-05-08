% code ceated by Xiaohang Feng

datatrain = xlsread('FileS1.xlsx');
datatest = xlsread('FileS2.xlsx');
% read data

sm10 = datatrain(3:3895,3);
st10 = datatrain(3:3895,4);
sm30 = datatrain(3:3895,5);
st30 = datatrain(3:3895,6);
sm90 = datatrain(3:3895,7);
st90 = datatrain(3:3895,8);
D1 = LagOp({1 -1},'Lags',[0,1]);D12 = LagOp({1 -1},'Lags',[0,24]);D = D1*D12;
dsm10 = filter(D,sm10);dst10 = filter(D,st10);dsm30 = filter(D,sm30);
dst30 = filter(D,st30);dsm90 = filter(D,sm90);dst90 = filter(D,st90);
% seasonal difference of the original data
% (note: to facilitate the readers who want to testify our results, the data we
% procided are data after seasonal difference)

[h,pValue,stat,cValue,reg] = adftest(dsm10,'model','TS','lags',0:1);
[h,pValue,stat,cValue,reg] = adftest(dst10,'model','TS','lags',0:1);
[h,pValue,stat,cValue,reg] = adftest(dsm30,'model','TS','lags',0:1);
[h,pValue,stat,cValue,reg] = adftest(dst30,'model','TS','lags',0:1);
[h,pValue,stat,cValue,reg] = adftest(dsm90,'model','TS','lags',0:1);
[h,pValue,stat,cValue,reg] = adftest(dst90,'model','TS','lags',0:1);
% conduct the ADF tests on the original data to ensure that all time series
% are stationary
% (The method can vary accoridng to the pattern of the time series)

Y1 = [dsm10,dst10];Y2 = [dsm30,dst30];Y3 = [dsm90,dst90];
% "sm10" denotes soil moisture at the depth of 10 cm, and "st10" denotes soil temperature at teh depth of 10 cm.

for ii = 1:1:12
    Mdl = varm(2,ii);EstMdl = estimate(Mdl,Y1);summarize(EstMdl)
end
for ii = 1:1:12
    Mdl = varm(2,ii);
    EstMdl = estimate(Mdl,Y2);
    summarize(EstMdl)
end
for ii = 1:1:12
    Mdl = varm(2,ii);
    EstMdl = estimate(Mdl,Y3);
    summarize(EstMdl)
end
% obtain the AIC and BIC for different order in the summary, and select the optimal order
% with the lowest AIC and BIC value
% build vector autoregression models
% (note that when building the model the lag order should be set according to the comparison results of AIC and BIC)

sizesho = 0;IRhoriz=120;IRtype='g';
[Bcomp, cvec, dvec, Bpl_ev, VC_eps] = estim(Y1,5,1,0);IR = irf(Bcomp,VC_eps,120,'g',0,0);
[IRpoint,IRup,IRlo,IRmed,IRmean,Bevmean,Bevstd] = ...
    confid_int(Bcomp,VC_eps,cvec,dvec,IRhoriz,IRtype,sizesho,NMC,Tbig,CIperc,ispl,burnin);

[Bcomp, cvec, dvec, Bpl_ev, VC_eps] = estim(Y2,8,1,0);IR = irf(Bcomp,VC_eps,120,'g',0,0);
[IRpoint,IRup,IRlo,IRmed,IRmean,Bevmean,Bevstd] = ...
    confid_int(Bcomp,VC_eps,cvec,dvec,IRhoriz,IRtype,sizesho,NMC,Tbig,CIperc,ispl,burnin);

[Bcomp, cvec, dvec, Bpl_ev, VC_eps] = estim(Y3,7,1,0);IR = irf(Bcomp,VC_eps,120,'g',0,0);
[IRpoint,IRup,IRlo,IRmed,IRmean,Bevmean,Bevstd] = ...
    confid_int(Bcomp,VC_eps,cvec,dvec,IRhoriz,IRtype,sizesho,NMC,Tbig,CIperc,ispl,burnin);
% build the impulse response functions
% (note that the lag order after "Yi" should be set according to the
% optimal lag order obtained for the AIC and BIC tests)

[F,c_v] = granger_cause(dst10,dsm10,0.05,8)
[F,c_v] = granger_cause(dsm10,dst10,0.05,8)
[F,c_v] = granger_cause(dst30,dsm30,0.05,8)
[F,c_v] = granger_cause(dsm30,dst30,0.05,8)
[F,c_v] = granger_cause(dst90,dsm90,0.05,8)
[F,c_v] = granger_cause(dsm90,dst90,0.05,8)
% Granger causality test
% (note that the lag order should be set according to the comparison results of AIC and BIC)

% the process for the test data is the same as that fpr training data
