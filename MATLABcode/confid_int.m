


function [IRpoint, IRup, IRlo, IRmed, IRmean, Bevmean, Bevstd] = ...
    confid_int(Bcomp, VC_eps, cvec, dvec, IRhoriz, IRtype, sizesho, NMC, Tbig, CIperc, ispl, burnin);


%--------------------------------------------------------------------------
%
% DESCRIPTION:
%
% This routine computes the impulse response functions, their confidence 
% intervals and associated statistics
%
%--------------------------------------------------------------------------
%
% INPUTS:
%
% Bcomp:   matrix with structure of estimated reduced-form coefficients in 
%           companion form excluding all the deterministic terms. 
%           This includes the estimated parameters in the partition 
%           n x (n x (n x vlag)). The remaining partition is made up of a 
%           diagonalic identity matrix that pins down the lead-lag relation 
%           between variables. The size is (n x vlag) x (n x vlag).     
%
% VC_eps:  covariance matrix of reduced-form residuals, size n x n.
%
% cvec:    matrix with estimated constants; this returns the scalar 0 
%           if no constant is included in the model; 
%           otherwise, this is a matrix of size n x 1.
%
% dvec:    matrix with estimated parameters on the deterministic trends;
%           this returns the scalar 0 if no constant is included 
%           otherwise, this is a matrix of size n x 1.
%
% IRhoriz: number of periods for which point impulse responses are computed
%
%
% IRtype:  this string variable can be assigned the values'c' or 'g'
%
% sizesho: this is a singleton that is used to normalize the Choleski 
%           transformation of the variance-covariance matrix of 
%           reduced-form residuals. If it has a value equal to 0, 
%           alternative normalizations of the Choleski are used,
%           according to the following table
%
%
% Table 1: structural shock vector as a function of the two key inputs
% ---------------------------------------------------------------------------------------------------------
% |            |                                                                                          |
% |            |              IRtype=c                                    IRtype=g                        |
% |            |------------------------------------------------------------------------------------------|
% |            |                                                                                          |
% | sizesho=0  |         VC_eps_chol*eye(Nbig)                    VC_eps*diag(stdepsvec.^(-1))            |
% |            |                                                                                          |
% |sizesho~=0  |   VC_eps_chol*diag(sizesho./stdepsvec)   VC_eps_chol*diag(stdepsvec.^(-2))*diag(sizesho) |
% |            |                                                                                          |
% |            |                                                                                          |
% ---------------------------------------------------------------------------------------------------------
% Notes: 
% A. VC_eps denotes the variance-covariance decomposition of shocks of
% reduced-form model;
% B. Nbig denotes the number of variables;
% C. VC_eps_chol denotes the Choleski decomposition of the
% variance-covariance matrix of shocks of reduced-form model;
% D. stdepsvec denotes the standard deviations of shocks of reduced-form
% model.
%
%
% NMC:     number of Monte Carlo simulations for the density of simulated
%           impulse responses
%
% Tbig:    number of simulation periods after burn-in
%
% CIperc:  this is a 2 x 1 vector with lower and upper percentiles for 
%           the confidence itervals that are computed and, eventually,
%           plotted
%
% ispl:    1 to plot the impulse responses; 0 for no plotting 
%
% burnin:  number of observations used for the burn-in/initialization of 
%           the estimation of the VAR models
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% IRpoint: point estimate of impulse responses 
%
% IRup:    upper bound of impulse responses 
%
% IRlo:    lower bound of impulse responses
%
% IRmed:   median of impulse response distribution (at each horizon point)
%
% IRmean:  mean of impulse response distribution (at each horizon point)
%
% Bevmean: mean of companion matrix coefficients over the Monte Carlo runs
%
% Bevstd:  standard deviation of companion matrix coefficients over 
%           the Monte Carlo runs
%
%--------------------------------------------------------------------------
%
% Author:  Xiaohang Feng, 2019
%
%--------------------------------------------------------------------------



trendh = (1-burnin):Tbig;
Nbig   = length(VC_eps);
vlag   = length(Bcomp)/Nbig;

istr  = sum(abs(dvec))>0;
iscon = sum(abs(cvec))>0;

Bevsum   = zeros(Nbig*vlag+iscon+istr, Nbig);
Bevsqsum = zeros(Nbig*vlag+iscon+istr, Nbig);

allIRs = zeros(IRhoriz+1, Nbig, Nbig, NMC);
IRup   = zeros(IRhoriz+1, Nbig, Nbig);
IRlo   = IRup;
IRmed  = IRup;
IRmean = IRup;

IRpoint=irf(Bcomp, VC_eps, IRhoriz, IRtype, sizesho, 0);


if istr==0;
    dvec = zeros(Nbig, 1);
end

if iscon==0;
    cvec = zeros(Nbig, 1);
end

if vlag>1;
    cvec = [cvec; zeros((vlag-1)*Nbig, 1)];
    dvec = [dvec; zeros((vlag-1)*Nbig, 1)];
end

Idmatcomp      = eye(vlag*Nbig);
uncondmeancomp = inv(Idmatcomp-Bcomp)*cvec;
VC_epscompchol = zeros(Nbig*vlag, Nbig*vlag);
VC_epscompchol(1:Nbig, 1:Nbig) = chol(VC_eps)';



for mm=1:NMC;

   Ymat     = zeros(Tbig+burnin, Nbig);
   Ycomplag = uncondmeancomp;
   
   for tt=1:Tbig+burnin;
       Ycomp    = cvec+Bcomp*Ycomplag+dvec.*trendh(tt)+VC_epscompchol*randn(vlag*Nbig,1);
       Ycomplag = Ycomp;
       Ymat(tt, 1:Nbig) = Ycomp(1:Nbig)';
   end 
   
   Ymat = Ymat(burnin+1:end,:);
   [BcompMC, cvecMC, dvecMC, Bpl_evMC, VC_epsMC] = estim(Ymat, vlag, iscon, istr);
   Bevsum   = Bevsum   + Bpl_evMC;
   Bevsqsum = Bevsqsum + Bpl_evMC.^2;
   IRloc    = irf(BcompMC, VC_epsMC, IRhoriz, IRtype, sizesho, 0);   
   allIRs(:, :, :, mm) = IRloc; 
   
end


Bevmean = Bevsum/NMC;
Bevstd  = (Bevsqsum/NMC-Bevmean.^2).^0.5;


for hh=1:IRhoriz+1;
    
    for ii=1:Nbig;
        
        for jj=1:Nbig;
            sliceloc   = squeeze(allIRs(hh,ii,jj,:));
            IRpercsloc = quantile(sliceloc,[CIperc(1);0.5;CIperc(2)]);
            IRmeanloc  = mean(sliceloc);
            IRlo(hh,ii,jj)   = IRpercsloc(1);
            IRup(hh,ii,jj)   = IRpercsloc(2);
            IRmed(hh,ii,jj)  = IRpercsloc(3);
            IRmean(hh,ii,jj) = IRmeanloc;
        end
        
    end
    
end





if ispl==1;
    
    xaxis=0:IRhoriz;
    zerol=zeros(IRhoriz+1,1);
    plotcount=1;
    
    for ii=1:Nbig;
        
        yaxmin=min(min(IRlo(:,ii,:)));
        yaxmin=yaxmin-0.1*abs(yaxmin);
        yaxmax=max(max(IRup(:,ii,:)));
        yaxmax=yaxmax+0.1*abs(yaxmax);
        
        for jj=1:Nbig;
            subplot(Nbig,Nbig,plotcount);  plot(xaxis, IRpoint(:,ii,jj), xaxis, IRlo(:,ii,jj), xaxis, IRup(:,ii,jj), xaxis, IRmed(:,ii,jj), xaxis, zerol);  axis([0 (IRhoriz) yaxmin yaxmax]); title(['shock' int2str(jj) ' to var ' int2str(ii)]);
            plotcount=plotcount+1;
        end
        
    end

end



