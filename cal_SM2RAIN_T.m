function [X,R,RMSE,Psim,SM]=cal_SM2RAIN_T(name,NN,X_ini,FIG,namefig)
if nargin==2,X_ini=[0.1,0.05,0.1,0.1]';end
if nargin<4,FIG=0;end
[RES,FVAL,EXITFLAG,OUTPUT]=fmincon(@calibOK,X_ini,[],[],[],[],...
     zeros(4,1),ones(4,1),[],...
     optimset('Display','off','MaxIter',300,'MaxFunEvals',500,...
     'TolFun',1E-20,'TolCon',4,'Largescale','off','Algorithm','active-set'),...
     name,NN);
% [RES,FVAL,EXITFLAG,OUTPUT]=fmincon(@calibOK,X_ini,[],[],[],[],...
%      zeros(4,1),ones(4,1),[],...
%      optimset('Display','off','MaxIter',300,'MaxFunEvals',500,...
%      'TolFun',1E-20,'TolCon',4,'Largescale','off','Algorithm','interior-point'),...
%      name,NN);
X=convert_adim(RES);
[R,RMSE,NS,KGE,Psim,SM]=SM2RAIN_T(name,X,NN,FIG,namefig);

%---------------------------------------------------------------------------------
function [err]=calibOK(X_0,name,NN)

X=convert_adim(X_0);
[R,RMSE,NS,KGE]=SM2RAIN_T(name,X,NN,0);
% err=RMSEcum+2*(1-R^2);
err=1-NS;
% save X_PAR

%---------------------------------------------------------------------------------
function X=convert_adim(X_0)
% LOW=[   1,   0.0,  1,  0.0, 0.000]';
% UP =[ 500, 800.0, 50, 20.0, 0.040]';
LOW=[ 20,   0.01,   1,  0.0]';
UP =[ 20, 160.0,  50,  0.0]';
X=LOW+(UP-LOW).*X_0;
