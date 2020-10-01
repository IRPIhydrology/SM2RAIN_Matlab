function [X,R,RMSE,Psim,SM]=cal_SM2RAIN_Tpot2_E(name,NN,X_ini,FIG,namefig)
if nargin==2,X_ini=[0.1,0.05,0.1,0.1,0.1,.5]';end
if nargin<4,FIG=0;end
[RES,FVAL,EXITFLAG,OUTPUT]=fmincon(@calibOK,X_ini,[],[],[],[],...
     zeros(6,1),ones(6,1),[],...
     optimset('Display','off','MaxIter',300,'MaxFunEvals',500,...
     'TolFun',1E-20,'Largescale','off','Algorithm','active-set'),...
     name,NN);
X=convert_adim(RES);
[R,RMSE,NS,KGE,Psim,SM]=SM2RAIN_Tpot2_E(name,X,NN,FIG,namefig);

%---------------------------------------------------------------------------------
function [err]=calibOK(X_0,name,NN)

X=convert_adim(X_0);
[R,RMSE,NS,KGE]=SM2RAIN_Tpot2_E(name,X,NN,0);
% err=RMSEcum+2*(1-R^2);
err=1-KGE;
% save X_PAR

%---------------------------------------------------------------------------------
function X=convert_adim(X_0)
% LOW=[   1,   0.0,  1,  0.0, 0.000]';
% UP =[ 500, 800.0, 50, 20.0, 0.040]';
LOW=[  10,   0.0,   1,  0.05, 0.05, 0.02]';
UP =[ 400, 160.0,  50,  0.75, 3.00, 1.98]';
X=LOW+(UP-LOW).*X_0;
