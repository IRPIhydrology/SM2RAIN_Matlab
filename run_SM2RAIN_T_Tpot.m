clc,clear,close all
currdir=pwd;
a=dir('*.txt');
ele={a(:,1).name};
NN=24;FIG=1;
X_ini=[0.1,0.05,0.1,0.1,.1]';
% X_ini=[0.1,0.05,0.1,0.1]';
for ii=1%:length(ele)
    data=load(ele{ii});
    namefig=ele{ii};
    namefig=namefig(1:end-4);
    [~,R(ii),KGE(ii)]=cal_SM2RAIN_Tpot2(data,NN,X_ini,FIG,[namefig,'_Tpot']);
%     [X,R(ii)]=cal_SM2RAIN_T(data,NN,X_ini,FIG,[namefig]);
end

%%
% NN=24
% name='ING_1hour_0912_T.txt';
% data=load(name);
% data=sortrows(data);
% namefig=name;
% namefig=namefig(1:end-4);
% [~,R(ii),KGE(ii)]=cal_SM2RAIN_T(data,NN,X_ini,FIG,namefig);