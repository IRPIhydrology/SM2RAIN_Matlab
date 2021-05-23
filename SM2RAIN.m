function [NS,R,RMSE,PPsim,Psim,Pobs,D]=SM2RAIN(data,PAR,NN,FIG,namefig)

D=data(1:end-1,1); Pobs=data(1:end-1,3); SM=data(:,2);
Z=PAR(1);
a=PAR(2);
b=PAR(3);
deltaSM=0.0005;

N=length(Pobs);
Psim=zeros(N,1);
for t=2:N
    if (SM(t)-SM(t-1))>deltaSM
%         SMav=(SM(t)+SM(t-1))./2;
%         Psim(t-1)=Z.*(SM(t)-SM(t-1))+a.*SMav.^b;
        Psim(t-1)=Z.*(SM(t)-SM(t-1))+(a.*SM(t).^b+a.*SM(t-1).^b)./2;
    end
end

Psim(isnan(SM(2:end)))=NaN;
Pobs(isnan(Psim))=NaN;
PPsim=Psim;
% Temporal aggregation
if NN>1
    MM=length(Psim);
    L=floor(MM/NN);
    Psim=reshape(Psim(1:L*NN),NN,L);
    Psim=sum(Psim)';
    Pobs=reshape(Pobs(1:L*NN),NN,L);
    Pobs=sum(Pobs)';
    SM=reshape(SM(1:L*NN),NN,L);
    SM=mean(SM)';
    D=reshape(D(1:L*NN),NN,L);
    D=mean(D)';
end

% Calculation of model performance
fRMSE=(nanmean((Psim-Pobs).^2).^0.5)./nanstd(Pobs);
RMSE=(nanmean((Psim-Pobs).^2).^0.5);
NS=1-nansum((Psim-Pobs).^2)./nansum((Pobs-nanmean(Pobs)).^2);
R=corr(Psim,Pobs,'rows','complete');
BIAS=nanmean(Psim-Pobs);
KGE=klinggupta(Psim,Pobs);
ANSE=1-nansum((Pobs+nanmean(Pobs)).*(Psim-Pobs).^2)./...
    nansum((Pobs+nanmean(Pobs)).*(Pobs-nanmean(Pobs)).^2);
X=0.5;
THR=Pobs>X; PIO=Psim>X;
H_SIM=sum(THR==1 & PIO==1); M_SIM=sum(THR==1 & PIO==0);
F_SIM=sum(THR==0 & PIO==1); C_SIM=sum(THR==0 & PIO==0);
FAR=(F_SIM)./(H_SIM+F_SIM);
TS=(H_SIM)./(H_SIM+M_SIM+F_SIM);
POD=(H_SIM)./(H_SIM+M_SIM);

% Figure
if FIG==1
    clf
    set(gcf,'paperpositionmode','manual','paperposition',[1 1 20 10],'Color','white')
    set(gcf,'position',[100   100   1200   400])
    
    axes('Position',[0.1 0.3 0.8 0.40]);
    set(gca,'Fontsize',13)
    s=([' R= ',num2str(R,'%4.3f'),...
        ' NS= ',num2str(NS,'%4.3f'),...
        ' KGE= ',num2str(KGE,'%4.3f'),...
        ' BIAS= ',num2str(BIAS,'%4.3f'),...
        ' RMSE= ',num2str(RMSE,'%4.3f'),...
        ' POD= ',num2str(POD,'%4.3f'),...
        ' FAR= ',num2str(FAR,'%4.3f'),...
        ' TS= ',num2str(TS,'%4.3f')]);
    nname=namefig; nname(nname=='_')='-';
    title(['\bf',nname,' ',s]);
    hold on
    plot(D,Pobs,'g','Linewidth',3)
    plot(D,Psim,'r--o','Linewidth',2,'Markersize',2);
    grid on, box on
    ylabel('rain [mm]')
    legend('P_o_b_s','P_s_i_m','orientation','horizontal');
    datetick('x',12)
    axis([data(1,1)-1 data(end,1)+1 0 nanmax(Pobs)+0.05])
    set(gca,'Xticklabel','')
    
    axes('Position',[0.1 0.1 0.8 0.20]);
    set(gca,'Fontsize',13)
    hold on
    plot(data(:,1),data(:,2),'color',[.5 .5 .5],'Linewidth',3);
    grid on, box on
    ylabel('saturation [-]')
    datetick('x',12)
    axis([data(1,1)-1 data(end,1)+1 -0.05 1.1])
    s=num2str(NN);
    if NN<10,s=['0',num2str(NN)];end
    
    export_fig([namefig,'_SM2R_',s], '-png','-q60','-r150')
    %     save res
end
