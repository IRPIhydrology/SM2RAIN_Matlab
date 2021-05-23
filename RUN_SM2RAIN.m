% ... few lines to run the different version of SM2RAIN code
% Application with Cerbara dataset, 24-hour rainfall estimation from hourly observations
%% -----------------------
% SM2RAIN.m: first version of SM2RAIN (3 parameters Z, a, b)
% -----------------------
name='CER_1hour_2011';                    % name of the input file
AGGR=24;                                  % aggregation period: 24 hour=1 day
data=load(['Test_data\',name,'.txt']);    % load of input data
X=cal_SM2RAIN(data,AGGR);                 % SM2RAIN calibration
FIG=1;                                    % 0 no figure, 1 create the figure
namefig=name;                             % name of the figure
SM2RAIN(data,X,AGGR,FIG,namefig);         % RUN SM2RAIN

%% -----------------------
% SM2RAIN_T.m: SM2RAIN with the exponential filter (4 parameters Z, a, b, T)
% -----------------------
name='CER_1hour_2011';                        % name of the input file
AGGR=24;                                      % aggregation period: 24 hour=1 day
X_ini=[0.1,0.05,0.1,0.1]';                    % initial conditions for the parameter values (dimensionless)
data=load(['Test_data\',name,'.txt']);        % load of input data
FIG=0;                                        % 0 no figure, 1 create the figure
namefig=name;                                 % name of the figure
X=cal_SM2RAIN_T(data,AGGR,X_ini,FIG,namefig); % SM2RAIN calibration
FIG=1;                                        % 0 no figure, 1 create the figure
SM2RAIN_T(data,X,AGGR,FIG,namefig);           % RUN SM2RAIN

%% -----------------------
% SM2RAIN_Tpot2.m: SM2RAIN with the advanced exponential filter (5 parameters Z, a, b, Tpot, Tbase)
% -----------------------
name='CER_1hour_2011';                            % name of the input file
AGGR=24;                                          % aggregation period: 24 hour=1 day
X_ini=[0.1,0.05,0.1,0.1,.1]';                     % initial conditions for the parameter values (dimensionless)
data=load(['Test_data\',name,'.txt']);            % load of input data
FIG=0;                                            % 0 no figure, 1 create the figure
namefig=name;                                     % name of the figure
X=cal_SM2RAIN_Tpot2(data,AGGR,X_ini,FIG,namefig); % SM2RAIN calibration
FIG=1;                                            % 0 no figure, 1 create the figure
SM2RAIN_Tpot2(data,X,AGGR,FIG,namefig);           % RUN SM2RAIN
