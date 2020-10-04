function [INFLOW] = GenerateInflow(PATHS,INFLOW,Plots_REF)
%%	Reference cardiac parameters
INFLOW.T			= 60/INFLOW.HR;	% Cardiac period [s]
INFLOW.CO           = INFLOW.SV*INFLOW.HR/60;  % Cardiac output [mL/s]
INFLOW.mult_period  = 1.0;          % Multiply systole period duration
INFLOW.mult_sysPeak = 1.0;          % Multiply distribution of times of systole
INFLOW.mult_sys     = 1.0;          % Multiply distribution of times of systole
INFLOW.age          = 'old';        % Select type ('young', 'old', 'avg') of flow waves based on O-Rourke's paper
INFLOW.reverseflow  = 1;            % 1: Reverse flow activated; 0: Without reverse flow
INFLOW.Plots		= Plots_REF;

%% Remarkable points in the inflow curve
INFLOW.Qmax = 1;    %Systolic peak flow
INFLOW.Qi   = 0.32; %Flow at inflection point
INFLOW.Qi2  = 0.4;  %Flow at inflection point in systolic decrease (old)
INFLOW.Qmin = -0.1; %Minimum flow during dicrotic notch

INFLOW.ti_0   = 0.020*INFLOW.T; %Time of inflection point in systolic increase
INFLOW.tmax_0 = 0.085*INFLOW.T; %Time of systolic peak
INFLOW.ti2_0  = 0.220*INFLOW.T; %Time of inflection point in systolic decrease (old only)
INFLOW.ts_0   = 0.290*INFLOW.T; %Time of start of dicrotic notch
INFLOW.tmin_0 = 0.300*INFLOW.T; %Time of minimum flow during dicrotic notch
INFLOW.td_0   = 0.330*INFLOW.T; %Time of start of diastole


%%	Ejection time calculation
ET_from_HR		= .266 - .0021*(INFLOW.HR - 73);
ET_from_SV		= .266 + .0017*(INFLOW.SV - 82);
ET_from_HR_SV	= .266 + .0011*(INFLOW.SV - 82) - .0009*(INFLOW.HR - 73);
ET_factor		= ET_from_HR_SV/INFLOW.td_0;

INFLOW.ti_0		= ET_factor*INFLOW.ti_0;
INFLOW.tmax_0	= ET_factor*INFLOW.tmax_0;
INFLOW.ti2_0	= ET_factor*INFLOW.ti2_0;
INFLOW.ts_0		= ET_factor*INFLOW.ts_0;
INFLOW.tmin_0	= ET_factor*INFLOW.tmin_0;
INFLOW.td_0		= ET_factor*INFLOW.td_0;


%% INFLOW - Reflective flow as sum of harmonics in bcs file
INFLOW = Create_parametric_inflow(PATHS,INFLOW);


end
