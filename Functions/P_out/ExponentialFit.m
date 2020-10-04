%% Exponential fit  - calculate P_out, tau from pressure waveform
%   Input
%	-P: pressure waveform
%	-T: cardiac period
%	-LVET: time of start of diastole (if missing it is found by fsg721; if ==0 it is done interactively
%	-Param_num: reference parameter value (no need to estimate)
%	-Param_str: 'P_out' or 'tau'
%	-Plot_on: 1 = plot on
%
%	Output
%	-tau: decay constant
%	-P_out: outflow pressure (theoretical asymptotic value)
%	-P_d_1: pressure value at the 
%	-LVET: detected time of the dichrotic notch
%	-Pn: pressure at dichrotic notch
%
%	uses function fsg721 and assumes data are smooth enough to find derivative using
%	7-pt SG filter
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (28/05/18) - Based on kreservoir_v10.m (Copyright 2008 Kim H Parker)
%	v1.1 (01/07/19) - clean up
%
%==========================================================================

function [tau,P_out,P_d_1,LVET,Fit]= ExponentialFit(P,T,LVET,Param_num,Param_str,Plot_on)

%	P must be a row vector
P = P(:)';
t=linspace(0,T,length(P));

DBP = min(P);

%	Read or determine LVET
if (nargin >= 3)
	if (LVET~=0)	%LVET given as input
		[~,nn]=min(abs(t-LVET));
	else	%determine LVET interactively
		figure;
		plot(P);
		disp('click on start of diastole');
		inp=ginput(1);
		nn=round(inp(1,1));
	end
else	% determine LVET automatically from the minimum dP
	warning('Calculating LVET from min(dP); can you provide LVET instead?')
	dp=fsg721(P);
	[~,nn]=min(dp);
end

LVET	= t(nn);
P_d_1	= P(nn);

%	Calculate pressure during diastole using: (P_sys - P_out)*exp(-t/tau) + P_out
P_d		= P(nn:end);
t_d		= t(nn:end)-LVET;

switch Param_str
	case 'P_out'
		P_out = Param_num;
		tau = 1;
		fun_exp = @(x) (x(1) - P_out) * exp(-t_d./x(2)) + P_out;
		fun_error = @(x) sum(abs(P_d - ((x(1) - P_out) * exp(-t_d./x(2)) + P_out)));
		x0 = [P_d_1,tau];
		
	case 'tau'
		P_out = DBP/2;
		tau = Param_num;
		fun_exp = @(x) (x(1) - x(2)) * exp(-t_d./tau) + x(2);
		fun_error = @(x) sum(abs(P_d - ((x(1) - x(2)) * exp(-t_d./tau) + x(2))));
		x0 = [P_d_1,P_out];
		
	case 'none'
		P_out = DBP/2;
		tau = 1;
		fun_exp = @(x) (x(1) - x(2)) * exp(-t_d./x(3)) + x(2);
		fun_error = @(x) sum(abs(P_d - ((x(1) - x(2)) * exp(-t_d./x(3)) + x(2))));
		x0 = [P_d_1,P_out,tau];
		
end

options = optimset('MaxFunEvals',1000);
x = fminsearch(fun_error,x0,options);

if (nargin == 6)
	switch Plot_on
		case 1
			switch Param_str
				case 'P_out'
					disp('Plotting diastolic fit for a pressure waveform, P_out fixed')
				case 'tau'
					disp('Plotting diastolic fit for a pressure waveform, tau fixed')
			end
			figure
			hold on
			plot(t,P,'k','LineWidth',2)
			plot(t_d+LVET,fun_exp(x0),'--b')
			plot(t_d+LVET,fun_exp(x),'-.r')
			hold off
			
			title('Exponential fit')
			xlabel('Time [s]')
			ylabel('Pressure [mmHg]')
			legend('Pressure','initial fit','optimisation')
			set(gca,'FontSize',24)
			
% 			addpath('~/Haemodynamic_Tools/Version6/Windkessel/Matlab/')
% 			SaveEstimationPlots('~/Haemodynamic_Tools/Version6/Windkessel/Parameter_study/',...
% 				['Exp_fit_fixed_',Param_str])
	end
end

% return variables
P_d_1 = x(1);
switch Param_str
	case 'P_out'
		tau		= x(2);
	case 'tau'
		P_out	= x(2);
	case 'none'
		P_out	= x(2);
		tau		= x(3);
end

if tau <= 0
	warning('tau < 0, P_out set to 0, recalculating tau')
	P_out = 0;
	tau = ExponentialFit(P,T,LVET,P_out,'P_out',Plot_on);
end

% global P_out_negative_counter

if P_out >= DBP
	warning('P_out >= DBP; P_out set to DBP/2, recalculating tau')
	P_out = DBP/2;
	tau = ExponentialFit(P,T,LVET,P_out,'P_out',Plot_on);
	
% 	P_out_negative_counter = P_out_negative_counter+1;
	
elseif P_out < 0
	warning('P_out < 0; P_out set to 0, recalculating tau')
	P_out = 0;
	tau = ExponentialFit(P,T,LVET,P_out,'P_out',Plot_on);
	
% 	P_out_negative_counter = P_out_negative_counter+1;
	
end

Fit = fun_exp(x);

end


%% Savitsky-Golay first derivative filter function
% dx=fsg721(x)
% 	7 point SavGol filter, 2nd order polynomial, 1st derivative
%	input 	x
%	output	dx
%	corrected for time shift

function dx=fsg721(x)

% 2nd order polynomial
C=[0.107143,0.071429,0.035714];

B=zeros(1,7);
for i=1:3 
	B(i)=C(i); 
end	
B(4)=0.0;
for i=5:7
	B(i)=-C(8-i);
end
A=[1,0];
	
s=size(x,2);
dx=filter(B,A,x);
dx=[dx(7),dx(7),dx(7),dx(7:s),dx(s),dx(s),dx(s)];

end
