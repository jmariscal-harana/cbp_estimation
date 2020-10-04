%%	Calculates the reference 2-element or 3-element Windkessel pressure
%	Input
%	-REF: parameters used to generate virtual patients
%	-PATHS: folders where required data is stored
%
%	Output
%	-REF_out: virtual patient data
%
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (17/10/17) 
%==========================================================================

function [REF_out] = Wk_ReferencePressure(PATHS,REF,Windkessel)
%%	INPUT
Inflow_data	= load([PATHS.Wk_study,'Inflow/',REF.INFO.Inflow,'_t.txt']);	% Parametric inflow waveform

REF.t_in = Inflow_data(:,1);		% s
REF.Q_in = Inflow_data(:,2);		% mL/s

switch Windkessel
	case 'Wk2'
		REF.tau		= REF.R_T*REF.C_T;
	case 'Wk3'
		REF.tau		= (REF.R_T - REF.Z_0)*REF.C_T;
end

%%	0-D Windkessel calculation
%	Initial quasi-steady state (QSS) estimation
REF	= P_0_iteration(REF,0,Windkessel);


%%	OUTPUT
REF_out.tau			= REF.tau;
REF_out.t_in		= REF.t_in;
REF_out.Q_in		= REF.Q_in;
REF_out.Pressure	= REF.Pressure;

%	Not saved at the moment to save HDD memory:
% REF_out.t_in_cycles	= t_in_temp;

%	Not saved at the moment to save HDD memory:
% REF_out.P_Wk2_exponential = P_Wk2_exponential;
% REF_out.P_Wk2_integral = P_Wk2_integral;
% REF_out.P_Wk2_exp_int = P_Wk2_exp_int;


end

