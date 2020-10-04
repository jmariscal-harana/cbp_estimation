%%	Generate virtual inflow waveform based on SV and HR
%	Input
%	-PATHS: folders where required data is stored
%	-INFLOW: parameters used to generate virtual patients
%
%	Output
%	-INFLOW: virtual flow waveform
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (30/07/19) - based on Wk_CBP_study
%
%==========================================================================
clear

%%	Paths
PATHS.Wk_study		= [cd,'/'];
PATHS.Data_inflow	= 'ORourke_inflows/';	% Path to the inflow data
PATHS.Matlab		= 'Matlab/';	% Matlab path
addpath(PATHS.Matlab)

Plots_REF			= 1;	%1: show plots


%%	Input parameters
INFLOW.Inflow	= 'Test_1';	%Inflow name or code
INFLOW.SV		= 100;	%Stroke volume [mL]
INFLOW.HR		= 60;	%Heart rate [bpm]
INFLOW.dt		= 10^-3;	%Waveform time step


%%	Generate inflow
INFLOW			= GenerateInflow(PATHS,INFLOW,Plots_REF);

%Flow and time vectors
Q				= INFLOW.Q;
t 				= INFLOW.t;

%Updated parameter values after interpolation
HR				= INFLOW.HR;
SV				= INFLOW.SV;

%Calculated during the interpolation
LVET			= INFLOW.LVET;


%%
restoredefaultpath