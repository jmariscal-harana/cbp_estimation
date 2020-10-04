%%	Main code for the 2-element and 3-element Windkessel datasets
%	Input
%	-REF: parameters used to generate virtual patients
%	-PATHS: folders where required data is stored
%
%	Output
%	-REF: virtual patient data
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (30/07/19) - based on Wk_CBP_study
%
%==========================================================================
clear;
close all;
dbstop if error


%%	REF input parameters
Dataset				= 'Wk2';				% Reference dataset (virtual): 'Wk2', 'Wk3'
Inflow_baseline		= 'Inflow_old';			% Use Inflow_old as it has a marked inflection point
Plots_REF			= 1;					% Plot virtual inflow waveforms
Windkessel_REF		= Dataset;				% Only used when generating Wk2 or Wk3 datasets
								
									
%%	PATHS
PATHS.Root							= '~/Haemodynamic_Tools/Version6/';	%Root folder containing all other folders
PATHS.Wk_study						= [PATHS.Root,'Windkessel/'];	% Path to the Windkessel study folder
PATHS.Figures						= [PATHS.Wk_study,'Figures/'];	% Path to the figures folder
PATHS.Matlab						= [PATHS.Wk_study,'Matlab/'];	% Path to the main Matlab folder
PATHS.Data_inflow					= [PATHS.Root,'GenerateInputFiles/ParametricInflow/ORourke_inflows/'];	% Path to the inflow data
	
PATHS.Input							= [PATHS.Wk_study,Dataset,'_input/'];

addpath(PATHS.Matlab)
addpath([PATHS.Root,'Others/P_0_iteration'])


%%	REF data: Virtual or Clinical reference dataset
%	Extract CV parameters from the clinical literature
REF = Wk_VirtualParameters(PATHS);

%	Generate Wk2/Wk3 virtual dataset
REF	= Wk_ReferencePopulation(PATHS,REF,Inflow_baseline,Plots_REF,Windkessel_REF);
save([PATHS.Input,Dataset,'_reference_v2.mat'],'REF')

				
%%	
restoredefaultpath