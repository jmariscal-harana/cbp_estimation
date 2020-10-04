%%	cBP estimation error analysis
%	Input
%	-REF: reference data (in silico or in vivo)
%	-EST: estimated data using Aortic1D model
%	-PATHS: folders for data storage
%	-PLOTS: choose the type of plots which are generated
%
%	Output
%	-Saves plots to PATHS.Figures
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (15/11/18) - Based on Wk_DataAnalysis from Windkessel study
%
%==========================================================================
format compact;
close all;
clear;
% clc;
% dbstop if error
% dbclear all


%%	Scaling factors for pressure, flow and area
Scale_P = 133.32;
Scale_U = 1;
Scale_A = 1;
Scale_Q = 1e6;


%%	PATHS
PATHS.Root				= '~/Haemodynamic_Tools/Version6/';


%%	Load REF dataset
REF_Dataset				= 'pwdb_4374';	%'Isra_AoCoarctation', 'Sam_Normotensive', 'Sam_Hypertensive', 'pwdb_6', 'pwdb_78', 'pwdb_4374'
Scenario				= 'Waveform';	% 'Waveform': pressure waveform available; 'BP' blood pressure values
EST_method				= '3419_0_17';
% Scenario				= 'BP';	% 'Waveform': pressure waveform available; 'BP' blood pressure values
% EST_method				= '5226_0_3';

switch Scenario
case 'Waveform'
	Sc = 'Sc1';
case 'BP'
	Sc = 'Sc2';
end

%%%
disp('Reference dataset: ')
disp(REF_Dataset)

warning('Pressure should be in mmHg')

switch REF_Dataset
	case 'Isra_AoCoarctation'	%Clinical: aortic coarctation
		PATHS.REF_data	= '~/Clinical_Data/Isra_AoCoarctation/';
		PATHS.Figures	= [PATHS.REF_data,'Figures/'];
		load([PATHS.REF_data,REF_Dataset,'_reference.mat'],'REF')
	
		%	Store waveform pressure data into REF_temp to save space
		for jj = 1:length(REF)
			REF_temp(jj).ID				= REF(jj).INFO.ID;
			REF_temp(jj).Pressure		= REF(jj).Desc_Ao_I.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).t				= linspace(0,REF(jj).HAEMO.T,length(REF_temp(jj).Pressure_CBP))';
			t							= linspace(0,REF(jj).HAEMO.T,length(REF_temp(jj).Pressure))';
			
% 			figure, hold on
% 			plot(t,REF_temp(jj).Pressure,'--k')
% 			plot(REF_temp(jj).t,REF_temp(jj).Pressure_CBP,'k')
		end
		
	case {'Sam_Normotensive','Sam_Hypertensive'}	%Clinical: normotensive or hypertensive
		PATHS.REF_data	= ['~/Clinical_data/',REF_Dataset,'/'];
		PATHS.Figures	= [PATHS.REF_data,'Figures/'];
		load([PATHS.REF_data,REF_Dataset(5:end),'_jmh_v2.mat'],'REF')

		%	Store waveform pressure data into REF_temp to save space
		for jj = 1:length(REF)
			REF_temp(jj).ID				= REF(jj).INFO.ID;
			REF_temp(jj).Pressure		= REF(jj).Carotid.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;			
			REF_temp(jj).t				= linspace(0,REF(jj).HAEMO.T,length(REF_temp(jj).Pressure))';
		end
		
	case {'pwdb_6','pwdb_78','pwdb_4374'}	%In silico: 116-artery (78 or all 4374 subjects)
		PATHS.REF_data	= [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Simulations/'];
		PATHS.Figures	= [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Figures/'];
		load([PATHS.REF_data,REF_Dataset,'_reference_v2.mat'],'REF')
		
		%	Store waveform pressure data into REF_temp to save space
		for jj = 1:length(REF)
			REF_temp(jj).ID				= REF(jj).INFO.ID;
			REF_temp(jj).Pressure		= REF(jj).L_Carotid.ONE_CYCLE.P;
			REF_temp(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;
			REF_temp(jj).t				= REF(jj).Asc_Ao.ONE_CYCLE.t;
		end
		
end

REF = REF_temp; clear REF_temp

%	Physiological filter based on the Normotensive and Hypertensive datasets
switch REF_Dataset
	case {'Wk3','pwdb_78','pwdb_4374'}
		addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
		REF = Physiological_BP_Filter(PATHS,REF);
		for jj = 1:length(REF)
			REF(jj).Pressure		= REF(jj).Pressure/Scale_P;
			REF(jj).Pressure_CBP	= REF(jj).Pressure_CBP/Scale_P;
		end
end


%%	Load EST dataset(s)
disp('Estimation model(s): ')

%	1-D aortic model
switch REF_Dataset
	case 'Isra_AoCoarctation'	%Clinical: aortic coarctation
		%	1-D aortic model
		EST_Dataset	= 'Aortic1D_AoCo_estimation';
		switch Scenario
			case 'Waveform'
				PATHS.EST_data	= [PATHS.Root,'Aortic1D/Nektar_outputFiles/Simulations/2019_09_10/Sc1/',EST_method,'/'];
			case 'BP'
				PATHS.EST_data	= [PATHS.Root,'Aortic1D/Nektar_outputFiles/Simulations/2019_09_19/Sc2/',EST_method,'/'];
		end
		load([PATHS.EST_data,EST_Dataset,'.mat'],'EST')	
		disp(EST_Dataset)
		
% 		%	Compare the right studies (not all simulations run)
		for jj = 1:length(EST)
			for kk = jj:length(REF)
				if contains(EST(jj).INFO.ID,REF(kk).ID)
					REF_index(jj) = kk;
					break
				end
			end
		end
		
		for jj = 1:length(REF_index)
			temp(jj) = REF(REF_index(jj));
		end

		REF = temp; clear temp
		for jj = 1:length(REF), disp(REF(jj).ID), disp(EST(jj).INFO.ID), end
		
		%	Store waveform pressure data into EST_temp to save space
		for jj = 1:length(EST)
			EST_Aortic1D(jj).Pressure	= EST(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			EST_Aortic1D(jj).t			= EST(jj).Asc_Ao.ONE_CYCLE.t;
		end
		clear EST
		
	case {'pwdb_6','pwdb_78','pwdb_4374'}	%In silico: 116-artery (78 or all 4374 subjects)
		%	1-D aortic model
		switch Scenario
			case 'Waveform'
				PATHS.EST_data	= [PATHS.Root,'Aortic1D/Nektar_outputFiles/Simulations/2019_08_30/Sc1/'];
			case 'BP'
				PATHS.EST_data	= [PATHS.Root,'Aortic1D/Nektar_outputFiles/Simulations/2019_09_04/Sc2/'];
		end
		EST_Dataset	= ['Aortic1D_',REF_Dataset,'_estimation'];
		load([PATHS.EST_data,EST_Dataset,'.mat'],'EST')
		disp(EST_Dataset)
		
		%	Store waveform pressure data into EST_temp to save space
		for jj = 1:length(EST)
			EST_Aortic1D(jj).ID			= EST(jj).INFO.ID;
			EST_Aortic1D(jj).Pressure	= EST(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			EST_Aortic1D(jj).t			= EST(jj).Asc_Ao.ONE_CYCLE.t;
		end
		clear EST
		
end

%	Wk2 and Wk3 models
switch REF_Dataset
	case {'Isra_AoCoarctation','Sam_Normotensive','Sam_Hypertensive','pwdb_6','pwdb_78','pwdb_4374'}
		%	3-Wk
		switch Scenario
			case 'Waveform'
				EST_Dataset	= [REF_Dataset,'_ref_Wk3_',EST_method,'_Sc1'];	%CV alternative
			case 'BP'
				EST_Dataset	= [REF_Dataset,'_ref_Wk3_',EST_method,'_Sc2'];	%CV alternative
		end
		PATHS.EST_data	= [PATHS.Root,'Windkessel/Output/Wk3/'];
		load([PATHS.EST_data,EST_Dataset,'.mat'],'EST')
		disp(EST_Dataset)
		
		%	Store waveform pressure data into REF_temp and EST_temp to save space
		for jj = 1:length(REF)
			EST_Wk3(jj).ID			= EST(jj).ID;
			EST_Wk3(jj).Pressure	= EST(jj).Pressure/Scale_P;
			EST_Wk3(jj).t			= linspace(0,REF(jj).t(end),length(EST_Wk3(jj).Pressure))';
		end
		clear EST
		
		%	2-Wk
		switch Scenario
			case 'Waveform'
				EST_Dataset	= [REF_Dataset,'_ref_Wk2_',EST_method(1:4),'_Sc1'];	%CV alternative
			case 'BP'
				EST_Dataset	= [REF_Dataset,'_ref_Wk2_',EST_method(1:4),'_Sc2'];	%CV alternative
		end
		PATHS.EST_data	= [PATHS.Root,'Windkessel/Output/Wk2/'];
		load([PATHS.EST_data,EST_Dataset,'.mat'],'EST')
		disp(EST_Dataset)
		
		%	Store waveform pressure data into REF_temp and EST_temp to save space
		for jj = 1:length(REF)
			EST_Wk2(jj).ID			= EST(jj).ID;
			EST_Wk2(jj).Pressure	= EST(jj).Pressure/Scale_P;
			EST_Wk2(jj).t			= linspace(0,REF(jj).t(end),length(EST_Wk2(jj).Pressure))';
		end
		clear EST
end


%%	PLOTS general format
addpath([PATHS.Root,'/Others/PlotFormat'])
addpath([PATHS.Root,'/Others/PlotSave'])

PLOTS.Format				= 1;
PLOTS.BA					= 0;

%	Fontsize
PLOTS.FontSize				= 50;
PLOTS.LegendSize			= 40;
PLOTS.FontSizeAnno			= 40;
PLOTS.LineWidth				= 2;
PLOTS.MarkerSize			= 20;

%	Title and other text
PLOTS.Title				= 0;	% 1: Display title on each figure
PLOTS.Legend			= 1;	% 1: Display Legend_text on each figure
PLOTS.Legend_box		= 0;
PLOTS.XLabel			= 1;	% 1: Display x-axis labels
PLOTS.YLabel			= 1;	% 1: Display y-axis labels
PLOTS.Annotation		= 2;	% 1: Display annotations

%	Axis properties	
PLOTS.XLim				= 2;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.YLim				= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.Grid				= 1;	% 1: grid on

%	Error box
PLOTS.ErrorBox			= 0;	% 1: Display absolute error box
								% 2: Display relative error box


PLOTS.Legend_text		= {'Reference','Estimated'};
PLOTS.XLabel_text		= {'Time [s]'};
PLOTS.YLabel_text		= {'Pressure [mmHg]'};
PLOTS.XLim_values		= [0 200];
PLOTS.YLim_values		= [70 110];
PLOTS.XTick				= 10;
PLOTS.YTick				= 10;
PLOTS.Markers			= {'bx'};


%%	Plot individual waveform comparisons
% for jj = 1%:length(EST_Aortic1D)
% 	figure, hold on
% 	plot(REF(jj).t,REF(jj).Pressure_CBP,'k','LineWidth',6)
% 	plot(EST_Aortic1D(jj).t,EST_Aortic1D(jj).Pressure,'-.k','LineWidth',4)
% 	hold off
% 	FormatPlots(PLOTS)
% % 	PlotSave(PATHS.Figures,[REF_Dataset,'_vs_1-D_aortic_waveform_',num2str(jj)])
% 	PlotSave('/Users/joh15/PhD/PAPER - CBP estimation/Figures/AlgorithmSummary/cBP/',[REF_Dataset,'_vs_1-D_aortic_waveform'])	%for cBP paper
% end
% 
% for jj = 1%:length(EST_Wk3)
% 	figure, hold on
% 	plot(REF(jj).t,REF(jj).Pressure_CBP,'k','LineWidth',6)
% 	plot(EST_Wk3(jj).t,EST_Wk3(jj).Pressure,'-.k','LineWidth',4)
% 	hold off
% 	FormatPlots(PLOTS)
% % 	PlotSave(PATHS.Figures,[REF_Dataset,'_vs_3-Wk_waveform'])
% 	PlotSave('/Users/joh15/PhD/PAPER - CBP estimation/Figures/AlgorithmSummary/cBP/',[REF_Dataset,'_vs_3-Wk_waveform'])	%for cBP paper
% end
% 
% for jj = 1%:length(EST_Wk2)
% 	figure, hold on
% 	plot(REF(jj).t,REF(jj).Pressure_CBP,'k','LineWidth',6)
% 	plot(EST_Wk2(jj).t,EST_Wk2(jj).Pressure,'-.k','LineWidth',4)
% 	hold off
% 	FormatPlots(PLOTS)
% % 	PlotSave(PATHS.Figures,[REF_Dataset,'_vs_2-Wk_waveform'])
% 	PlotSave('/Users/joh15/PhD/PAPER - CBP estimation/Figures/AlgorithmSummary/cBP/',[REF_Dataset,'_vs_2-Wk_waveform'])	%for cBP paper
% end


%%	Plot multiple waveform comparisons
%	Number and arrangement of subplots per figure
PLOTS.Horizontal		= 2;
PLOTS.Vertical			= 5;
PLOTS.FontSize			= 10;

switch REF_Dataset
	case 'pwdb_4374'
		PLOTS.YLim_values		= [60 160];
		PLOTS.YTick				= 20;
		PLOTS.Markers			= {'k--','b--','r--'};
		Representative_plots	= [1:50:4000];
		REF_A4 = REF(Representative_plots);
		EST_A4 = {EST_Aortic1D(Representative_plots),EST_Wk3(Representative_plots),EST_Wk2(Representative_plots)};
	case 'Isra_AoCoarctation'
		PLOTS.YLim_values		= [30 130];
		PLOTS.YTick				= 20;
		PLOTS.Markers			= {'k--','b--','r--'};
		REF_A4 = REF;
		EST_A4 = {EST_Aortic1D,EST_Wk3,EST_Wk2};
	case 'Sam_Normotensive'
		PLOTS.YLim_values		= [30 190];
		PLOTS.YTick				= 50;
		PLOTS.Markers			= {'b--','r--'};
		REF_A4 = REF;
		EST_A4 = {EST_Wk3,EST_Wk2};
	case 'Sam_Hypertensive'
		PLOTS.YLim_values		= [40 210];
		PLOTS.YTick				= 50;
		PLOTS.Markers			= {'b--','r--'};
		REF_A4 = REF;
		EST_A4 = {EST_Wk3,EST_Wk2};
end

addpath([PATHS.Root,'/Others/PlotWaveformComparison_A4/'])
% PlotWaveformComparison_A4_v2(PATHS,PLOTS,REF_A4,EST_A4,REF_Dataset,Sc)
clear REF_A4 EST_A4


%%	Extract DBP, MBP, SBP, PP from REF and EST
addpath([PATHS.Root,'/Others/ExtractReferenceBP'])

for jj = 1:length(REF)
	REF(jj).DBP	= min(REF(jj).Pressure_CBP);	% DBP
	REF(jj).MBP	= mean(REF(jj).Pressure_CBP);	% MBP
	REF(jj).SBP	= max(REF(jj).Pressure_CBP);	% SBP
	REF(jj).PP = REF(jj).SBP - REF(jj).DBP;
	switch REF_Dataset
		case {'Isra_AoCoarctation','pwdb_6','pwdb_78','pwdb_4374'}
			EST_Aortic1D_temp(jj)	= ExtractReferenceBP('Waveform',EST_Aortic1D(jj));
	end
	EST_Wk3_temp(jj)		= ExtractReferenceBP('Waveform',EST_Wk3(jj));
	EST_Wk2_temp(jj)		= ExtractReferenceBP('Waveform',EST_Wk2(jj));
end

switch REF_Dataset
	case {'Isra_AoCoarctation','pwdb_6','pwdb_78','pwdb_4374'}
		EST_Aortic1D	= EST_Aortic1D_temp; clear EST_Aortic1D_temp
end
EST_Wk3			= EST_Wk3_temp; clear EST_Wk3_temp
EST_Wk2			= EST_Wk2_temp; clear EST_Wk2_temp

%	Multiple estimation B-A plots
switch REF_Dataset
	case {'Isra_AoCoarctation','pwdb_6','pwdb_78','pwdb_4374'}
		EST_temp.SBP = {[EST_Aortic1D.SBP],[EST_Wk3.SBP],[EST_Wk2.SBP]};
		EST_temp.MBP = {[EST_Aortic1D.MBP],[EST_Wk3.MBP],[EST_Wk2.MBP]};
		EST_temp.DBP = {[EST_Aortic1D.DBP],[EST_Wk3.DBP],[EST_Wk2.DBP]};
		EST_temp.PP = {[EST_Aortic1D.PP],[EST_Wk3.PP],[EST_Wk2.PP]};
		
% 		EST_temp.SBP = {[EST_Aortic1D.SBP]};
% 		EST_temp.MBP = {[EST_Aortic1D.MBP]};
% 		EST_temp.DBP = {[EST_Aortic1D.DBP]};
% 		EST_temp.PP = {[EST_Aortic1D.PP]};
	case {'Sam_Normotensive','Sam_Hypertensive'}
		EST_temp.SBP = {[EST_Wk2.SBP],[EST_Wk3.SBP]};
		EST_temp.MBP = {[EST_Wk2.MBP],[EST_Wk3.MBP]};
		EST_temp.DBP = {[EST_Wk2.DBP],[EST_Wk3.DBP]};
		EST_temp.PP = {[EST_Wk2.PP],[EST_Wk3.PP]};
end


%%	Calculate and display ERRORS
addpath([PATHS.Root,'/Others/CalculateBPErrors'])
addpath([PATHS.Root,'/Others/P_0_iteration/'])
addpath([PATHS.Root,'/Others/InterpolateSpline/'])
addpath([PATHS.Root,'/Others/Quartiles/'])

switch REF_Dataset
	case {'Isra_AoCoarctation','pwdb_6','pwdb_78','pwdb_4374'}
% 		Model = 'Wk2';
% 		ERRORS = CalculateBPErrors(REF,EST_Wk2);
% 		ERRORS = ERRORS_RMSE(REF,EST_Wk2,ERRORS,0);
% 		DisplayErrors_abs(ERRORS,Model)
% % 		DisplayErrors_rel(ERRORS,Model)
% 		DisplayR2(REF,EST_Wk2);
% % 		Quartiles([[REF.PP];[ERRORS.PP_abs]],3);
% 		
% 		Model = 'Wk3';
% 		ERRORS = CalculateBPErrors(REF,EST_Wk3);
% 		ERRORS = ERRORS_RMSE(REF,EST_Wk3,ERRORS,0);
% 		DisplayErrors_abs(ERRORS,Model)
% % 		DisplayErrors_rel(ERRORS,Model)
% 		DisplayR2(REF,EST_Wk3);
% % 		Quartiles([[REF.PP];[ERRORS.PP_abs]],3);
% 
% 		Model = 'Aortic1D';
		ERRORS = CalculateBPErrors(REF,EST_Aortic1D);
% 		ERRORS = ERRORS_RMSE(REF,EST_Aortic1D,ERRORS,0);
% 		DisplayErrors_abs(ERRORS,Model)
% % 		DisplayErrors_rel(ERRORS,Model)
% 		DisplayR2(REF,EST_Aortic1D);
% % 		Quartiles([[REF.PP];[ERRORS.PP_abs]],3);
		
	case {'Sam_Normotensive','Sam_Hypertensive'}		
		Model = 'Wk2';		
		ERRORS = CalculateBPErrors(REF,EST_Wk2);
		ERRORS = ERRORS_RMSE(REF,EST_Wk2,ERRORS,0);
		DisplayErrors_abs(ERRORS,Model)
% 		DisplayErrors_rel(ERRORS,Model)
		DisplayR2(REF,EST_Wk2);

		Model = 'Wk3';
		ERRORS = CalculateBPErrors(REF,EST_Wk3);
		ERRORS = ERRORS_RMSE(REF,EST_Wk3,ERRORS,0);
		DisplayErrors_abs(ERRORS,Model)
% 		DisplayErrors_rel(ERRORS,Model)
		DisplayR2(REF,EST_Wk3);
end


%%	Temporary code to find patients corresponding to largest BP errors
load('/Users/joh15/Haemodynamic_Tools/Version6/ParameterEstimation/Output/pwdb_4374_carotid/PWV_outliers.mat')
PWV_outliers = PWV_outliers';

% Scenario 1
DBP_outliers = find(abs([ERRORS.DBP_abs]) > 8)';
MBP_outliers = find(abs([ERRORS.MBP_abs]) > 5)';
SBP_outliers = find(abs([ERRORS.SBP_abs]) > 10)';
PP_outliers = find(abs([ERRORS.PP_abs]) > 10)';

% Scenario 2
DBP_outliers = find(abs([ERRORS.DBP_abs]) > 5)';
MBP_outliers = find(abs([ERRORS.MBP_abs]) > 5)';
SBP_outliers = find(([ERRORS.SBP_abs] > 10)|([ERRORS.SBP_abs] < -18))';
PP_outliers = find(([ERRORS.PP_abs] > 12)|([ERRORS.PP_abs] < -20))';

[val_1,pos_1] = intersect(PWV_outliers,DBP_outliers);
[val_2,pos_2] = intersect(PWV_outliers,MBP_outliers);
[val_3,pos_3] = intersect(PWV_outliers,SBP_outliers);
[val_4,pos_4] = intersect(PWV_outliers,PP_outliers);


%%	Bland-Altman PLOTS
addpath([PATHS.Root,'/Others/PlotBlandAltman'])
PLOTS.BA				= 1;
PLOTS.FontSize			= 50;


%	Title, legend, axis labels, and annotations
PLOTS.Title				= 0;	% 1: Display title on each figure
PLOTS.Legend			= 0;	% 1: Display Legend_text on each figure
PLOTS.XLabel			= 0;	% 1: Display x-axis labels
PLOTS.YLabel			= 0;	% 1: Display y-axis labels
PLOTS.Annotation		= 2;	% 1: Display annotations

%	Axis properties	
PLOTS.XLim				= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.YLim				= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.Grid				= 1;	% 1: grid on

%	Error box
PLOTS.ErrorBox			= 0;	% 1: Display absolute error box
								% 2: Display relative error box

%	Text
PLOTS.Legend_text		= {''};
PLOTS.XLabel_text		= {'Reference [mmHg]'};
PLOTS.YLabel_text		= {'Estimate - Reference [mmHg]'};
switch REF_Dataset
	case {'Isra_AoCoarctation','pwdb_6','pwdb_78','pwdb_4374'}
		EST_Name				= {'Aortic1D','Wk3','Wk2'};
% 		EST_Name				= {'Aortic1D'};
		Markers					= {'ko','bx','r^'};
	case {'Sam_Normotensive','Sam_Hypertensive'}
		EST_Name				= {'Wk3','Wk2'};
		Markers					= {'bx','r^'};
end

BP						= {'DBP','MBP','SBP','PP'};
Title_text				= BP;

switch REF_Dataset
	case {'pwdb_6','pwdb_78','pwdb_4374'}
		XLim_values				= {[40 100],	[70 110],	[80 170],	[10 110]};
		XTick					= {20,			10,			20,			20};
		YLim_values				= {[-20 20],	[-10 10],	[-30 30],	[-40 40]};
		YTick					= {10,			5,			10,			20};
	case 'Isra_AoCoarctation'
		XLim_values				= {[30 70],		[50 90],	[70 130],	[20 80]};
		XTick					= {10,			10,			20,			20};
		YLim_values				= {[-20 20],	[-15 15],	[-35 35],	[-40 40]};
		YTick					= {10,			5,			10,			20};
	case 'Sam_Normotensive'
		XLim_values				= {[40 100],	[50 130],	[70 170],	[20 80]};
		XTick					= {20,			20,			30,			20};
		YLim_values				= {[-20 20],	[-8 8],		[-50 50],	[-60 60]};
		YTick					= {10,			5,			20,			20};
	case 'Sam_Hypertensive'
		XLim_values				= {[40 140],	[60 170],	[50 220],	[0 110]};
		XTick					= {20,			20,			50,			20};
		YLim_values				= {[-20 20],	[-8 8],	[-50 50],	[-60 60]};
		YTick					= {10,			5,			20,			20};
end

PLOTS.BA_sf_all			= [1, 1, 1, 1];

for jj = 1:length(BP)
	
	PLOTS.XLim_values	= XLim_values{jj};
	PLOTS.YLim_values	= YLim_values{jj};
	PLOTS.XTick			= XTick{jj};
	PLOTS.YTick			= YTick{jj};
	PLOTS.BP			= BP{jj};
	
	PLOTS.BA_sf			= PLOTS.BA_sf_all(jj);
	
	for kk = 1:length(EST_Name)
		EST_Name_temp	= EST_Name{kk};
		PLOTS.Markers	= {Markers{kk}};
		eval(['EST		= EST_',EST_Name_temp,';'])
		
		BlandAltman_cBP_models(PATHS,PLOTS,REF,EST,[EST_method,'_',REF_Dataset],EST_Name_temp,Sc)
	end
end

clear


%%	Linear correlation PLOTS - OLD
% addpath([PATHS.Root,'/Others/PlotLinearCorrelation')
% 
% %	Title and other text
% PLOTS.Title				= 0;	% 1: Display title on each figure
% PLOTS.Legend			= 0;	% 1: Display Legend_text on each figure
% PLOTS.XLabel			= 1;	% 1: Display x-axis labels
% PLOTS.YLabel			= 1;	% 1: Display y-axis labels
% PLOTS.Annotation		= 2;	% 1: Display annotations
% 
% %	Axis properties	
% PLOTS.XLim			= 1;	% 0: no limits, no axis
% 								% 1: 'custom'
% 								% 2: 'custom tight'
% PLOTS.YLim			= 1;	% 0: no limits, no axis
% 								% 1: 'custom'
% 								% 2: 'custom tight'
% PLOTS.Grid			= 1;	% 1: grid on
% 
% %	Error box
% PLOTS.ErrorBox		= 0;	% 1: Display absolute error box
% 								% 2: Display relative error box
% 
% PLOTS.Legend_text	= {'Physical','Optimised'};
% PLOTS.XLabel_text	= {'Reference (mmHg)'};
% PLOTS.YLabel_text	= {'Estimation (mmHg)'};
% PLOTS.XLim_values	= [0 160];  
% PLOTS.YLim_values	= PLOTS.XLim_values;
% PLOTS.XTick			= 50;
% PLOTS.YTick			= PLOTS.XTick;	% They should be the same
% PLOTS.Markers		= {'ko','bx','r^'};

% PLOTS.BA			= 0;
% PLOTS.Title_text	= 'Systolic CBP';
% % PLOTS.PlotName		= [REF_Dataset,'_',EST_Name,'_SBP'];	%[SBP, MBP, DBP, PP]
% PLOTS.XLabel_text	= {'Reference [mmHg]'};
% PLOTS.YLabel_text	= {'Estimation [mmHg]'};
% PLOTS.XLim_values	= [60 140];
% PLOTS.YLim_values	= [60 140];
% PLOTS.XTick			= 20;
% PLOTS.YTick			= 20;
% 
% %	LINEAR CORRELATION
% R_coeff							= corrcoef(REF_value.Pressure,EST_value.Pressure{1});
% PLOTS.Annotation_text{1}		= ['R = ',num2str(round(R_coeff(1,2),2))];
% PlotLinearCorrelation([REF_temp.SBP],[EST_temp.SBP],PLOTS)
% 
% %	Aortic1D
% PLOTS.Markers		= {'ko'};
% PlotLinearCorrelation([REF_temp.SBP],[EST_temp.SBP(1)],PLOTS)
% PLOTS.PlotName		= [REF_Dataset,'_Aortic1D_SBP'];
% PlotSave(PATHS.Figures,[PLOTS.PlotName,'_corr'])
% 
% %	3-Wk
% PLOTS.Markers		= {'bx'};
% PlotLinearCorrelation([REF_temp.SBP],[EST_temp.SBP(2)],PLOTS)
% PLOTS.PlotName		= [REF_Dataset,'_Wk3_SBP'];
% PlotSave(PATHS.Figures,[PLOTS.PlotName,'_corr'])
% 
% %	2-Wk
% PLOTS.Markers		= {'r^'};
% PlotLinearCorrelation([REF_temp.SBP],[EST_temp.SBP(3)],PLOTS)
% PLOTS.PlotName		= [REF_Dataset,'_Wk2_SBP'];
% PlotSave(PATHS.Figures,[PLOTS.PlotName,'_corr'])


%%%
restoredefaultpath


%%	Function definitions
function DisplayErrors_abs(ERRORS,Model)
%	For tables
DBP_mean_SD = [mean([ERRORS.DBP_abs]), std([ERRORS.DBP_abs])];
MBP_mean_SD = [mean([ERRORS.MBP_abs]), std([ERRORS.MBP_abs])];
SBP_mean_SD = [mean([ERRORS.SBP_abs]), std([ERRORS.SBP_abs])];
PP_mean_SD = [mean([ERRORS.PP_abs]), std([ERRORS.PP_abs])];
RMS_mean_SD = [mean([ERRORS.RMS_wave_abs]), std([ERRORS.RMS_wave_abs])];

fprintf('\nAbsolute errors for: %s (DBP, MBP, SBP, PP, RMSE)\n',Model)
fprintf(['& ',num2str(DBP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(DBP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(MBP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(MBP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(SBP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(SBP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(PP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(PP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(RMS_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(RMS_mean_SD(2),'%.1f'),' '])
disp('\\')

%	For reference
SBP_range	= [min([ERRORS.SBP_abs]), max([ERRORS.SBP_abs])];
MBP_range	= [min([ERRORS.MBP_abs]), max([ERRORS.MBP_abs])];
DBP_range	= [min([ERRORS.DBP_abs]), max([ERRORS.DBP_abs])];
PP_range	= [min([ERRORS.PP_abs]), max([ERRORS.PP_abs])];
RMS_range	= [min([ERRORS.RMS_wave_abs]), max([ERRORS.RMS_wave_abs])];
end

function DisplayErrors_rel(ERRORS,Model)
%	For tables
DBP_mean_SD = [mean([ERRORS.DBP_rel]), std([ERRORS.DBP_rel])];
MBP_mean_SD = [mean([ERRORS.MBP_rel]), std([ERRORS.MBP_rel])];
SBP_mean_SD = [mean([ERRORS.SBP_rel]), std([ERRORS.SBP_rel])];
PP_mean_SD = [mean([ERRORS.PP_rel]), std([ERRORS.PP_rel])];
RMS_mean_SD = [mean([ERRORS.RMS_wave_rel]), std([ERRORS.RMS_wave_rel])];

fprintf('\nRelative errors for: %s (DBP, MBP, SBP, PP, RMSE)\n',Model)
fprintf(['& ',num2str(DBP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(DBP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(MBP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(MBP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(SBP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(SBP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(PP_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(PP_mean_SD(2),'%.1f'),' '])
fprintf(['& ',num2str(RMS_mean_SD(1),'%.1f'),' $\\pm$ ',num2str(RMS_mean_SD(2),'%.1f'),' '])
disp('\\')

%	For reference
SBP_range	= [min([ERRORS.SBP_rel]), max([ERRORS.SBP_rel])];
MBP_range	= [min([ERRORS.MBP_rel]), max([ERRORS.MBP_rel])];
DBP_range	= [min([ERRORS.DBP_rel]), max([ERRORS.DBP_rel])];
PP_range	= [min([ERRORS.PP_rel]), max([ERRORS.PP_rel])];
RMS_range	= [min([ERRORS.RMS_wave_rel]), max([ERRORS.RMS_wave_rel])];
end

function [ERR] = ERRORS_RMSE(REF,EST,ERR,Plots)
%	Calculate pressure errors
for jj = 1:length(EST)
	%	RMS: the square root of the arithmetic mean of the squares of the values
	REF_P = REF(jj).Pressure_CBP;	
	EST_P = EST(jj).Pressure;
	
	[REF_P, EST_P, Time] = InterpolateSpline(REF(jj).t,REF_P,EST_P);
	
	if Plots == 1
		figure, hold on
		plot(Time,REF_P,'k','LineWidth',4), plot(Time,EST_P,'--b','LineWidth',4)
		hold off
		legend('Reference cBP','Estimated cBP')
		xlabel('Time [s]')
		ylabel('Pressure [mmHg]')
		set(gca,'FontSize',30)
	end
	
	P_diff = EST_P - REF_P;
	
	ERR(jj).RMS_wave_abs = sqrt( sum(P_diff.^2) / (length(P_diff)-1) );
	ERR(jj).RMS_wave_rel = ERR(jj).RMS_wave_abs / (max(REF_P) - min(REF_P));
end

end
		
function DisplayR2(REF,EST)
%r is calculated using the covariance as cov(REF,EST)/sqrt(cov(REF)*cov(EST))
[r,p] = corrcoef([REF.SBP],[EST.SBP]);
SBP_r = r(1,2);
SBP_p = p(1,2);
[r,p] = corrcoef([REF.MBP],[EST.MBP]);
MBP_r = r(1,2);
MBP_p = p(1,2);
[r,p] = corrcoef([REF.DBP],[EST.DBP]);
DBP_r = r(1,2);
DBP_p = p(1,2);
[r,p] = corrcoef([REF.PP],[EST.PP]);
PP_r = r(1,2);
PP_p = p(1,2);

disp(['R2 values (DBP, MBP, SBP, PP): ',num2str(SBP_r^2,'%.3f'),' ',num2str(MBP_r^2,'%.3f'),' ',...
	num2str(DBP_r^2,'%.3f'),' ',num2str(PP_r^2,'%.3f')])
disp(['p values (DBP, MBP, SBP, PP): ',num2str(SBP_p,'%.3f'),' ',num2str(MBP_p,'%.3f'),' ',...
	num2str(DBP_p,'%.3f'),' ',num2str(PP_p,'%.3f')])
end

function BlandAltman_cBP_models(PATHS,PLOTS,REF,EST,REF_Dataset,EST_Name,Sc)
PLOTS.PlotName	= [REF_Dataset,'_',EST_Name,'_',PLOTS.BP];
eval(['PlotBlandAltman_v2([REF.',PLOTS.BP,'],[EST.',PLOTS.BP,'],PLOTS);']);
PlotSave(PATHS.Figures,[PLOTS.PlotName,'_B-A_',Sc])
end

