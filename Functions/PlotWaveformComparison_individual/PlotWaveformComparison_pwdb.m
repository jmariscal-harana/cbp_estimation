close all;
clear;
clc;


%%	Scaling factors for pressure, flow and area
% Scale_P = 133.32;
% Scale_U = 1;
% Scale_A = 1;
% Scale_Q = 1e6;
Scale_P = 1;
Scale_U = 1;
Scale_A = 1;
Scale_Q = 1e6;


%%	Load REF dataset
PATHS.Root		= '~/Haemodynamic_Tools/Version6/';
REF_Dataset		= 'pwdb_4374';	%'Isra_AoCoarctation', 'Sam_Hypertensive', 'Sam_Normotensive', 'pwdb_78', 'pwdb_4374'

%%%
disp('Reference dataset: ')
disp(REF_Dataset)

switch REF_Dataset
	case 'Isra_AoCoarctation'	%Clinical: aortic coarctation
		PATHS.REF_data	= '~/Clinical_Data/Isra_AoCoarctation/';
		PATHS.Figures	= [PATHS.REF_data,'Figures/'];
		load([PATHS.REF_data,REF_Dataset,'_reference.mat'],'REF')
		Scale_P_filt = 1;
	
		%	Store waveform pressure data into REF_temp to save space
		for jj = 1:length(REF)
			REF_temp(jj).Pressure	= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).t			= [0: REF(jj).HAEMO.T/(length(REF_temp(jj).Pressure) - 1) :REF(jj).HAEMO.T]';
		end
		
	case {'Sam_Normotensive','Sam_Hypertensive'}	%Clinical: normotensive or hypertensive
		PATHS.REF_data	= ['~/Clinical_data/',REF_Dataset,'/'];
		PATHS.Figures	= [PATHS.REF_data,'Figures/'];
		load([PATHS.REF_data,REF_Dataset(5:end),'_jmh.mat'],'REF')
		Scale_P_filt = 1;

		%	Store waveform pressure data into REF_temp to save space
		for jj = 1:length(REF)
			REF_temp(jj).ID				= REF(jj).INFO.ID{1};
			REF_temp(jj).Pressure		= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;			
			REF_temp(jj).t				= [0: REF(jj).HAEMO.T/(length(REF_temp(jj).Pressure) - 1) :REF(jj).HAEMO.T];
		end
		
	case {'pwdb_78','pwdb_4374'}	%In silico: 116-artery (78 or all 4374 subjects)
		PATHS.REF_data	= [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Simulations/'];
		PATHS.Figures	= [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Figures/'];
		load([PATHS.REF_data,REF_Dataset,'_reference.mat'],'REF')
		Scale_P_filt = 1;
		
		%	Store waveform pressure data into REF_temp to save space
		for jj = 1:length(REF)
			REF_temp(jj).ID				= REF(jj).INFO.ID;
			REF_temp(jj).Pressure		= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P/Scale_P;
			REF_temp(jj).t				= REF(jj).Asc_Ao.ONE_CYCLE.t;
		end
		
end

REF = REF_temp; clear REF_temp

%	Physiological filter based on the Normotensive and Hypertensive datasets
switch REF_Dataset
	case {'Isra_AoCoarctation','pwdb_78','pwdb_4374'}
		addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
		REF = Physiological_BP_Filter(PATHS,REF,Scale_P_filt);
end


%%	Load EST dataset(s)
% %	1-D aortic
% PATHS.EST_data	= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_05_09/Theoretical_6_subjects/';
% EST_Population	= 'Aortic1D_pwdb_6_estimation';
% 
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% 
% %	Store waveform pressure data into REF_temp and EST_temp to save space
% for jj = 1:length(EST)
% 	EST_temp_1(jj).Pressure		= EST(jj).Asc_Ao.ONE_CYCLE.P;
% 	EST_temp_1(jj).Flow			= EST(jj).Asc_Ao.ONE_CYCLE.Q;
% 	EST_temp_1(jj).Area			= EST(jj).Asc_Ao.ONE_CYCLE.A;
% 	EST_temp_1(jj).t			= EST(jj).Asc_Ao.ONE_CYCLE.t;
% end
% clear EST
% 
% %	1-D aortic
% PATHS.EST_data	= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_05_10/Viscoelastic_6_subjects_opt/';
% EST_Population	= 'Aortic1D_pwdb_6_estimation';
% 
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% 
% %	Store waveform pressure data into REF_temp and EST_temp to save space
% for jj = 1:length(EST)
% 	EST_temp_1_visc(jj).Pressure		= EST(jj).Asc_Ao.ONE_CYCLE.P;
% 	EST_temp_1_visc(jj).Flow			= EST(jj).Asc_Ao.ONE_CYCLE.Q;
% 	EST_temp_1_visc(jj).Area			= EST(jj).Asc_Ao.ONE_CYCLE.A;
% 	EST_temp_1_visc(jj).t				= EST(jj).Asc_Ao.ONE_CYCLE.t;
% end
% clear EST

% %	1-D aortic
% PATHS.EST_data	= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_06_05/';
% EST_Population	= 'Aortic1D_pwdb_4374_estimation';
% 
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% 
% %	Store waveform pressure data into REF_temp and EST_temp to save space
% for jj = 1:length(EST)
% 	EST_temp_1(jj).Pressure		= EST(jj).Asc_Ao.ONE_CYCLE.P;
% 	EST_temp_1(jj).Flow			= EST(jj).Asc_Ao.ONE_CYCLE.Q;
% 	EST_temp_1(jj).Area			= EST(jj).Asc_Ao.ONE_CYCLE.A;
% 	EST_temp_1(jj).t			= EST(jj).Asc_Ao.ONE_CYCLE.t;
% end
% clear EST

% 	Visualise REF vs EST
% for jj = 1:length(REF_temp)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k')
% 	plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure)
% 	hold off
% end
% 
% for jj = 1:length(REF_temp)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Area,'k')
% 	plot(EST_temp_1(jj).t,EST_temp_1(jj).Area)
% 	hold off
% end


%	3-Wk
PATHS.EST_data	= [PATHS.Root,'Windkessel/Wk3_output/Clinical/',REF_Dataset,'/'];
EST_Dataset	= [REF_Dataset,'_ref_Wk3_3311_0_3'];
load([PATHS.EST_data,EST_Dataset,'.mat'],'EST')
		
%	Store waveform pressure data into REF_temp and EST_temp to save space
for jj = 1:length(REF)
	EST_temp_2(jj).Pressure		= EST(jj).Pressure;
	EST_temp_2(jj).t			= EST(jj).t_in;
end
clear EST

% 	Visualise REF vs EST
% for jj = 1:length(REF_temp)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k')
% 	plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure)
% 	hold off
% end


%	2-Wk
PATHS.EST_data	= [PATHS.Root,'Windkessel/Wk2_output/Clinical/',REF_Dataset,'/'];
EST_Dataset	= [REF_Dataset,'_ref_Wk2_3311'];
load([PATHS.EST_data,EST_Dataset,'.mat'],'EST')

%	Store waveform pressure data into REF_temp and EST_temp to save space
for jj = 1:length(REF)
	EST_temp_3(jj).Pressure		= EST(jj).Pressure;
	EST_temp_3(jj).t			= EST(jj).t_in;
end
clear EST

% 	Visualise REF vs EST
% for jj = 1:length(REF_temp)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k')
% 	plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure)
% 	hold off
% end


%%	Plot waveforms
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotFormat')
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotSave')

PLOTS.Format				= 1;

%	Fontsize
PLOTS.FontSize				= 40;
PLOTS.LegendSize			= 30;
PLOTS.LineWidth				= 2;
PLOTS.MarkerSize			= 10;

%	Title and other text
PLOTS.Title				= 0;	% 1: Display title on each figure
PLOTS.Legend			= 1;	% 1: Display Legend_text on each figure
PLOTS.Legend_box		= 0;	% 1: Legend box on
PLOTS.XLabel			= 1;	% 1: Display x-axis labels
PLOTS.YLabel			= 1;	% 1: Display y-axis labels
PLOTS.Annotation		= 0;	% 1: Display annotations

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

PLOTS.Legend_text		= {'cBP: reference','cBP: estimation'};
PLOTS.XLabel_text		= {'Time [s]'};
PLOTS.YLabel_text		= {'Pressure [mmHg]'};
PLOTS.XLim_values		= [0,	1];
PLOTS.YLim_values		= [60,	130];
PLOTS.XTick				= 0.2;
PLOTS.YTick				= 20;	
PLOTS.Markers			= {'bx'};


%%	Plot individual estimations for the same subject
% Subject = 68;	%Choose pwdb patient
% jj = Subject;
% 
% %	pwdb vs 2-Wk
% figure, hold on
% plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
% plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure,'--g','LineWidth',6)
% if PLOTS.Format == 1
% 	FormatPlots(PLOTS)
% end
% PlotSave(PATHS.Figures,[REF_Population,'_vs_2-Wk_waveform_',num2str(jj)])
% 
% %	pwdb vs 3-Wk
% figure, hold on
% plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
% plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure,'--r','LineWidth',6)
% if PLOTS.Format == 1
% 	FormatPlots(PLOTS)
% end
% PlotSave(PATHS.Figures,[REF_Population,'_vs_3-Wk_waveform_',num2str(jj)])
% 
% %	pwdb vs 1-D aortic
% % %	flow
% % PLOTS.YLim_values		= [-50, 350];
% % PLOTS.YTick				= 100;	
% % PLOTS.YLabel_text		= {'Flow [mL/s]'};
% % PLOTS.Legend			= 0;	% 1: Display Legend_text on each figure
% % 
% % figure, hold on
% % plot(REF_temp(jj).t,REF_temp(jj).Flow,'k','LineWidth',6)
% % if PLOTS.Format == 1
% % 	FormatPlots(PLOTS)
% % end
% % PlotSave(PATHS.Figures,[REF_Population,'_vs_1-D_aortic_waveform_flow_',num2str(jj)])	
% 
% %	pressure
% figure, hold on
% plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
% plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure,'--b','LineWidth',6)
% if PLOTS.Format == 1
% 	FormatPlots(PLOTS)
% end
% PlotSave(PATHS.Figures,[REF_Population,'_vs_1-D_aortic_waveform_',num2str(jj)])


%%	Plot all estimations for the same subject
% Subjects = [13 25 44 68];	%Choose pwdb patient(s)
% 
% PLOTS.Legend			= 0;	% 1: Display Legend_text on each figure
% PLOTS.Legend_box		= 1;	% 1: Legend box on
% PLOTS.Legend_text		= {'Reference cBP','Estimation: 1-D aortic','Estimation: 2-el Windkessel','Estimation: 3-el Windkessel'};
% PLOTS.XLabel			= 0;	% 1: Display x-axis labels
% PLOTS.YLabel			= 0;	% 1: Display y-axis labels
% 
% for jj = Subjects
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
% 	plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure,'--b','LineWidth',6)
% 	plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure,'--g','LineWidth',6)
% 	plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure,'--r','LineWidth',6)
% 	if PLOTS.Format == 1
% 		FormatPlots(PLOTS)
% 	end
% 	PlotSave(PATHS.Figures,[REF_Population,'_vs_all_waveforms_',num2str(jj)])
% end


%%	Plot individual estimations for different subjects
% 	Visualise REF vs EST
PLOTS.Legend			= 1;	% 1: Display Legend_text on each figure
PLOTS.Legend_box		= 1;	% 1: Legend box on
PLOTS.Legend_text		= {'REF: SphygmoCor','EST: 3-el Windkessel','EST: 2-el Windkessel'};
% PLOTS.Legend_text		= {'REF: 116-artery','EST: 1-D aortic','EST: 3-el Windkessel','EST: 2-el Windkessel'};
PLOTS.XLabel			= 0;	% 1: Display x-axis labels
PLOTS.YLabel			= 0;	% 1: Display y-axis labels

%	Pressure
PLOTS.YLabel_text		= {'Pressure [Pa]'};
PLOTS.XLim_values		= [0,	1];
PLOTS.YLim_values		= [60,	140];
PLOTS.XTick				= 0.2;
PLOTS.YTick				= 20;	
for jj = round(rand(1,5)*length(REF))%1:length(EST_temp_1)
	figure, hold on
	plot(REF(jj).t,REF(jj).Pressure/Scale_P,'k','LineWidth',4)
% 	plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure/Scale_P,'--k','LineWidth',3)
	plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure/Scale_P,'--b','LineWidth',3)
	plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure/Scale_P,'--r','LineWidth',3)
% 	plot(EST_temp_1_visc(jj).t,EST_temp_1_visc(jj).Pressure/Scale_P,'--r','LineWidth',2)
	hold off
	
	if PLOTS.Format == 1
		FormatPlots(PLOTS)
	end
	PlotSave(PATHS.Figures,[REF_Dataset,'_vs_Wk_P_',num2str(jj)])
end

%	Area
PLOTS.YLabel_text		= {'Area [cm^2]'};
PLOTS.XLim_values		= [0,	0.85];
PLOTS.YLim_values		= [10,	16];
PLOTS.XTick				= 0.2;
PLOTS.YTick				= 2;	
for jj = 1:length(EST_temp_1)
	figure, hold on
	plot(REF_temp(jj).t,REF_temp(jj).Area*10^4,'k','LineWidth',2)
	plot(EST_temp_1(jj).t,EST_temp_1(jj).Area*10^4,'--b','LineWidth',2)
	plot(EST_temp_1_visc(jj).t,EST_temp_1_visc(jj).Area*10^4,'--r','LineWidth',2)
	hold off
	
	if PLOTS.Format == 1
		FormatPlots(PLOTS)
	end
	PlotSave(PATHS.Figures,[REF_Population,'_vs_1-D_aortic_waveform_A_',num2str(jj)])
end

%	Flow
PLOTS.YLabel_text		= {'Flow [mL/s]'};
PLOTS.XLim_values		= [0,	0.85];
PLOTS.YLim_values		= [-75,	450];
PLOTS.XTick				= 0.2;
PLOTS.YTick				= 100;	
for jj = 1:length(EST_temp_1)
	figure, hold on
	plot(REF_temp(jj).t,REF_temp(jj).Flow*10^6,'k','LineWidth',2)
	plot(EST_temp_1(jj).t,EST_temp_1(jj).Flow*10^6,'--b','LineWidth',2)
	plot(EST_temp_1_visc(jj).t,EST_temp_1_visc(jj).Flow*10^6,'--r','LineWidth',2)
	hold off
	
	if PLOTS.Format == 1
		FormatPlots(PLOTS)
	end
	PlotSave(PATHS.Figures,[REF_Population,'_vs_1-D_aortic_waveform_Q_',num2str(jj)])
end