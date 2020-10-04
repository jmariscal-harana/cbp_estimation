close all;
clear;
clc;


%%	Load REF dataset
%	Aortic coarctation
PATHS.REF_data	= '/Users/joh15/Clinical_Data/Isra_AoCoarctation/';
PATHS.Figures	= [PATHS.REF_data,'Figures/'];
REF_Population	= 'Isra_AoCoarctation';
load([PATHS.REF_data,REF_Population,'_reference.mat'],'REF')

%	Store waveform pressure data into REF_temp to save space
for jj = 1:length(REF)
	REF_temp(jj).Pressure	= REF(jj).Asc_Ao.ONE_CYCLE.P/133.322365;
	REF_temp(jj).t			= [0: REF(jj).HAEMO.T/(length(REF_temp(jj).Pressure) - 1) :REF(jj).HAEMO.T]';
end


%%	Load EST dataset(s)
%	1-D aortic
PATHS.EST_data	= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_03_21/';
EST_Population	= 'Aortic1D_AoCo_estimation';
load([PATHS.EST_data,EST_Population,'.mat'],'EST')

%	Compare the right studies (not all simulations run)
for jj = 1:length(EST)
	for kk = jj:length(REF)
		if contains(EST(jj).INFO.ID,REF(kk).INFO.ID)
			REF_index(jj) = kk;
			break
		end
	end
end

REF_temp = REF_temp(REF_index);

%	Store waveform pressure data into REF_temp and EST_temp to save space
for jj = 1:length(EST)
	EST_temp_1(jj).Pressure		= EST(jj).Asc_Ao.ONE_CYCLE.P;
	EST_temp_1(jj).t			= EST(jj).Asc_Ao.ONE_CYCLE.t;
end
clear EST

% 	Visualise REF vs EST
for jj = 1:length(EST_temp_1)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k')
% 	plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure)
% 	hold off
end


%	3-Wk
PATHS.EST_data	= '~/Haemodynamic_Tools/Version6/Windkessel/Wk3_output/Clinical/Isra_AoCoarctation/';
EST_Population	= 'Wk3';
Parameters		= '1511_0_3';
load([PATHS.EST_data,REF_Population,'_ref_',EST_Population,'_',Parameters,'.mat'],'EST')

%	Compare the right studies (not all simulations run)
EST = EST(REF_index);

%	Store waveform pressure data into REF_temp and EST_temp to save space
for jj = 1:length(EST)
	EST_temp_2(jj).Pressure		= EST(jj).Pressure;
	EST_temp_2(jj).t			= EST(jj).t_in;
end
clear EST

% 	Visualise REF vs EST
for jj = 1:length(EST_temp_2)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k')
% 	plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure)
% 	hold off
end


%	2-Wk
PATHS.EST_data	= '~/Haemodynamic_Tools/Version6/Windkessel/Wk2_output/Clinical/Isra_AoCoarctation/';
EST_Population	= 'Wk2';
Parameters		= '1511';
load([PATHS.EST_data,REF_Population,'_ref_',EST_Population,'_',Parameters,'.mat'],'EST')

%	Compare the right studies (not all simulations run)
EST = EST(REF_index);

%	Store waveform pressure data into REF_temp and EST_temp to save space
for jj = 1:length(EST)
	EST_temp_3(jj).Pressure		= EST(jj).Pressure;
	EST_temp_3(jj).t			= EST(jj).t_in;
end

% 	Visualise REF vs EST
for jj = 1:length(EST_temp_3)
% 	figure, hold on
% 	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k')
% 	plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure)
% 	hold off
end

clear REF EST


%%	Plot parameters
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotFormat')
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotSave')

PLOTS.Format				= 1;

%	Fontsize
PLOTS.FontSize				= 40;
PLOTS.LegendSize			= 50;
PLOTS.LineWidth				= 2;
PLOTS.MarkerSize			= 10;

%	Title and other text
PLOTS.Title				= 0;	% 1: Display title on each figure
PLOTS.Legend			= 1;	% 1: Display Legend_text on each figure
PLOTS.XLabel			= 1;	% 1: Display x-axis labels
PLOTS.YLabel			= 1;	% 1: Display y-axis labels
PLOTS.Annotation		= 0;	% 1: Display annotations

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

PLOTS.Legend_text		= {'cBP: invasive','cBP: estimation'};
PLOTS.XLabel_text		= {'Time [s]'};
PLOTS.YLabel_text		= {'Pressure [mmHg]'};
PLOTS.XLim_values		= [0 1];
PLOTS.YLim_values		= [40,130];
PLOTS.XTick				= 50;
PLOTS.YTick				= 20;	% They should be the same
PLOTS.Markers			= {'bx'};


%%	Plot individual estimations for the same subject
jj = 4;	%Choose CoAo patient

%	Aortic coarctation vs 1-D aortic
figure, hold on
plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure,'--b','LineWidth',6)
if PLOTS.Format == 1
	FormatPlots(PLOTS)
end
PlotSave(PATHS.Figures,[REF_Population,'_vs_1-D_aortic_waveform'])

%	Aortic coarctation vs 3-Wk
figure, hold on
plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure,'--r','LineWidth',6)
if PLOTS.Format == 1
	FormatPlots(PLOTS)
end
PlotSave(PATHS.Figures,[REF_Population,'_vs_3-Wk_waveform'])

%	Aortic coarctation vs 2-Wk
figure, hold on
plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure,'--g','LineWidth',6)
if PLOTS.Format == 1
	FormatPlots(PLOTS)
end
PlotSave(PATHS.Figures,[REF_Population,'_vs_2-Wk_waveform'])


%%	Plot all estimations for the same subject
Subjects = [4 6 9 10];	%Choose pwdb patient(s)

PLOTS.Legend			= 0;	% 1: Display Legend_text on each figure
PLOTS.Legend_box		= 1;	% 1: Legend box on
PLOTS.Legend_text		= {'Reference cBP','Estimation: 1-D aortic','Estimation: 2-el Windkessel','Estimation: 3-el Windkessel'};
PLOTS.XLabel			= 0;	% 1: Display x-axis labels
PLOTS.YLabel			= 0;	% 1: Display y-axis labels

for jj = Subjects
	figure, hold on
	plot(REF_temp(jj).t,REF_temp(jj).Pressure,'k','LineWidth',6)
	plot(EST_temp_1(jj).t,EST_temp_1(jj).Pressure,'--b','LineWidth',6)
	plot(EST_temp_2(jj).t,EST_temp_2(jj).Pressure,'--r','LineWidth',6)
	plot(EST_temp_3(jj).t,EST_temp_3(jj).Pressure,'--g','LineWidth',6)
	if PLOTS.Format == 1
		FormatPlots(PLOTS)
	end
	PlotSave(PATHS.Figures,[REF_Population,'_vs_all_waveforms_',num2str(jj)])
end
