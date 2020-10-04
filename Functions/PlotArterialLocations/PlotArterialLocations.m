%%	Plot waveforms at different locations for the same patient
format compact;
close all;
clear;
clc;


%%	PATHS
PATHS.REF_data			= '~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_12_06/';


%%	Load REF - 55_art
Population = '55_art';
load([PATHS.REF_data,Population,'_reference.mat'],'REF')


%%	Select .HIS to plot
Artery_name					= {'Asc_Ao'; 'Desc_Ao'}';

%	Patients whose R_T (aortic root) < R_T_Wk (parallel sum of Windkessel resistances)
Patient_ID = find([REF.R_T]-[REF.R_T_Wk] < 0);


%%	Plot all waveforms simultaneously, accounting for time-delay among locations
for jj = Patient_ID
	
% 	figure
	hold on
	for kk = 1:length(Artery_name)
		eval(['x_0 = REF(jj).',Artery_name{kk},'.ONE_CYCLE.start - REF(jj).',Artery_name{1},'.ONE_CYCLE.start;']);
		eval(['if x_0 < 0, x_0 = REF(jj).',Artery_name{kk},'.ONE_CYCLE.start - REF(jj).',Artery_name{1},'.ONE_CYCLE.start + length(REF(jj).',Artery_name{kk},'.ONE_CYCLE.P);, end'])
		eval(['x_end = x_0 + length(REF(jj).',Artery_name{kk},'.ONE_CYCLE.P) - 1;'])
		eval(['plot([x_0:x_end],REF(jj).',Artery_name{kk},'.ONE_CYCLE.P)'])
		Artery_name_plot{kk} = Artery_name{kk};
		Artery_name_plot{kk}(strfind(Artery_name{kk},'_')) = '-';
	end
% 	hold off
	legend(Artery_name_plot)
	set(gca,'FontSize',20)
	
end