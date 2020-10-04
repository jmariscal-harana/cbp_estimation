%%	Run CBP for multiple patients
format compact;
close all;
clear;
clc;


%%	PATHS
PATHS.REF_data			= '~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';
PATHS.EST_data_1D		= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_08_30/Sc1/';
% PATHS.EST_data_1D		= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_09_04/Sc2/';
PATHS.EST_data_Wk2		= '~/Haemodynamic_Tools/Version6/Windkessel/Output/Wk2/';
PATHS.EST_data_Wk3		= '~/Haemodynamic_Tools/Version6/Windkessel/Output/Wk3/';
PATHS.EST_figure		= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Figures/';
PATHS.Matlab			= '~/Haemodynamic_Tools/Version6/Aortic1D/Matlab/';						%Matlab functions

addpath(PATHS.Matlab)


%%	55_art vs Aortic1D
%	Load REF - 55_art
% Population = '55_art';
% load([PATHS.REF_data,Population,'_reference.mat'],'REF')

%	Load EST - Aortic1D
% Population = 'Aortic1D';
% load([PATHS.EST_data_1D,Population,'.mat'],'EST')
% EST_1D = EST;
% 
% clear EST
% 
% %	Number and arrangement of subplots per figure
% Plots_horizontal = 7;
% Plots_vertical = 10;
% Plots_page = Plots_horizontal*Plots_vertical;
% 
% figure
% kk = 0;
% for jj = 1:length(REF)
% 	if REF(jj).INFO.ID{1}(7:end) ~= EST_1D(jj).INFO.ID{1}(8:end), error('Wrong patient comparison'), end
% 	kk = kk + 1;
% 	
% 	subplot(Plots_vertical,Plots_horizontal,kk);
% 	hold on
% 	plot(REF(jj).Asc_Ao.ONE_CYCLE.t, REF(jj).Asc_Ao.ONE_CYCLE.P, 'k','LineWidth',2)
% 	plot(EST_1D(jj).Asc_Ao.ONE_CYCLE.t, EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25, '--b','LineWidth',1)
% 	hold off
% 	
% 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P)
% 		min(EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25)]);
% 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P) 
% 		max(EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25)]);
% 	YTick = 20;	%[mmHg]
% 	
% 	XMin = 0;
% 	XMax = round(max(REF(jj).Asc_Ao.ONE_CYCLE.t)+0.049,1);
% 	XTick = [XMin, XMax];
% 	
% 	xlim([XMin, XMax])
% 	ylim([YMin, YMax])
% 	set(gca,'XTick',XTick)
% 	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
% 	set(gca,'FontSize',10)

% 	if mod(jj,Plots_horizontal) == 1
% 		ylabel('P [mmHg]')
% 	end
% 	if jj > length(REF) - Plots_horizontal || mod(jj-1,Plots_page) >= Plots_page - Plots_horizontal
% 		xlabel('t [s]')
% 	end

% 	if jj == length(REF)
% 		print('-fillpage',[PATHS.EST_figure,'WaveformComparison55-art_',num2str(ceil(jj/Plots_page))],'-dpdf')
% 	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
% 		print('-fillpage',[PATHS.EST_figure,'WaveformComparison55-art_',num2str(jj/Plots_page)],'-dpdf')
% 		figure
% 		kk = 0;
% 	end	
% end
% 
% clear


%%	55_art vs (2-Wk, 3-Wk and Aortic1D)
% %	Load REF - 55_art
% Population = '55_art';
% load([PATHS.REF_data,Population,'_reference.mat'],'REF')
% 
% %	Load EST - Aortic1D
% Population = 'Aortic1D';
% load([PATHS.EST_data_1D,Population,'.mat'],'EST')
% EST_1D = EST;
% 
% %	Load EST - Wk2
% Population = '55_art_ref_Wk2_est_6515';
% load([PATHS.EST_data_Wk2,Population,'.mat'],'EST')
% EST_Wk2 = EST;
% 
% %	Load EST - Wk3
% Population = '55_art_ref_Wk3_6515_0_5';
% load([PATHS.EST_data_Wk3,Population,'.mat'],'EST')
% EST_Wk3 = EST;
% 
% clear EST
% 
% %	Number and arrangement of subplots per figure
% Plots_horizontal = 2;
% Plots_vertical = 2;
% Plots_page = Plots_horizontal*Plots_vertical;
% 
% figure
% kk = 0;
% for jj = 1:length(REF)
% 	if sum(REF(jj).INFO.ID{1}(7:end) ~= EST_1D(jj).INFO.ID{1}(8:end)) ~= 0 || sum(REF(jj).INFO.ID{1} ~= EST_Wk2(jj).ID) ~= 0 || sum(REF(jj).INFO.ID{1} ~= EST_Wk3(jj).ID) ~= 0
% 		error('Wrong patient comparison')
% 	end
% 	kk = kk + 1;
% 	
% 	subplot(Plots_vertical,Plots_horizontal,kk);
% 	hold on
% 	plot(REF(jj).Asc_Ao.ONE_CYCLE.t, REF(jj).Asc_Ao.ONE_CYCLE.P, 'k','LineWidth',2)
% 	plot(EST_1D(jj).Asc_Ao.ONE_CYCLE.t, EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25, '--b','LineWidth',1)
% 	plot(EST_Wk2(jj).t_in, EST_Wk2(jj).Pressure, '--r','LineWidth',1)
% 	plot(EST_Wk3(jj).t_in, EST_Wk3(jj).Pressure, '--m','LineWidth',1)
% 	hold off
% 	
% % 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P)
% % 		min(EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25)
% % 		min(EST_Wk2(jj).Pressure)
% % 		min(EST_Wk3(jj).Pressure)]);
% % 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P) 
% % 		max(EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25)
% % 		max(EST_Wk2(jj).Pressure)
% % 		max(EST_Wk3(jj).Pressure)]);
% 
% 	YMin = 40;
% 	YMax = 150;
% 	YTick = 40;	%[mmHg]
% 	
% 	XMin = 0;
% 	XMax = round(max(REF(jj).Asc_Ao.ONE_CYCLE.t)+0.049,1);
% 	XMax = 1.2;
% 	XTick = [XMin, XMax];
% 	
% 	xlim([XMin, XMax])
% 	ylim([YMin, YMax])
% 	set(gca,'XTick',XTick)
% 	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
% 	set(gca,'FontSize',20)
% 	
% if mod(jj,Plots_horizontal) == 1
% 	ylabel('P [mmHg]')
% end
% if jj > length(REF) - Plots_horizontal || mod(jj-1,Plots_page) >= Plots_page - Plots_horizontal
% 	xlabel('t [s]')
% end
% 
% 	if jj == length(REF)
% 		print('-fillpage',['~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Figures/55_art_vs_3_models_',num2str(ceil(jj/Plots_page))],'-dpdf')
% 	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
% % 		legend('Reference: 55-artery', '1-D aortic', '2-element Windkessel', '3-element Windkessel')
% 		print('-fillpage',['~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Figures/55_art_vs_3_models_',num2str(jj/Plots_page)],'-dpdf')
% 		figure
% 		kk = 0;
% 	end	
% end
% 
% clear


%%	116_art vs (2-Wk, 3-Wk and Aortic1D)
%	Load REF - 116_art
Dataset = 'pwdb_4374';
load([PATHS.REF_data,Dataset,'_reference_v2_filt.mat'],'REF')

%	Load EST - Aortic1D
Dataset = 'Aortic1D_pwdb_4374_estimation';
load([PATHS.EST_data_1D,Dataset,'.mat'],'EST')
EST_1D = EST;

%	Load EST - Wk2
Dataset = 'pwdb_4374_ref_Wk2_3419_Sc1';
% Dataset = 'pwdb_4374_ref_Wk2_5228_Sc2';
load([PATHS.EST_data_Wk2,Dataset,'.mat'],'EST')
EST_Wk2 = EST;

%	Load EST - Wk3
Dataset = 'pwdb_4374_ref_Wk3_3419_0_17_Sc1';
% Dataset = 'pwdb_4374_ref_Wk3_5228_0_3_Sc2';
load([PATHS.EST_data_Wk3,Dataset,'.mat'],'EST')
EST_Wk3 = EST;

clear EST

%	Number and arrangement of subplots per figure
Plots_horizontal = 2;
Plots_vertical = 2;
Plots_page = Plots_horizontal*Plots_vertical;

figure
kk = 0;
for jj = 1:length(REF)
	kk = kk + 1;
	subplot(Plots_vertical,Plots_horizontal,kk);
	hold on
	plot(REF(jj).Asc_Ao.ONE_CYCLE.t, REF(jj).Asc_Ao.ONE_CYCLE.P, 'k','LineWidth',2)
	plot(EST_1D(jj).Asc_Ao.ONE_CYCLE.t, EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25, '--b','LineWidth',1)
	plot(EST_Wk2(jj).t_in, EST_Wk2(jj).Pressure, '--r','LineWidth',1)
	plot(EST_Wk3(jj).t_in, EST_Wk3(jj).Pressure, '--m','LineWidth',1)
	hold off
	
% 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P)
% 		min(EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25)
% 		min(EST_Wk2(jj).Pressure)
% 		min(EST_Wk3(jj).Pressure)]);
% 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P)
% 		max(EST_1D(jj).Asc_Ao.ONE_CYCLE.P/133.25)
% 		max(EST_Wk2(jj).Pressure)
% 		max(EST_Wk3(jj).Pressure)]);
	
	YMin = 40;
	YMax = 150;
	YTick = 40;	%[mmHg]
	
	XMin = 0;
	XMax = round(max(REF(jj).Asc_Ao.ONE_CYCLE.t)+0.049,1);
	XMax = 1.2;
	XTick = [XMin, XMax];
	
	xlim([XMin, XMax])
	ylim([YMin, YMax])
	set(gca,'XTick',XTick)
	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
	set(gca,'FontSize',20)
	
	if mod(jj,Plots_horizontal) == 1
		ylabel('P [mmHg]')
	end
	if jj > length(REF) - Plots_horizontal || mod(jj-1,Plots_page) >= Plots_page - Plots_horizontal
		xlabel('t [s]')
	end
	
	if jj == length(REF)
		print('-fillpage',['~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Figures/116_art_vs_3_models_',num2str(ceil(jj/Plots_page))],'-dpdf')
	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
% 		legend('Reference: 116-artery', '1-D aortic', '2-element Windkessel', '3-element Windkessel')
		print('-fillpage',['~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Figures/116_art_vs_3_models_',num2str(jj/Plots_page)],'-dpdf')
		figure
		kk = 0;
	end
end

clear


%%	Normotensive vs (2Wk and 3Wk)
% %	Load REF - Normotensive
% REF_Population = 'Sam_Normotensive';
% PATHS.REF_data	= ['~/Clinical_data/',REF_Population,'/'];
% load([PATHS.REF_data,'Normotensive_jmh.mat'],'REF')
% 
% %	Load EST - Wk2
% Windkessel_EST	= 'Wk2';
% Parameters		= '1511';
% 
% EST_Population = [REF_Population,'_ref_',Windkessel_EST,'_',Parameters];
% PATHS.EST_data = ['~/Haemodynamic_Tools/Version6/Windkessel/',Windkessel_EST,'_output/Clinical/',REF_Population,'/'];
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% EST_Wk2 = EST;
% 
% %	Load EST - Wk3
% Windkessel_EST	= 'Wk3';
% Parameters		= '1511_0_3';
% 
% EST_Population = [REF_Population,'_ref_',Windkessel_EST,'_',Parameters];
% PATHS.EST_data = ['~/Haemodynamic_Tools/Version6/Windkessel/',Windkessel_EST,'_output/Clinical/',REF_Population,'/'];
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% EST_Wk3 = EST;
% 
% %	Number and arrangement of subplots per figure
% Plots_horizontal = 7;
% Plots_vertical = 10;
% Plots_page = Plots_horizontal*Plots_vertical;
% 
% figure
% kk = 0;
% for jj = 1:length(REF)
% 	kk = kk + 1;
% 	
% 	subplot(Plots_vertical,Plots_horizontal,kk);
% 	hold on
% 	plot(REF(jj).Asc_Ao.ONE_CYCLE.t, REF(jj).Asc_Ao.ONE_CYCLE.P, 'k','LineWidth',2)
% 	plot(EST_Wk2(jj).t_in, EST_Wk2(jj).Pressure, '--b','LineWidth',1)
% 	plot(EST_Wk3(jj).t_in, EST_Wk3(jj).Pressure, '--r','LineWidth',1)
% 	hold off
% 		
% 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P)
% 		min(EST_Wk2(jj).Pressure)
% 		min(EST_Wk3(jj).Pressure)]);
% 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P) 
% 		max(EST_Wk2(jj).Pressure)
% 		max(EST_Wk3(jj).Pressure)]);
% % 
% 	YMin = 50;
% 	YMax = 170;
% 	YTick = 50;	%[mmHg]
% 	
% 	XMin = 0;
% 	XMax = round(max(REF(jj).Asc_Ao.ONE_CYCLE.t)+0.049,1);
% % 	XMax = 1.2;
% 	XTick = [XMin, XMax];
% % 	XTick = [XMin:0.5:XMax];
% 
% 	xlim([XMin, XMax])
% 	ylim([YMin, YMax])
% 	set(gca,'XTick',XTick)
% 	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
% 	set(gca,'FontSize',10)
% 	
% 	if mod(jj,Plots_horizontal) == 1
% 		ylabel('P [mmHg]')
% 	end
% 	if jj > length(REF) - Plots_horizontal || mod(jj-1,Plots_page) >= Plots_page - Plots_horizontal
% 		xlabel('t [s]')
% 	end
% 	
% 	if jj == length(REF)
% 		print('-fillpage',['/Users/joh15/Clinical_Data/Sam_Normotensive/Figures/Normotensive_vs_Wk_models_',num2str(ceil(jj/Plots_page))],'-dpdf')
% 	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
% % 		legend('Reference: 55-artery', '1-D aortic', '2-element Windkessel', '3-element Windkessel')
% 		print('-fillpage',['/Users/joh15/Clinical_Data/Sam_Normotensive/Figures/Normotensive_vs_Wk_models_',num2str(jj/Plots_page)],'-dpdf')
% 		figure
% 		kk = 0;
% 	end	
% end
% 
% clear


%%	Hypertensive vs (2Wk, 3Wk)
% %	Load REF - Hypertensive
% REF_Population = 'Sam_Hypertensive';
% PATHS.REF_data	= ['~/Clinical_data/',REF_Population,'/'];
% load([PATHS.REF_data,'Hypertensive_jmh.mat'],'REF')
% 
% %	Load EST - Wk2
% Windkessel_EST	= 'Wk2';
% Parameters		= '1511';
% 
% EST_Population = [REF_Population,'_ref_',Windkessel_EST,'_',Parameters];
% PATHS.EST_data = ['~/Haemodynamic_Tools/Version6/Windkessel/',Windkessel_EST,'_output/Clinical/',REF_Population,'/'];
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% EST_Wk2 = EST;
% 
% %	Load EST - Wk3
% Windkessel_EST	= 'Wk3';
% Parameters		= '1511_0_3';
% 
% EST_Population = [REF_Population,'_ref_',Windkessel_EST,'_',Parameters];
% PATHS.EST_data = ['~/Haemodynamic_Tools/Version6/Windkessel/',Windkessel_EST,'_output/Clinical/',REF_Population,'/'];
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% EST_Wk3 = EST;
% 
% %	Number and arrangement of subplots per figure
% Plots_horizontal = 7;
% Plots_vertical = 10;
% Plots_page = Plots_horizontal*Plots_vertical;
% 
% figure
% kk = 0;
% for jj = 1:length(REF)
% 	kk = kk + 1;
% 	
% 	time = linspace(0,REF(jj).Asc_Ao.ONE_CYCLE.t(end),length(REF(jj).Asc_Ao.ONE_CYCLE.P));
% 	
% 	subplot(Plots_vertical,Plots_horizontal,kk);
% 	hold on
% 	plot(time, REF(jj).Asc_Ao.ONE_CYCLE.P, 'k','LineWidth',2)
% 	plot(EST_Wk2(jj).t_in, EST_Wk2(jj).Pressure, '--b','LineWidth',1)
% 	plot(EST_Wk3(jj).t_in, EST_Wk3(jj).Pressure, '--r','LineWidth',1)
% 	hold off
% 		
% 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P)
% 		min(EST_Wk2(jj).Pressure)
% 		min(EST_Wk3(jj).Pressure)]);
% 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P) 
% 		max(EST_Wk2(jj).Pressure)
% 		max(EST_Wk3(jj).Pressure)]);
% % 
% 	YMin = 50;
% 	YMax = 210;
% 	YTick = 50;	%[mmHg]
% 	
% 	XMin = 0;
% 	XMax = round(max(REF(jj).Asc_Ao.ONE_CYCLE.t)+0.049,1);
% % 	XMax = 1.2;
% 	XTick = [XMin, XMax];
% % 	XTick = [XMin:0.5:XMax];
% 	
% 	xlim([XMin, XMax])
% 	ylim([YMin, YMax])
% 	set(gca,'XTick',XTick)
% 	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
% 	set(gca,'FontSize',10)
% 	
% 	if mod(jj,Plots_horizontal) == 1
% 		ylabel('P [mmHg]')
% 	end
% 	if jj > length(REF) - Plots_horizontal || mod(jj-1,Plots_page) >= Plots_page - Plots_horizontal
% 		xlabel('t [s]')
% 	end
% 	
% 	if jj == length(REF)
% 		print('-fillpage',['/Users/joh15/Clinical_Data/Sam_Hypertensive/Figures/Hypertensive_vs_Wk_models_',num2str(ceil(jj/Plots_page))],'-dpdf')
% 	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
% % 		legend('Reference: 55-artery', '1-D aortic', '2-element Windkessel', '3-element Windkessel')
% 		print('-fillpage',['/Users/joh15/Clinical_Data/Sam_Hypertensive/Figures/Hypertensive_vs_Wk_models_',num2str(jj/Plots_page)],'-dpdf')
% 		figure
% 		kk = 0;
% 	end	
% end
% 
% clear


%%	Aortic coarctation vs (Aortic1D, Wk2, Wk3)
% %	Load REF - Hypertensive
% REF_Population = 'Isra_AoCoarctation';
% PATHS.REF_data	= ['~/Clinical_data/',REF_Population,'/'];
% load([PATHS.REF_data,'Isra_AoCoarctation_reference.mat'],'REF')
% 
% Exclude_studies = [5 14 15];
% 
% REF = REF([1:4, 6:13]);
% 
% %	Vectors t and P should have the same length
% for jj = 1:length(REF)
% 	if length(REF(jj).Asc_Ao.ONE_CYCLE.t) ~= length(REF(jj).Asc_Ao.ONE_CYCLE.P)
% 		T	= REF(jj).Asc_Ao.ONE_CYCLE.t(end);
% 		dt	= T / (length(REF(jj).Asc_Ao.ONE_CYCLE.P) - 1);
% 		REF(jj).Asc_Ao.ONE_CYCLE.t = [0 : dt : REF(jj).Asc_Ao.ONE_CYCLE.t(end)]';
% 	end
% end
% 
% %	Load EST - Aortic1D
% PATHS.EST_data_1D		= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_03_21/';
% load([PATHS.EST_data_1D,'Aortic1D_AoCo_estimation.mat'],'EST')
% EST_1D = EST;
% 
% %	Load EST - Wk2
% % Windkessel_EST	= 'Wk2';
% % Parameters		= '1511';
% % 
% % EST_Population = [REF_Population,'_ref_',Windkessel_EST,'_',Parameters];
% % PATHS.EST_data = ['~/Haemodynamic_Tools/Version6/Windkessel/',Windkessel_EST,'_output/Clinical/',REF_Population,'/'];
% % load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% % EST_Wk2 = EST;
% 
% %	Load EST - Wk3
% Windkessel_EST	= 'Wk3';
% Parameters		= '1511_0_3';
% 
% EST_Population = [REF_Population,'_ref_',Windkessel_EST,'_',Parameters];
% PATHS.EST_data = ['~/Haemodynamic_Tools/Version6/Windkessel/',Windkessel_EST,'_output/Clinical/',REF_Population,'/'];
% load([PATHS.EST_data,EST_Population,'.mat'],'EST')
% EST_Wk3 = EST;
% 
% EST_Wk3 = EST_Wk3([1:4, 6:13]);
% 
% %	Number and arrangement of subplots per figure
% Plots_horizontal = 3;
% Plots_vertical = 5;
% Plots_page = Plots_horizontal*Plots_vertical;
% 
% figure
% kk = 0;
% 
% for jj = 1:length(REF)
% 	kk = kk + 1;
% 	
% 	subplot(Plots_vertical,Plots_horizontal,kk);
% 	hold on
% 	plot(REF(jj).Asc_Ao.ONE_CYCLE.t, REF(jj).Asc_Ao.ONE_CYCLE.P/133.25, 'k','LineWidth',2)
% 	plot(EST_1D(jj).Asc_Ao.ONE_CYCLE.t, EST_1D(jj).Asc_Ao.ONE_CYCLE.P, '--b','LineWidth',1)
% % 	plot(EST_Wk2(jj).t_in, EST_Wk2(jj).Pressure, '--g','LineWidth',1)
% 	plot(EST_Wk3(jj).t_in, EST_Wk3(jj).Pressure, '--r','LineWidth',1)
% 	hold off
% 		
% 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P/133.25)
% 		min(EST_1D(jj).Asc_Ao.ONE_CYCLE.P);
% % 		min(EST_Wk2(jj).Pressure)
% 		min(EST_Wk3(jj).Pressure)]);
% 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P/133.25)
% 		max(EST_1D(jj).Asc_Ao.ONE_CYCLE.P);
% % 		max(EST_Wk2(jj).Pressure)
% 		max(EST_Wk3(jj).Pressure)]);
% 
% 	YMin = 30;
% 	YMax = 130;
% 	YTick = 25;	%[mmHg]
% 	
% 	XMin = 0;
% 	XMax = round(max(REF(jj).Asc_Ao.ONE_CYCLE.t)+0.049,1);
% % 	XMax = 1.2;
% 	XTick = [XMin, XMax];
% % 	XTick = [XMin:0.5:XMax];
% 
% 	xlim([XMin, XMax])
% 	ylim([YMin, YMax])
% 	set(gca,'XTick',XTick)
% 	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
% 	set(gca,'FontSize',16)
% 
% if mod(jj,Plots_horizontal) == 1
% 	ylabel('P [mmHg]')
% end
% if mod (jj,Plots_page) > Plots_page - Plots_horizontal
% 	xlabel('t [s]')
% end
%
% 	if jj == length(REF)
% 		print('-fillpage',['/Users/joh15/Clinical_Data/Isra_AoCoarctation/Figures/AoCoarctation_vs_3_models_',num2str(ceil(jj/Plots_page))],'-dpdf')
% 	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
% % 		legend('Reference: 55-artery', '1-D aortic', '2-element Windkessel', '3-element Windkessel')
% 		print('-fillpage',['/Users/joh15/Clinical_Data/Isra_AoCoarctation/Figures/AoCoarctation_vs_3_models_',num2str(jj/Plots_page)],'-dpdf')
% 		figure
% 		kk = 0;
% 	end	
% end
% 
% clear
