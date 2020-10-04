%%	Plot waveforms at different locations for a given subject
%	Input
%	-REF: parameters used to generate virtual patients
%	-EST: methods for parameter estimation for virtual patients
%	-PATHS: folders where required data is stored
%
%	Output
%	-REF: virtual patient data
%	-EST: estimated CV parameters
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (15/09/19) - based on my (lack of) imagination
%
%==========================================================================
format compact;
close all;
clear;
clc;


%%	Load data
PathRoot = [cd,'/'];
PathData = '~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';
NameData = 'pwdb_data_w_aorta_foot_path_p';
load([PathData,NameData])
Data_aorta = data;
clear data
NameData = 'pwdb_data_w_aorta_brain_path';
load([PathData,NameData])
Data_brain = data;
clear data


%%	Plot pressure every ~5 cm, starting at 0 for the first plot and adding 0.1*MBP from previous waveform
%	Aortic path (no shift)
Subject_aorta = Data_aorta.path_waves.aorta_foot(3457);

figure, hold on
MBP = 0;
DBP = min(Subject_aorta.P{1});
k = 0;
% plot([0.196 0.18],[52.66 152],'-.r','LineWidth',3)
plot([0.289 0.294],[39.21, 121.7],'--r','LineWidth',3)
for j = [1,3,9,12,16,18,21,29,32,34,40]
	time = linspace(0,0.002*(length(Subject_aorta.P{j})-1),length(Subject_aorta.P{j}));
	plot(time,Subject_aorta.P{j} - DBP + 0.1*k*MBP,'k','LineWidth',2)
	MBP = mean(Subject_aorta.P{j});
	k = k+1;
end
legend('Dicrotic notch')
xlim([0,0.65])
set(gca,'XTick',0:0.2:0.6)
set(gca,'YTick',[0:2*8.6:86])
yticklabels({'0','5','10','15','20','25','30','35','40','45','50'})
yticklabels({'0','10','20','30','40','50'})
ylabel('Distance from aortic root [cm]')
xlabel('Time [s]')
set(gca,'FontSize',30)
addpath('~/Haemodynamic_Tools/Version6/Others/PlotSave/')
PlotSave(PathRoot,'Subject_arterial_path_aorta')

%	Aortic path
Subject_aorta = Data_aorta.path_waves.aorta_foot(3457);

figure, hold on
MBP = 0;
DBP = min(Subject_aorta.P{1});
k = 0;
plot([0.196 0.18],[52.66 152],'-.r','LineWidth',3)
plot([0.29 0.342],[39.21, 121.7],'--r','LineWidth',3)
for j = [1,3,9,12,16,18,21,29,32,34,40]
	time = Subject_aorta.onset_time(j) + linspace(0,0.002*(length(Subject_aorta.P{j})-1),length(Subject_aorta.P{j}));
	plot(time,Subject_aorta.P{j} - DBP + 0.1*k*MBP,'k','LineWidth',2)
	MBP = mean(Subject_aorta.P{j});
	k = k+1;
end
legend('Systolic peak','Dicrotic notch')
xlim([0,0.65])
set(gca,'XTick',0:0.2:0.6)
set(gca,'YTick',[0:2*8.6:86])
yticklabels({'0','5','10','15','20','25','30','35','40','45','50'})
yticklabels({'0','10','20','30','40','50'})
ylabel('Distance from aortic root [cm]')
xlabel('Time [s]')
set(gca,'FontSize',30)
addpath('~/Haemodynamic_Tools/Version6/Others/PlotSave/')
PlotSave(PathRoot,'Subject_arterial_path_aorta')

%	Aortic root - carotid path
Subject_brain = Data_brain.path_waves.aorta_brain(1);

figure, hold on
MBP = 0;
DBP = min(Subject_brain.P{1});
k = 0;
% plot([0.196 0.18],[52.66 152],'-.r','LineWidth',3)
% plot([0.29 0.342],[39.21, 121.7],'--r','LineWidth',3)
for j = [1,3,7,10,12,15,17,20]
	time = Subject_brain.onset_time(j) + linspace(0,0.002*(length(Subject_brain.P{j})-1),length(Subject_brain.P{j}));
	plot(time,Subject_brain.P{j} - DBP + 0.1*k*MBP,'k','LineWidth',2)
	MBP = mean(Subject_brain.P{j});
	k = k+1;
end
% legend('Systolic peak','Dicrotic notch')
xlim([0,0.65])
set(gca,'XTick',0:0.2:0.6)
set(gca,'YTick',[0:2*8.6:86])
yticklabels({'0','5','10','15','20','25','30','35','40','45','50'})
yticklabels({'0','10','20','30','40','50'})
ylabel('Distance from aortic root [cm]')
xlabel('Time [s]')
set(gca,'FontSize',30)
addpath('~/Haemodynamic_Tools/Version6/Others/PlotSave/')
PlotSave(PathRoot,'Subject_arterial_path_carotid')


