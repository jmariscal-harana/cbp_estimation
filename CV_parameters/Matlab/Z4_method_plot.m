%%	Plot pressure and flow waveforms to be used to illustrate the method Z4 for Z_0 estimation
%	Input
%	-REF: parameters used to generate virtual patients
%
%	Output
%	-REF: virtual patient data
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (17/09/19) - based on (mainly) despair
%
%==========================================================================
format compact;
close all;
clear;
clc;


%%	Load data
PathRoot = [cd,'/'];
PathData = '~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';
NameData = 'pwdb_data_w_aorta_finger_path';
load([PathData,NameData])
Data_waves = data.path_waves.aorta_finger;
clear data


%%	Save data for subject X
Subject_ID = 3;
Subject = Data_waves(Subject_ID);


%%	Plot pressure and flow at central and peripheral locations
addpath('~/Haemodynamic_Tools/Version6/Others/PlotSave/')

%	Pressure
figure(1)
j = 1;	%Aortic root, inlet
time = Subject.onset_time(j) + linspace(0,0.002*(length(Subject.P{j})-1),length(Subject.P{j}));
plot(time,Subject.P{j},'k','LineWidth',4)
% xlim([0,0.8])
ylim([70,130])
set(gca,'XTick',0:0.2:0.8)
set(gca,'YTick',[80:20:120])
set(gca,'FontSize',40)
PlotSave(PathRoot,['Central_BP_',num2str(Subject_ID)])

figure(2)
j = 20;	%Brachial artery, middle
time = Subject.onset_time(j) + linspace(0,0.002*(length(Subject.P{j})-1),length(Subject.P{j}));
plot(time,Subject.P{j},'k','LineWidth',4)
xlim([0,0.8])
ylim([70,130])
set(gca,'XTick',0:0.2:0.8)
set(gca,'YTick',[80:20:120])
set(gca,'FontSize',40)
PlotSave(PathRoot,['Peripheral_BP_',num2str(Subject_ID)])

%	Flow
figure(3)
j = 1;	%Aortic root, inlet
time = Subject.onset_time(j) + linspace(0,0.002*(length(Subject.P{j})-1),length(Subject.P{j}));
plot(time,Subject.U{j}.*Subject.A{j}*10^6,'k','LineWidth',4)
xlim([0,0.8])
% ylim(1e-4*[-1,4])
set(gca,'XTick',0:0.2:0.8)
% set(gca,'YTick',[])
set(gca,'FontSize',40)
PlotSave(PathRoot,['Central_Q_',num2str(Subject_ID)])

% figure(4)
% j = 20;	%Brachial artery, middle
% time = Subject.onset_time(j) + linspace(0,0.002*(length(Subject.P{j})-1),length(Subject.P{j}));
% plot(time,Subject.U{j}.*Subject.A{j},'k','LineWidth',4)
% xlim([0,0.8])
% % ylim(1e-4*[-1,4])
% set(gca,'XTick',0:0.2:0.8)
% set(gca,'YTick',[])
% set(gca,'FontSize',40)
% PlotSave(PathRoot,['Peripheral_Q_',num2str(Subject_ID)])


