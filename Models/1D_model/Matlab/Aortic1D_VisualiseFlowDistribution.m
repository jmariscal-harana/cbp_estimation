%%	Visualise Aortic1D flow distribution for the 55-artery dataset
format compact;
close all;
clear;
clc;


%%	PATHS
PATHS.REF_data			= '~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_11_07/';
PATHS.EST_data			= '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_inputFiles/2018_11_21/';
PATHS.Matlab			= '~/Haemodynamic_Tools/Version6/Aortic1D/Matlab/';						%Matlab functions

addpath(PATHS.Matlab)

Population = '55_art';
load([PATHS.REF_data,Population,'_reference.mat'],'REF')

%	EST
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/HisFiles/')

figure
hold on
for jj = 1:length([REF.INFO])
	EST(jj).INFO.ID = {['Aortic_',REF(jj).INFO.ID{1}(7:end)]};
	load([PATHS.EST_data,EST(jj).INFO.ID{1},'.mat'],'BCS');
% 	outflow{jj} = BCS.Outflow_distribution;
	outflow{jj} = BCS.Area_distribution;
	plot(4,outflow{jj}(4),'ok','MarkerSize',20)
	plot(5,outflow{jj}(5),'ok','MarkerSize',20)
	plot(6,outflow{jj}(6),'ok','MarkerSize',20)
	plot(7,outflow{jj}(7),'ok','MarkerSize',20)	
end

title('Flow distribution across 55-artery dataset')
ylabel('Outflow fraction')
xlabel('4=DescAo, 5=Brachio, 6=Carotid, 7=Subclav')
set(gca,'FontSize',20)

clear