%%	Plot Nektar1D histograms (.his)
%	Input
%	-FilePath: folders where required data is stored
%	-FileName: subject name
%	-Art_loc: arterial location (domain); e.g. 1 (Ascending Aorta)
%	-His_loc: history point location along a given domain; e.g. 1 (x = 0.0 [m]), 2 (x = 0.1 [m]), etc.
%
%	Output
%	-Plot comparing different .his files from any folder
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (26/07/18)
%
%==========================================================================

function CompareHis(FilePath,FileName,Art_loc,Dom_loc)
close all

%%	Read VD patient data, check array size, and plot pressure and flow
N_subjects = length(FileName);	%number of VD patients

% Markers = {'k-', 'r--', 'g--', 'b--', 'c--'};
Scale_P = 1/133.322365;
Scale_Q = 10^6;

figure;
title('Pressure')
hold on
for jj = 1:N_subjects
	%	Read .his files for each subject
	File = [FilePath{jj},'/',FileName{jj},'_',num2str(Art_loc),'.his'];
	[Data] = ReadHis(File);
	His_Skip	= max(Data(:,end));
	
	%	Plot pressure
% 	plot(Data(Dom_loc:His_Skip:end,1),Scale_P*Data(Dom_loc:His_Skip:end,2),Markers{jj},'LineWidth',2);
	plot(Data(Dom_loc:His_Skip:end,1),Scale_P*Data(Dom_loc:His_Skip:end,2),'LineWidth',2);
		
	clear Data
end
hold off

% figure;
% title('Flow')
% hold on
% for jj = 1:N_subjects
% 	%	Read .his files for each subject
% 	File = [FilePath{jj},'/',FileName{jj},'_',num2str(Art_loc),'.his'];
% 	[Data] = ReadHis(File);
% 	His_Skip	= max(Data(:,end));
% 	
% 	%	Plot flow
% % 	plot(Data(Dom_loc:His_Skip:end,1),Scale_Q*Data(Dom_loc:His_Skip:end,4),Markers{jj},'LineWidth',2);
% 	plot(Data(Dom_loc:His_Skip:end,1),Scale_Q*Data(Dom_loc:His_Skip:end,4),'LineWidth',2);
% 	
% 	clear Data
% end
% hold off

% addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/FormatPlots')
% FormatPlots(PLOTS)

end