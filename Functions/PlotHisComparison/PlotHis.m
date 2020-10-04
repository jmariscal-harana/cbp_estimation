%%	Plot Nektar1D histograms (.his)
%	Input
%	-PATHS: folders where required data is stored
%	-PLOTS: choose the type of plots which are generated
%	-Art_loc: arterial location (domain); e.g. 1 (Ascending Aorta)
%	-Dom_loc: location of history point along a given domain; e.g. 1 (x =
%		0.0 [m]), 2 (x = 0.1 [m]), etc.
%
%	Output
%	-Plot comparing different .his files
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (23/05/18)
%
%==========================================================================

function PlotHis(PATHS,PLOTS,Art_loc,Dom_loc,Waveform)
%	Paths to subjects' .his data
SimulationFolder	= PATHS.Simulations;
VariationFolder		= PATHS.Variations;

addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotFormat')

%	Subjects' names
Subjects			= PATHS.FileNames;


%%	Read VD patient data, check array size, and plot pressure and flow
N_subjects = length(Subjects);	%number of VD patients

figure;
hold on
for jj = 1:N_subjects
	%	Read .his files for each subject
	if length(SimulationFolder) == 1
		Folder = [SimulationFolder,VariationFolder{jj},Subjects{jj},'_',num2str(Art_loc{jj}),'.his'];
	else
		Folder = [SimulationFolder{jj},VariationFolder{jj},Subjects{jj},'_',num2str(Art_loc{jj}),'.his'];
	end
	[Data]		= ReadHis(Folder);
	His_Skip	= max(Data(:,end));
	
	switch Waveform
		case 'P'
			plot(Data(Dom_loc:His_Skip:end,1),Data(Dom_loc:His_Skip:end,2)*PLOTS.PUnits,PLOTS.Markers{jj},'LineWidth',1);
		case 'Q'
			plot(Data(Dom_loc:His_Skip:end,1),Data(Dom_loc:His_Skip:end,4)*PLOTS.QUnits,PLOTS.Markers{jj},'LineWidth',1);
		case 'A'
			plot(Data(Dom_loc:His_Skip:end,1),Data(Dom_loc:His_Skip:end,5)*PLOTS.AUnits,PLOTS.Markers{jj},'LineWidth',1);
	end
	
	clear Data
end

FormatPlots(PLOTS)


end
