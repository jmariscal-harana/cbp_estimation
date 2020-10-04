%%	Example run of PlotHis.m
clear;
close all;


%%	PATHS
PlotID = 12;

switch PlotID
	case 1	%55-art vs Virtual Ubuntu
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_11_07/';'~/VirtualUbuntu/'};
		PATHS.Variations	= {''; ''};
		PATHS.FileNames		= {'55art_0000-000X'; 'Aortic_0000-000X'};
	case 2	%55-art vs Aortic
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/'};
		PATHS.Variations	= {'2018_11_07/';'2018_11_12/'};
		PATHS.FileNames		= {'55art_ZZZZ-000X';'Aortic_ZZZZ-000X'};	
	case 3	%55-art intra-subject
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/'};
		PATHS.Variations	= {'2018_12_12/';'2018_12_12/';'2018_12_12/'};
		PATHS.FileNames		= {'55art_AAZZ-ZZZX';'55art_AAZZ-ZZZX_long';'55art_AAZZ-ZZZX_Pext'};
	case 4
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/';'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_05_22/'};
		PATHS.Variations	= {'/'; '/'; '/'; '/'; '/'; '/'; '/'; '/'; '/'};
		PATHS.FileNames		= {'55art_DB_Baseline_0000-0000';
			'55art_DB_Baseline_0000-00A0'
			'55art_DB_Baseline_0000-00Z0'
			'55art_DB_Baseline_0000-0AA0'
			'55art_DB_Baseline_0000-0A00'
			'55art_DB_Baseline_0000-0AZ0'
			'55art_DB_Baseline_0000-0ZA0'
			'55art_DB_Baseline_0000-0Z00'
			'55art_DB_Baseline_0000-0ZZ0'};
	case 5
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_10_24/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_10_24/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_10_24/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_10_24/'};
		PATHS.Variations	= {'2018_09_12/';'Same_geo_PWV_Empirical/'; 'Tex_geo_same_PWV_Empirical/'; 'Tex_geo_PWV_Empirical/'; 'Tex_geo_f2f_PWV_Empirical/'};
		PATHS.FileNames		= {'55art_0000-000X';'Aortic_0000-000X';'Aortic_0000-000X';'Aortic_0000-000X';'Aortic_0000-000X'};
	case 6
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_10_24/';'/Users/joh15/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_10_24/'};
		PATHS.Variations	= {'2018_09_12/';'Tex_geo_f2f_PWV_Beta/'; 'Tex_geo_PWV_Beta/'};
		PATHS.FileNames		= {'55art_0000-000X';'Aortic_0000-000X';'Aortic_0000-000X'};		
	case 7	
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_11_07/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_11_15/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_11_15/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_11_15/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2018_11_15/'};
		PATHS.Variations	= {'','','R1_80/','R1_100/','R1_120/'};
		PATHS.FileNames		= {'55art_0000-000X';'Aortic_0000-000X';'Aortic_0000-000X_R1_80';'Aortic_0000-000X_R1_100';'Aortic_0000-000X_R1_120'};
	case 8
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_11_21/', '~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_11_07/'};
		PATHS.Variations	= {'',''};
		PATHS.FileNames		= {'M55art','55art_0000-000X'};
	case 9	%	dP study - David Nordsletten
		Cohort = 'Patients_normal';
		Path = ['~/PhD/Collaborations/David Nordsletten/Nektar_outputFiles/Simulations/2019_02_12/',Cohort,'/'];
		Cohort_size = 11;
		PATHS.Simulations	= repmat({Path},1,Cohort_size);
		PATHS.Variations	= repmat({''},1,Cohort_size);
		PATHS.FileNames		= importdata([Path,'Nektar_Success.txt'])';
		PATHS.Figures		= Path;
	case 10	%116-art vs 1-D Aortic
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';'~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/'};
		PATHS.Variations	= {'2019_05_09/';'2019_05_09/Theoretical_6_subjects/';'2019_05_09/';'2019_05_09/Theoretical_6_subjects/';'2019_05_09/';'2019_05_09/Theoretical_6_subjects/';'2019_05_09/';'2019_05_09/Theoretical_6_subjects/';'2019_05_09/';'2019_05_09/Theoretical_6_subjects/';'2019_05_09/';'2019_05_09/Theoretical_6_subjects/';};
		PATHS.FileNames		= {'sim_1';'Aortic1D_pwdb_1';'sim_2';'Aortic1D_pwdb_2';'sim_3';'Aortic1D_pwdb_3';'sim_4';'Aortic1D_pwdb_4';'sim_5';'Aortic1D_pwdb_5';'sim_6';'Aortic1D_pwdb_6'};
		PATHS.Figures		= ['~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Figures/',PATHS.Variations{1},'/'];
	case 11 %116-art baseline vs 1-D Aortic Scale 1
		PATHS.Simulations	= {'~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/', '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/'};
		PATHS.Variations	= {'2019_05_09/';'2019_05_10/';'2019_05_10/';'2019_05_10/';'2019_05_10/';'2019_05_10/'};
		PATHS.FileNames		= {'sim_1';'Aortic1D_pwdb_1_R1_A';'Aortic1D_pwdb_1_R1_B';'Aortic1D_pwdb_1_R1_0';'Aortic1D_pwdb_1_R1_Y';'Aortic1D_pwdb_1_R1_Z';};
		PATHS.Figures		= ['~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Figures/',PATHS.Variations{1}];
	case 12	%Virtual Ubuntu comparison
		PATHS.Simulations	= {'~/VirtualUbuntu/';'~/VirtualUbuntu/'};
		PATHS.Variations	= {''; ''};
		PATHS.FileNames		= {'Aorta_1_absorb'; 'Aorta_1_reflect'};

end


%%	PLOTS
PLOTS.WhichPlots			= 0;	% 0: Plot Wk reference
									% 1: Plot Wk estimation
									% 2: Plot reference vs estimation
PLOTS.OptimisedBCs			= 0;	% 1: Include 'optimised' simulations
PLOTS.FullSimulation		= 0;	% 0: Plot last simulation cycle (quasi-steady state)
									% 1: Plot simulations from the start of the first cycle

%	Title and other text
PLOTS.Title					= 1;	% 1: Display title on each figure
PLOTS.Legend				= 0;	% 1: Display legends on each figure
PLOTS.XLabel				= 1;	% 1: Display x-axis labels
PLOTS.YLabel				= 1;	% 1: Display y-axis labels
PLOTS.Annotation			= 0;	% 1: Display annotations

PLOTS.Title_text			= {'Theoretical study'};
PLOTS.Legend_text			= PATHS.FileNames;
PLOTS.YLabel_text			= {'Pressure [mmHg]'};
PLOTS.XLabel_text			= {'Time [s]'};

%	Axis properties	
PLOTS.XLim					= 2;	% 0: no limits, no axis
									% 1: 'custom'
									% 2: 'auto'
PLOTS.YLim					= 2;	% 0: no limits, no axis
									% 1: 'custom'
									% 2: 'auto'
PLOTS.Grid					= 1;	% 1: grid on
									
%	Error box
PLOTS.ErrorBox				= 0;	% 1: Display absolute error box
									% 2: Display relative error box
%	Others
PLOTS.PUnits				= 1/133.322368;	% From Pa to mmHg
PLOTS.QUnits				= 10^6;			% From m3/s to mL/s
PLOTS.AUnits				= 10^4;			% From m^2 to cm^2



% PLOTS.Markers				= {'k','--b','--r','--k','--g','--c'};
PLOTS.Markers				= {'k','--b','k','--b','k','--b','k','--b','k','--b','k','--b',};

%	Fontsize
PLOTS.FontSize				= 30;
PLOTS.LegendSize			= 10;


%%	Call function PlotHis.m
Waveform = 'A';
ArterialLoc = {1 1 1 1 1 1 1 1 1 1 1 1};			% 1: Aortic root
HistogramLoc = 1;			% 1: x = 0 [m]

addpath('~/Haemodynamic_Tools/Version6/Others/PlotHisComparison/')
PlotHis(PATHS,PLOTS,ArterialLoc,HistogramLoc,Waveform)
legend(PLOTS.Legend_text,'location','southeast','FontSize',PLOTS.LegendSize);

addpath('~/Haemodynamic_Tools/Version6/Others/PlotSave/')
PlotSave(PATHS.Figures,'CYCLES_all')


%%	Important - avoid using functions from different folder with same names
restoredefaultpath

