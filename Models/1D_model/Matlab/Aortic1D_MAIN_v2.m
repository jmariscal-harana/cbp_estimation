%%	Create Nektar1D input files (.IN) for the Aortic1D model
%	1. Extract REF. data for each subject or patient
%	2. Estimate cardiovascular parameters from REF. data
%	3. Create Aortic1D models from REF. data
%
%	Input
%	-PATHS: folders where I/O data is stored
%	-REF: reference parameters required to create an Aortic1D model
%
%	Output
%	-.IN and .BCS containing patient-specific files for Aortic1D simulations
%	-.MAT containing values required to generate .IN and .BCS files
%
%	For more information, open most recent READ_ME.DOCX file
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (23/01/19) - based on 'Aortic1D_55art.m'
%	v2.0 (14/03/19) - code cleanup and rearranging
%
%==========================================================================
format compact;
close all;
clear;
clc;
dbstop if error


%%	PATHS
PATHS.Root				= '~/Haemodynamic_Tools/Version6/';					%Root folder
PATHS.Aortic1D			= [PATHS.Root,'Aortic1D/'];							%Aortic1D folder
PATHS.Matlab			= [PATHS.Aortic1D,'Matlab/'];						%Matlab functions for Aortic1D
PATHS.Aortic1D_files	= [PATHS.Aortic1D,'Nektar_inputFiles/'];			%.BCS, .IN and .MAT aortic model files are saved here

addpath(PATHS.Matlab)
addpath([PATHS.Root,'Others/PlotShaded/']);


%%	1. Extract haemodynamic data as REF. with a standard format for every dataset
addpath([PATHS.Root,'Others/ExtractReferenceBP/']);

%	Dataset type
FormatREFData	= 1;	%0: format REF data from in silico or in vivo dataset;
						%1: load Aortic1D REF. data
						%2: load REF data for theoretical study
						
Scenario		= 'Waveform';	% 'Waveform': pressure waveform available
Plots			= 0;			% 1: show estimation method plots

% Dataset			= 'pwdb_4374_reference_v2';
% Dataset			= 'pwdb_6_reference';	%Theoretical study with 116-artery model baseline subjects
Dataset			= 'Isra_AoCoarctation_reference';

switch Dataset
	case 'Isra_AoCoarctation_reference'
		PATHS.REF_data = '~/Clinical_Data/Isra_AoCoarctation/';
	case {'pwdb_4374_reference_v2','pwdb_4374_reference','pwdb_6_reference'}
		PATHS.REF_data = '~/Haemodynamic_Tools/Version6/VirtualDB/pwdb/Nektar_outputFiles/Simulations/';
end

switch Scenario
	case 'Waveform'
		Sc = 'Sc1';
	case 'BP'
		Sc = 'Sc2';
end

if FormatREFData == 0	% In vivo model (same info as in vivo)
	load([PATHS.REF_data,Dataset,'.mat'],'REF')
	for jj = 1:length(REF)
		REF_temp(jj) = Aortic1D_ExtractREF_v2(PATHS,REF(jj),Dataset,Scenario);
	end
	REF = REF_temp;
	save([PATHS.REF_data,Dataset,'_Aortic1D_',Sc,'.mat'],'REF')
	
elseif FormatREFData == 1
	load([PATHS.REF_data,Dataset,'_Aortic1D_',Sc,'.mat'],'REF')
	
elseif FormatREFData == 2	% Theoretical model (some information not available in vivo)
	load([PATHS.REF_data,Dataset,'.mat'],'REF')
	for jj = 1:length(REF)
		REF_temp(jj) = Aortic1D_ExtractREF_theory(REF(jj),Dataset);
	end
	REF = REF_temp;
	save([PATHS.REF_data,Dataset,'_Aortic1D_theory.mat'],'REF')

end

if exist('REF_temp') == 1
	clear REF_temp
end


%%	2. Estimate CV parameters from REF. data
addpath([PATHS.Root,'ParameterEstimation/Matlab/']);
addpath([PATHS.Root,'Others/P_0_iteration/']);

CV_estimation_required	= 1;

switch CV_estimation_required
	case 0

	case 1				
		%	Estimation methods
		MethodLVET		= 5;
		MethodP_out		= 4;
		MethodR_T		= 1;
		MethodC_T		= 9;
		MethodPWV		= 0;
		MethodZ_0		= 6;
		
% 		global P_out_negative_counter
% 		P_out_negative_counter = 0;
		
		for jj = 1:length(REF)
			%	Estimation parameters			
			EST(jj).MethodLVET		= MethodLVET;
			EST(jj).MethodP_out		= MethodP_out;
			EST(jj).MethodR_T		= MethodR_T;
			EST(jj).MethodC_T		= MethodC_T;
			EST(jj).MethodPWV		= MethodPWV;
			EST(jj).MethodZ_0		= MethodZ_0;
		end
		
		for jj = 1:length(REF)
			%	Estimate CV parameters
			EST_temp = ParameterEstimation(PATHS,REF(jj),EST(jj),Plots);
			
			REF(jj).P_out = EST_temp.P_out;
			REF(jj).LVET = EST_temp.LVET;
			REF(jj).R_T = EST_temp.R_T;
			REF(jj).C_T = EST_temp.C_T;
% 			REF(jj).PWV = EST_temp.PWV;
			REF(jj).tau = EST_temp.tau;	
		end
		
% 		disp(P_out_negative_counter)
		
		save([PATHS.REF_data,Dataset,'_Aortic1D_',Sc,'.mat'],'REF')
		
end


%%	3. Create an Aortic1D model for each subject/patient
fileCommand = 'command_cloud.txt';
fidCommand  = fopen([PATHS.Aortic1D_files, fileCommand],'w');

% figure, hold on, for jj = 1:length(REF), plot(REF(jj).Rad_prox,'o'), end, legend(EST.ID)
% figure, hold on, for jj = 1:length(REF), plot(REF(jj).Rad_dist,'o'), end, legend(EST.ID)
% figure, hold on, for jj = 1:length(REF), plot(REF(jj).L,'o'), end, legend(EST.ID)

for jj = 1:length(REF)
	PATHS.PatientName = [REF(jj).INFO.ID,'_Aortic1D'];
	
	Aortic1D_Input_File(PATHS,REF(jj))
	
	FileName = PATHS.PatientName;
	commandLine = ['oneDbio ',FileName,'.in && echo ''', FileName,''' >> Nektar_Success.txt || echo ''', FileName,''' >> Nektar_Error.txt'];
	fprintf(fidCommand, [commandLine,'\n']);                          
end

%R_1 scaling study
% REF = REF(1);
% 
% %	Peripheral pressure measurement
% Scale_R1			= [0.5, 0.75, 1, 1.25, 1.5];
% Scale_R1_ID		= {'A', 'B', '0', 'Y', 'Z'};
% % REF(jj).Pressure = REF(jj).L_Carotid.ONE_CYCLE.P;
% 
% for jj = 1:length(REF)
% 	for kk = 1:length(Scale_R1)
% 		REF(jj).Scale_R1 = Scale_R1(kk);		
% 		PATHS.PatientName = ['Aortic1D_',REF(jj).INFO.ID,'_R1_',Scale_R1_ID{kk}];
% 		
% 		Aortic1D_Input_File(PATHS,REF(jj))
% 		
% 		FileName = PATHS.PatientName;
% 		commandLine = ['oneDbio ',FileName,'.in && echo ''', FileName,''' >> Nektar_Success.txt || echo ''', FileName,''' >> Nektar_Error.txt'];
% 		fprintf(fidCommand, [commandLine,'\n']);
% 	end
% end

fclose(fidCommand);
disp(['All input files generated in ',PATHS.Aortic1D_files])


%%	Generate script for parallel (multi-core) simulations
filepath = PATHS.Aortic1D_files;
filename = 'command_cloud';
Parallel_Simulations(filepath,filename)


restoredefaultpath
%%%%%%%%%%%%%%%%%%


%%	Function definitions
function Parallel_Simulations(filepath,filename)
%%	Number of cores
TotalCores = 8;

%% Initialize variables.
delimiter = {''};

%% Format for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

%% Open the text file.
fileID = fopen([filepath,filename,'.txt'],'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
commandcloud = [dataArray{1:end-1}];

TotalLength = length(commandcloud);
IndividualLength = ceil(TotalLength/TotalCores); 

for CurrentCore = 1:TotalCores
	fidCommand  = fopen([filepath,filename,num2str(CurrentCore),'.txt'],'w');
		for CurrentLine = 1:IndividualLength
			if ((CurrentCore-1)*IndividualLength + CurrentLine) > TotalLength
				break
			else
				fprintf(fidCommand, [commandcloud{(CurrentCore-1)*IndividualLength + CurrentLine},'\n']);
		end
	end

end
end
