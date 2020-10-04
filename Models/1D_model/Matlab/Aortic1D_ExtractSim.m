%%	Extract Nektar1D haemodynamic data
%	Input
%	-Nektar1D .his files
%
%	Output
%	-EST: formatted clinical data
%	?	Patient info: ID, Age, Gender, etc.
%	?	LOCATION. ? arterial measurement location: Ao_Root, Asc_Ao
%		o	CYCLES. ? local haemodynamic data for multiple cardiac cycles
%		o	ONE_CYCLE. ? local haemodynamic data for a single cycle
%			?	t., P., Q., etc.
%			?	SBP, MBP, DBP, etc.
%	?	Global haemodynamic parameters: SV, HR, CO, PWV, etc.
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (28/02/19)
%	v2.0 (08/06/19) - updated for compatibility with ConvertHistoryFiles,
%		applied physiological filter
%
%==========================================================================
clear;
close all;


%%	PATHS
PATHS.Root		= '~/Haemodynamic_Tools/Version6/';
PATHS.REF_data	= [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Simulations/'];
PATHS.EST_data	= [PATHS.Root,'Aortic1D/Nektar_outputFiles/Simulations/2019_09_19/Sc2/5226_0_3/'];

REF_Population	= 'Isra_AoCoarctation';


%%	Scaling factors for pressure, flow and area
% Scale_P = 1/133.322;
% Scale_A = 1e6;
% Scale_Q = 1e6;
Scale_P = 1;
Scale_A = 1;
Scale_Q = 1;


%%	Select .HIS to extract
Artery_name					= {'Asc_Ao'}';
Artery_segment				= [1];
Artery_element				= [{1}];

%	'up' parameters for ConvertHistoryFiles.m
up.all_beats				= false;
up.all_data					= false;
up.required_signals			= {'P', 'Q', 'A'};
up.dir						= PATHS.EST_data;
up.required_domains			= Artery_segment;	%[Ao Root, Brachio, L Common Carotid, Desc Ao, L Subclavian, L Brachial, R Femoral]
up.required_distance_els	= Artery_element;
up.ds_factor				= 2;				%Downsample factor, integer, > 1

fid		= fopen([PATHS.EST_data,'Nektar_Success.txt']);
ID		= textscan(fid,'%s');
ID		= ID{1};

for jj = 1:length(ID)
	EST(jj).INFO.ID		= ID{jj};

	for kk = 1:length(Artery_name)
		eval(['EST(jj).', Artery_name{kk}, '.Segment = ', num2str(Artery_segment(kk)), ';'])
	end
end

addpath([PATHS.Root,'/Others/ConvertHistoryFiles/'])
data = ExtractHisData(up,EST);


%%	Save data into EST structure
%	Arterial waveforms of interest
Artery_wave_name		= {'Asc_Ao'};

ID = fieldnames(data);

for jj = 1:length(EST)
	%%	Arterial waveforms
	for kk = 1:length(Artery_wave_name)
		eval(['fs												= data.',EST(jj).INFO.ID,'.fs;']);
		dt														= 1/fs;
		eval(['EST(jj).',Artery_wave_name{kk},'.ONE_CYCLE.P		= data.',EST(jj).INFO.ID,'.P*Scale_P;']);
		eval(['EST(jj).',Artery_wave_name{kk},'.ONE_CYCLE.A		= data.',EST(jj).INFO.ID,'.A*Scale_A;']);
		eval(['EST(jj).',Artery_wave_name{kk},'.ONE_CYCLE.Q		= data.',EST(jj).INFO.ID,'.Q*Scale_Q;']);
		eval(['EST(jj).',Artery_wave_name{kk},'.ONE_CYCLE.t		= [0:1/fs:(length(EST(jj).',Artery_wave_name{kk},'.ONE_CYCLE.P)-1)/fs]'';'])					%[s]
	end
	
end


%%	Physiological filter based on the Normotensive and Hypertensive datasets (IN SILICO ONLY)
% %	Load REF dataset
% load([PATHS.REF_data,REF_Population,'_reference_v2.mat'],'REF')
% 
% %	Store waveform pressure data into REF_temp to save space
% for jj = 1:length(REF)
% 	REF_temp(jj).Pressure		= REF(jj).Asc_Ao.ONE_CYCLE.P;
% 	REF_temp(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;
% 	REF_temp(jj).ID				= REF(jj).INFO.ID;
% end
% 
% addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
% REF_temp = Physiological_BP_Filter(PATHS,REF_temp);
% clear REF
% 
% %	Compare physiological (filtered) cBP waveforms only
% k = 1;
% for jj = 1:length(EST)
% 	if strcmp([REF_temp(k).ID,'_Aortic1D'],EST(jj).INFO.ID)
% 		EST_temp(k) = EST(jj);
% 		k = k+1;
% 	end
% end
% EST = EST_temp; clear EST_temp


%%	Choose save folder
% save([PATHS.EST_data,'Aortic1D_pwdb_4374_estimation'],'EST')
% save([PATHS.EST_data,'Aortic1D_pwdb_6_estimation'],'EST')
% save([PATHS.EST_data,'Aortic1D_pwdb_78_estimation'],'EST')
save([PATHS.EST_data,'Aortic1D_AoCo_estimation'],'EST')


%%
restoredefaultpath


%%	Function definitions
function [data] = ExtractHisData(up,EST)
if exist([up.dir,'history_files_data.mat'],'file') ~= 2
	counter = 1;
	for jj = 1:length(EST)
		for kk = up.required_domains
			up.filename{counter} = [EST(jj).INFO.ID, '_', num2str(kk), '.his'];
			counter = counter + 1;
		end
	end
	
	data = ConvertHistoryFiles(up);
else
	load([up.dir,'history_files_data.mat'],'data')
end
end
