%%	Main code for the CV parameter estimation study: individual analysis
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
%	v1.0 (22/05/19) - Based on Wk_Parameter_Study
%
%==========================================================================
format compact;
close all;
clear;
clc;
% dbstop if error


%%	REF input parameters
REF_Dataset = 'pwdb_4374';	% Reference dataset (virtual): 'Wk2', 'Wk3', '55_art', 'pwdb_78', 'pwdb_4374'
disp(['Reference dataset: ',REF_Dataset])


%%	EST input parameters
EST_Dataset		= 'pwdb_4374_carotid';
Scenario		= 'Waveform';	% 'Waveform': pressure waveform available
								% 'BP': only DBP/MBP/SBP values available
Plots_EST		= 0;			% 1: show estimation method plots

%	Estimation methods
MethodLVET		= 0;		% 0: LVET = LVET,reference
								% 1: LVET from P using Sam's algorithm
								% 2: LVET from P using Pete's PPG algorithm
								% 3: LVET from P using WF								
								% 4: LVET = 0.37*sqrt(T)
								% 5: Q=0 or max(Q) after min Q
MethodP_out		= 0;%[1:5];		% 0: P_out = P_out,reference
								% 1: P_out = 0
								% 2: P_out = P_dias/2
								% 3: P_out = 0.7 * DBP (Parragh 2014)
								% 4: Decay time (Jorge): t = LVET
								% 5: Decay time (Jorge): t = LVET + 1/3 t_dias
MethodR_T		= 0;			% 0: R_T = R_T,reference
								% 1: Scenario 1
								% 2: Scenario 2
MethodC_T		= 0;%[1:9];			% 0: C_T = C_T,reference
								% 1: 2-point exponential fit
								% 2: Decay time (Jorge): t = LVET
								% 3: Decay time (Jorge): t = LVET + 1/3 t_dias
								% 4: Area method
								% 5: Two-area method
								% 6: Diastolic pressure method
								% 7: Pulse pressure method
								% 8: SV/PP
								% 9: 3-Wk optimisation
MethodPWV		= 1;%[1:5];			% 0: PWV = PWV,reference
								% 1: TT: ascending-descending aorta Q
								% 2: TT: carotid-femoral P
								% 3: Least-squares: ascending-descending aorta Q
								% 4: Least-squares: carotid-femoral P
								% 5: Sum of squares
MethodZ_0		= 0;			% 0: Z_0 = Z_0,reference
								% 1: Z_0 = 0.1*R_T (Murgo1980)
								% 2: Z_0 = 0.05*R_T (Murgo1980)
								% 3: Z_0 = (P_mean - P_dias)/Q_peak (Jordi)
								% 4: rho*PWV/A_d
								% 5: Freq-domain: 2-12th harmonics (Nichols1977)
								% 6: Freq-domain: 6-10th harmonics (Segers2000)
								% 7: Freq-domain: 1-8th harmonics (Clarke et al (1978))
								% 8: Freq-domain: 1-9th harmonics (Peluso et al (1978))
								% 9: Freq-domain: 2-10th harmonics (Mitchell et al (1994))
								% 10: Freq-domain: 3-10th harmonics (Hughes and Parker (2009))
								% 11: Freq-domain: 4-10th harmonics (Tabima et al (2012))
								% 12: Freq-domain: 6-8th harmonics (Abel (1971))
								% 13: Freq-domain: 4-8th harmonics (Qureshi et al (2018))
								% 14: Early systole up-slope?(Dujardin and Stone (1981))
								% 15: Least-squares fit for early-systole
								% 16: P_es to Q_es ratio at time of maximum dQ
								% 17: P_es to Q_es ratio [t = 0, t = maximum dQ]
								% 18: 3-Wk optimisation

switch Scenario
	case 'Waveform'
		Sc = 'Sc1';
	case 'BP'
		Sc = 'Sc2';
end

								
%% PATHS
PATHS.Root					= '~/Haemodynamic_Tools/Version6/';	% Path to the Windkessel study folder
PATHS.ParameterEstimation	= [PATHS.Root,'ParameterEstimation/'];
PATHS.Matlab				= [PATHS.ParameterEstimation,'Matlab/'];
PATHS.Output				= [PATHS.ParameterEstimation,'Output/',EST_Dataset,'/'];
PATHS.Figures				= [PATHS.ParameterEstimation,'Figures/'];

addpath(PATHS.Matlab)


%%	REF data: choose virtual dataset
switch REF_Dataset
	case {'Wk2', 'Wk3'}
		PATHS.Input = [PATHS.Root,'Windkessel/',REF_Dataset,'_input/'];
		load([PATHS.Input,REF_Dataset,'_reference_v2.mat'],'REF')
		
		addpath([PATHS.Root,'Others/ExtractReferenceBP/']);
		
		%	Format 0-D dataset for CV parameter estimation
		for jj = 1:length(REF)
			REF(jj).Dataset			= REF_Dataset;
			REF(jj).Pressure		= REF(jj).Pressure*133.322368;		%[Pa]
			REF(jj).Pressure_CBP	= REF(jj).Pressure;					%[Pa]
			REF(jj).Q_in			= REF(jj).Q_in/10^6;				%[m3/s]
			REF(jj).DBP				= min(REF(jj).Pressure);			% DBP
			REF(jj).SBP				= max(REF(jj).Pressure);			% SBP
			REF(jj).PWV				= NaN;
			REF(jj).Z_0				= REF(jj).Z_0*(133.322368*10^6);
			REF(jj).P_out			= REF(jj).P_out*133.322368;
			REF(jj).R_T				= REF(jj).R_T*(133.322368*10^6);
			REF(jj).C_T				= REF(jj).C_T/(133.322368*10^6);
			REF(jj).SV				= REF(jj).SV/10^6;					%[m3]

			
			REF_temp(jj)			= ExtractReferenceBP(Scenario,REF(jj));
		end
		
		REF = REF_temp; clear REF_temp
				
		%	Physiological filter based on the Normotensive and Hypertensive datasets
		addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
		REF = Physiological_BP_Filter(PATHS,REF);
		
	case '55_art'
		PATHS.Input = [PATHS.Root,'VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_12_12/'];
		load([PATHS.Input,REF_Dataset,'_reference.mat'],'REF')
		
		addpath([PATHS.Root,'Others/ExtractReferenceBP/']);
		
		%	Format 55-artery dataset for CV parameter estimation
		for jj = 1:length(REF)
			REF(jj).Dataset			= REF_Dataset;
			REF(jj).Pressure		= REF(jj).L_Brachial.ONE_CYCLE.P;
			REF(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;
			REF(jj).t_in			= REF(jj).Asc_Ao.ONE_CYCLE.t;
			REF(jj).Q_in			= REF(jj).Asc_Ao.ONE_CYCLE.Q;
			REF(jj).SV				= REF(jj).HAEMO.SV;
			REF(jj).DBP				= min(REF(jj).Pressure);	% DBP
			REF(jj).SBP				= max(REF(jj).Pressure);	% SBP
			
			REF_temp(jj)			= ExtractReferenceBP(Scenario,REF(jj));
		end
		
		REF = REF_temp; clear REF_temp
		
	case {'pwdb_78', 'pwdb_4374'}
		PATHS.Input = [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Simulations/'];
		load([PATHS.Input,REF_Dataset,'_reference_v2.mat'],'REF')
		warning('Peripheral BP for Scenario 1 from carotid wave (previously brachial)')
		
		addpath([PATHS.Root,'Others/ExtractReferenceBP/'])
		
		%	Format Pete's dataset for CV parameter estimation
		for jj = 1:length(REF)
			REF(jj).INFO.ID			= REF(jj).INFO.ID;
			REF(jj).Dataset			= REF_Dataset;
			switch Scenario
				case 'Waveform'
					REF(jj).Pressure		= REF(jj).L_Carotid.ONE_CYCLE.P;	%[Pa]
				case 'BP'
					REF(jj).Pressure		= REF(jj).L_Brachial.ONE_CYCLE.P;	%[Pa]
			end
			REF(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;		%[Pa]
			REF(jj).t_in			= REF(jj).Asc_Ao.ONE_CYCLE.t;		%[s]
			REF(jj).Q_in			= REF(jj).Asc_Ao.ONE_CYCLE.Q;		%[m3/s]
			REF(jj).SV				= REF(jj).HAEMO.SV;					%[m3]
			REF(jj).LVET			= REF(jj).HAEMO.LVET;				%[s]
			REF(jj).PWV				= REF(jj).PWV;					%[m/s]
			REF(jj).P_out			= REF(jj).P_out;					
			REF(jj).R_T				= REF(jj).R_T;					
			REF(jj).C_T				= REF(jj).C_T;					
			REF(jj).Z_0				= REF(jj).Z_0;					
			REF(jj).DBP				= min(REF(jj).Pressure);	% DBP
			REF(jj).SBP				= max(REF(jj).Pressure);	% SBP
			
			REF_temp(jj)			= ExtractReferenceBP(Scenario,REF(jj));
		end
			
		REF = REF_temp; clear REF_temp
		
% 		for jj = 1:length(REF)
% 			MBP_c(jj) = mean(REF(jj).Pressure_CBP);
% 			MBP_p(jj) = mean(REF(jj).Pressure);
% 		end
% 		figure, hold on, plot(MBP_p - MBP_c,'o')
		
		%	Physiological filter based on the Normotensive and Hypertensive datasets
		addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
		REF = Physiological_BP_Filter(PATHS,REF);
		
	otherwise
		error('Load a valid dataset')
end

%TEMP
load('/Users/joh15/Haemodynamic_Tools/Version6/ParameterEstimation/Output/pwdb_4374_carotid/PWV_outliers.mat')
% REF = REF(PWV_outliers);


%%	EST data: estimate CV parameters
addpath([PATHS.Root,'Others/LVET/'])
addpath([PATHS.Root,'Others/PulseAnalyse/'])
addpath([PATHS.Root,'Others/P_out/'])
addpath([PATHS.Root,'Others/PWV_TT/'])
addpath([PATHS.Root,'Others/ImpedanceAnalysis/'])
addpath([PATHS.Root,'Others/export_fig'])
addpath([PATHS.Root,'Others/P_0_iteration/'])

EST = struct;
k = 1;


%	Loop through each estimation method
for a1 = MethodLVET
	for a2 = MethodP_out
		for a3 = MethodR_T
			for a4 = MethodC_T
				for a5 = MethodPWV
					for a6 = MethodZ_0
						%	Loop through every patient
						for jj = 1:length(REF)
							%	1. Estimation parameters
							EST(jj).ID				= REF(jj).INFO.ID;
							EST(jj).ID_num			= jj;
							
							%	2. Estimation methods
							EST(jj).MethodLVET		= a1;
							EST(jj).MethodP_out		= a2;
							EST(jj).MethodR_T		= a3;
							EST(jj).MethodC_T		= a4;
							EST(jj).MethodPWV		= a5;
							EST(jj).MethodZ_0		= a6;
						end
						
						Methods = [num2str(a1),num2str(a2),num2str(a3),num2str(a4),'_',num2str(a5),'_',num2str(a6)];
						fprintf('Current estimation method: %s\n',Methods)
						
						%	3.a Estimate parameters
						for jj = 1:length(REF)
							EST_temp(jj) = ParameterEstimation(PATHS,REF(jj),EST(jj),Plots_EST);
						end
						
						EST = EST_temp; clear EST_temp
						
% 							Save parameter data
						save([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')

						%	3.b Load estimated parameters
% 						load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
% 						%	Compare the right studies (post-physiological filter)					
% 						for jj = 1:length(REF)
% 							for kk = 1:length(EST)
% 								if contains(REF(jj).INFO.ID,EST(kk).ID)
% 									EST_index(jj) = kk;
% 									break
% 								end
% 							end
% 						end
% 						for jj = 1:length(EST_index)
% 							temp(jj) = EST(EST_index(jj));
% 						end
% 						
% 						EST = temp; clear temp
% 												
% 						%	Update parameter data
% 						save([PATHS.Output,'EST_',Dataset,'_',Methods],'EST')
						
						%	Calculate individual (intrasubject) errors
						ERR(k).Methods = Methods;
						
						if a1 ~= 0, Parameter = 'LVET'; end
						if a2 ~= 0, Parameter = 'P_out'; end
						if a3 ~= 0, Parameter = 'R_T';	end
						if a4 ~= 0, Parameter = 'C_T';	end
						if a5 ~= 0,	Parameter = 'PWV'; end
						if a6 ~= 0,	Parameter = 'Z_0'; end
						
						ERROR_POP(k) = ParameterError(REF,EST,Parameter);
						
						eval(['fprintf(''Mean +/- SD relative error: %.1f +/- %.1f\n'',','ERROR_POP(k).Rel.',Parameter,'.Mean,ERROR_POP(k).Rel.',Parameter,'.SD)'])
						
						k = k+1;
							
					end
				end
			end
		end
	end
end

figure, hold on, plot([1:length(EST)],[EST.PWV],'ko')
plot(PWV_outliers,[EST(PWV_outliers).PWV],'rx')

clear REF EST


%%	Save overall results
Parameters = {'LVET', 'P_out', 'R_T', 'C_T', 'PWV', 'Z_0'};
Parameters = Parameters(1,[a1 ~= 0, a2 ~= 0, a3 ~= 0, a4 ~= 0, a5 ~= 0, a6 ~= 0]);

if sum([a1 a2 a3 a4 a5 a6]) ~= 0
	save([PATHS.Output,'PARAM_',EST_Dataset,'_',Parameters{:},'_',Sc],'ERROR_POP')
end


%%	Print to text
if sum([a1 a2 a3 a4 a5 a6]) ~= 0
	for kk = 1:length(Parameters)
		fileID = fopen([PATHS.Output,'PARAM_',EST_Dataset,'_',Parameters{kk},'_',Sc,'.txt'],'w');
		
		%	Absolute errors
		fprintf(fileID,'Mean +/- SD parameter estimation absolute errors [mmHg]\n\n');
		
		for jj = 1:length(ERR)
			fprintf(fileID,'Estimation method: %s\n',ERR(jj).Methods);
			
			fprintf(fileID,'& %3.1f $\\pm$ %3.1f\n\n',...
				[eval(['ERROR_POP(jj).Abs.',Parameters{kk},'.Mean']);
				eval(['ERROR_POP(jj).Abs.',Parameters{kk},'.SD'])]');
		end
		
		fprintf(fileID,'---------------------------------------\n\n');
		
		%	Relative errors
		fprintf(fileID,'Mean +/- SD parameter estimation relative errors [%%]\n\n');
		
		for jj = 1:length(ERR)
			fprintf(fileID,'Estimation method: %s\n',ERR(jj).Methods);
			
			fprintf(fileID,'& %3.1f $\\pm$ %3.1f\n\n',...
				[eval(['ERROR_POP(jj).Rel.',Parameters{kk},'.Mean']);
				eval(['ERROR_POP(jj).Rel.',Parameters{kk},'.SD'])]');
		end
		fclose(fileID);
	end
end


%%
restoredefaultpath

