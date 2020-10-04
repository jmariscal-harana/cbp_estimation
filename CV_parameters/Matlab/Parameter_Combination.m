%%	Main code for the CV parameter estimation study: combined analysis
%	Input
%	-REF: parameters used to generate virtual patients
%	-EST: methods for parameter estimation for virtual patients
%	-PATHS: folders where required data is stored
%	-PLOTS: choose the type of plots which are generated
%
%	Output
%	-REF: virtual patient data
%	-EST: estimated pressure using the 2-Wk or 3-Wk model
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (22/05/19) - Based on Wk_CBP_Study
%
%==========================================================================
format compact;
close all;
clear;
clc;
dbstop if error


%%	REF input parameters
Dataset				= 'pwdb_4374';	% Reference dataset (virtual): 'Wk2', 'Wk3', '55_art', 'pwdb_78', 'pwdb_4374'


%%	EST input parameters
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
MethodC_T		= 9;%[1:9];			% 0: C_T = C_T,reference
								% 1: 2-point exponential fit
								% 2: Decay time (Jorge): t = LVET
								% 3: Decay time (Jorge): t = LVET + 1/3 t_dias
								% 4: Area method
								% 5: Two-area method
								% 6: Diastolic pressure method
								% 7: Pulse pressure method
								% 8: SV/PP
								% 9: 3-Wk optimisation
MethodPWV		= 0;%[1:5];			% 0: PWV = PWV,reference
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


%%	PATHS
PATHS.Root					= '~/Haemodynamic_Tools/Version6/';	% Path to the Windkessel study folder
PATHS.ParameterEstimation	= [PATHS.Root,'ParameterEstimation/'];
PATHS.Matlab				= [PATHS.ParameterEstimation,'Matlab/'];
PATHS.Output				= [PATHS.ParameterEstimation,'Output/',Dataset,'/'];
PATHS.Figures				= [PATHS.ParameterEstimation,'Figures/'];

addpath(PATHS.Matlab)


%%	REF data: choose virtual dataset
switch Dataset
	case {'Wk2', 'Wk3'}
		PATHS.Input = [PATHS.Root,'Windkessel/',Dataset,'_input/'];
		load([PATHS.Input,Dataset,'_reference.mat'],'REF')
		Scale_P_filt = 1;	%Dataset BP values given in mmHg
		
	case '55_art'
		PATHS.Input = [PATHS.Root,'VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_12_12/'];
		load([PATHS.Input,Dataset,'_reference.mat'],'REF')
		Scale_P_filt = 1;	%Dataset BP values given in mmHg
		
		addpath([PATHS.Root,'Others/ExtractReferenceBP/']);

		%	Format 55-artery dataset for Wk2/Wk3 estimation
		for jj = 1:length(REF)
			REF(jj).Dataset			= Dataset;
			REF(jj).Pressure		= REF(jj).L_Brachial.ONE_CYCLE.P;
			REF(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;
			REF(jj).t_in			= REF(jj).Asc_Ao.ONE_CYCLE.t;
			REF(jj).Q_in			= REF(jj).Asc_Ao.ONE_CYCLE.Q;
			REF(jj).SV				= REF(jj).HAEMO.SV;
			REF(jj).DBP				= min(REF(jj).Pressure);	% DBP
			REF(jj).SBP				= max(REF(jj).Pressure);	% SBP

			REF_temp(jj)			= ExtractReferenceBP(Scenario,REF(jj));
		end
		REF = REF_temp;
		clear REF_temp
		
	case {'pwdb_78', 'pwdb_4374'}
		PATHS.Input = [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Simulations/'];
		load([PATHS.Input,Dataset,'_reference.mat'],'REF')
		Scale_P_filt = 133.32;	%Dataset BP values given in Pa

		addpath([PATHS.Root,'Others/ExtractReferenceBP/'])
		
		%	Format Pete's dataset for Wk2/Wk3 estimation
		for jj = 1:length(REF)
			REF(jj).Dataset			= Dataset;
			REF(jj).INFO.ID			= REF(jj).INFO.ID;
			REF(jj).Pressure		= REF(jj).L_Brachial.ONE_CYCLE.P;	%[Pa]
			REF(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;		%[Pa]
			REF(jj).t_in			= REF(jj).Asc_Ao.ONE_CYCLE.t;
			REF(jj).Q_in			= REF(jj).Asc_Ao.ONE_CYCLE.Q;		%[m3/s]
			REF(jj).SV				= REF(jj).HAEMO.SV/10^6;			%[m3]
			REF(jj).LVET			= REF(jj).HAEMO.LVET;
			REF(jj).PWV				= REF(jj).PWV;
			REF(jj).C_T				= REF(jj).C_T_Wk3;
			REF(jj).DBP				= min(REF(jj).Pressure);	% DBP
			REF(jj).SBP				= max(REF(jj).Pressure);	% SBP

			REF_temp(jj)			= ExtractReferenceBP(Scenario,REF(jj));
		end
		REF = REF_temp; clear REF_temp
		
	otherwise
		error('Load a valid dataset')
end

%	Physiological filter based on the Normotensive and Hypertensive datasets
addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
REF = Physiological_BP_Filter(PATHS,REF,Scale_P_filt); clear REF_temp

		
%%	EST data: estimate CV parameters and calculate RMSE
addpath([PATHS.Root,'Others/LVET/'])
addpath([PATHS.Root,'Others/PulseAnalyse/'])
addpath([PATHS.Root,'Others/P_out/'])
addpath([PATHS.Root,'Others/PWV_TT/'])
addpath([PATHS.Root,'Others/ImpedanceAnalysis/'])
addpath([PATHS.Root,'Others/export_fig'])
addpath([PATHS.Root,'Others/P_0_iteration/'])

EST = struct;
k = 1;

tic

% LVET_list	= zeros(length(REF),max(MethodLVET));
load([PATHS.ParameterEstimation,'Output/LVET_4374'])

% PWV_list	= zeros(length(REF),max(MethodPWV));
load([PATHS.ParameterEstimation,'Output/PWV_4374'])

% P_out_list	= zeros(length(REF),max(MethodLVET),max(MethodP_out));
load([PATHS.ParameterEstimation,'Output/P_out_4374'])

% Z_0_list	= zeros(length(REF),max(MethodP_out),max(MethodPWV)+1,max(MethodZ_0));
load([PATHS.ParameterEstimation,'Output/Z_0_4374'])

%	Loop through each estimation method
for a1 = MethodLVET
	for a2 = MethodP_out
		for a3 = MethodR_T
			for a4 = MethodC_T
				for a6 = MethodZ_0
					if a6 == 4
						for a5 = MethodPWV
							%	Loop through every patient
							for jj = 1:length(REF)
								%	1. Estimation parameters
								EST(jj).ID				= REF(jj).INFO.ID;
								EST(jj).ID_num			= jj;
								EST(jj).Plots			= Plots_EST;
								
								%	2. Estimation methods
								EST(jj).MethodLVET		= a1;
								EST(jj).MethodP_out		= a2;
								EST(jj).MethodR_T		= a3;
								EST(jj).MethodC_T		= a4;
								EST(jj).MethodPWV		= a5;
								EST(jj).MethodZ_0		= a6;
								
% 								EST(jj).LVET			= LVET_list(jj,a1);
% 								EST(jj).PWV				= PWV_list(jj,a5);
% 								EST(jj).P_out			= P_out_list(jj,a1,a2);
% 								EST(jj).Z_0				= Z_0_list(jj,a2,a5,a6);
							end
							
							%	3. Estimate parameters
							for jj = 1:length(REF)
								EST_temp(jj) = ParameterEstimation(REF(jj),EST(jj),Scenario);
							end
							
								%	4. Extract parameter lists for LVET, PWV, P_out, and Z_0
% 							for jj = 1:length(REF)
% 								LVET_list(jj,a1)		= EST_temp(jj).LVET;
% 								PWV_list(jj,a5)			= EST_temp(jj).PWV;
% 								P_out_list(jj,a1,a2)	= EST_temp(jj).P_out;
% 								Z_0_list(jj,a2,a5,a6)	= EST_temp(jj).Z_0;	
% 							end
							
							EST = EST_temp;
							clear EST_temp
							
							%	Calculate pressure errors
							Methods = [num2str(a1),num2str(a2),num2str(a3),num2str(a4),'_',num2str(a5),'_',num2str(a6)];
							[ERR(k)] = Pressure_RMSE(REF,EST,Plots_EST,Methods);
							
							%	Update errors for every iteration
							save([PATHS.Output,'RMSE_',Dataset,'_',Methods], 'ERR')
							fprintf('Iteration %i complete; absolute and relative errors saved\n', k)
							k = k+1;
							
						end
						
					else
						a5 = 5;	%no PWV estimation
						
						%	Loop through every patient
						for jj = 1:length(REF)
							%	1. Estimation parameters
							EST(jj).ID				= REF(jj).INFO.ID;
							EST(jj).ID_num			= jj;
							EST(jj).Plots			= Plots_EST;
							
							%	2. Estimation methods
							EST(jj).MethodLVET		= a1;
							EST(jj).MethodP_out		= a2;
							EST(jj).MethodR_T		= a3;
							EST(jj).MethodC_T		= a4;
							EST(jj).MethodPWV		= a5;
							EST(jj).MethodZ_0		= a6;
							
% 							EST(jj).LVET			= LVET_list(jj,a1);
% 							EST(jj).P_out			= P_out_list(jj,a1,a2);
% 							EST(jj).Z_0				= Z_0_list(jj,a2,a5,a6);
						end
						
						%	3. Estimate parameters and pressure
						for jj = 1:length(REF)
							EST_temp(jj) = ParameterEstimation(REF(jj),EST(jj),Scenario);
						end
						
							%	4. Extract parameter lists for LVET, P_out, and Z_0
% 						for jj = 1:length(REF)
% 							LVET_list(jj,a1)		= EST_temp(jj).LVET;
% 							P_out_list(jj,a1,a2)	= EST_temp(jj).P_out;
% 							Z_0_list(jj,a2,a5,a6)	= EST_temp(jj).Z_0;
% 						end
						
						EST = EST_temp;
						clear EST_temp
						
						%	Calculate pressure errors
						Methods = [num2str(a1),num2str(a2),num2str(a3),num2str(a4),'_',num2str(a5),'_',num2str(a6)];
						[ERR(k)] = Pressure_RMSE(REF,EST,Plots_EST,Methods);
						
						%	Update errors for every iteration
						save([PATHS.Output,'RMSE_',Dataset,'_',Methods], 'ERR')
						fprintf('Iteration %i complete; absolute and relative errors saved\n', k)
						k = k+1;
						
					end
				end
			end
		end
	end
end

Comp_time = toc;

[~, b]			= min([ERR.RMS_abs_mean]);
best_method_abs	= ERR(b).ID;
best_RMS_abs	= ERR(b).RMS_abs_mean;

[~, b]			= min([ERR.RMS_rel_mean]);
best_method_rel	= ERR(b).ID;
best_RMS_rel	= ERR(b).RMS_rel_mean;
% 
% [~, RMS_index]	= sort([ERR.RMS_mean]);
% RMS_best_10		= [ERR(RMS_index(1:10)).RMS_mean];
% RMS_best_10_ID	= {ERR(RMS_index(1:10)).ID};


%%	Important - avoid using functions from different folder with same names
restoredefaultpath


%%	Function definitions
%	Length of 'ShortVector' becomes that of the other vector via spline interpolation
function [ShortVector] = SameVectorLength(Vector_1,Vector_2)

if length(Vector_1) < length(Vector_2)
	N_old = length(Vector_1);
	N_new = length(Vector_2);
	ShortVector = spline(linspace(0,1,N_old),Vector_1,linspace(0,1,N_new))';
elseif length(Vector_1) > length(Vector_2)
	N_old = length(Vector_2);
	N_new = length(Vector_1);
	ShortVector = spline(linspace(0,1,N_old),Vector_2,linspace(0,1,N_new))';
end

end

function [ERR] = Pressure_RMSE(REF,EST,Plots,Methods)

for jj = 1:length(EST)
	EST_temp(jj) = P_0_iteration(EST(jj),0,'Wk3');
end
EST = EST_temp;
clear EST_temp

%	Calculate pressure errors
for jj = 1:length(EST)
	%	RMS: the square root of the arithmetic mean of the squares of the values
	REF_P = REF(jj).Pressure_CBP;	
	EST_P = EST(jj).Pressure;
	
	if length(REF_P) > length(EST_P)
		%upsample EST_P
		t_short		= linspace(REF(jj).t_in(1),REF(jj).t_in(end),length(EST_P));
		EST_P		= spline(t_short,EST_P,REF(jj).t_in);
	elseif length(REF_P) < length(EST_P)
		%upsample REF_P
		t_short		= linspace(REF(jj).t_in(1),REF(jj).t_in(end),length(REF_P));
		REF_P		= spline(t_short,REF_P,REF(jj).t_in);
	end
	
	if Plots == 1
		figure, hold on
		plot(REF(jj).t_in,REF_P,'k','LineWidth',4), plot(EST(jj).t_in,EST_P,'--b','LineWidth',4)
		hold off
		legend('Reference cBP','Estimated cBP')
		xlabel('Time [s]')
		ylabel('Pressure [mmHg]')
		set(gca,'FontSize',30)
	end
	
	P_diff = EST_P - REF_P;
	PP_ref = (max(REF_P) - min(REF_P));
	
	RMS_wave_abs(jj) = sqrt( sum(P_diff.^2) / (length(P_diff)-1) );
	RMS_wave_rel(jj) = RMS_wave_abs(jj) / PP_ref;
end

ERR.ID									= Methods;

ERR.RMS_abs_mean						= mean(RMS_wave_abs);
ERR.RMS_abs_SD							= std(RMS_wave_abs);
[ERR.RMS_abs_max, ERR.RMS_abs_max_ID]	= max(RMS_wave_abs);
[ERR.RMS_abs_min, ERR.RMS_abs_min_ID]	= min(RMS_wave_abs);

ERR.RMS_rel_mean						= mean(RMS_wave_rel);
ERR.RMS_rel_SD							= std(RMS_wave_rel);
[ERR.RMS_rel_max, ERR.RMS_rel_max_ID]	= max(RMS_wave_rel);
[ERR.RMS_rel_min, ERR.RMS_rel_min_ID]	= min(RMS_wave_rel);

end

