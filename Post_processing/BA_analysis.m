%%	Generate Bland-Altman plots for the CV parameter estimation study
%	Input
%	-REF: parameters used to generate virtual patients
%	-EST: methods for parameter estimation for virtual patients
%	-PATHS: folders where required data is stored
%
%	Output
%	B-A plots
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (20/07/19)
%
%==========================================================================
clear;
close all;


%%	REF input parameters
REF_Dataset			= 'pwdb_4374';	% Reference dataset (virtual): 'Wk2', 'Wk3', '55_art', 'pwdb_78', 'pwdb_4374'
EST_Dataset			= 'pwdb_4374_carotid';


%%	EST input parameters
Scenario		= 'Waveform';	% 'Waveform': pressure waveform available
								% 'BP': only DBP/MBP/SBP values available
Plots_EST		= 0;			% 1: show estimation method plots

%	Estimation methods
MethodLVET		= 5;		% 0: LVET = LVET,reference
								% 1: LVET from P using Sam's algorithm
								% 2: LVET from P using Pete's PPG algorithm
								% 3: LVET from P using WF								
								% 4: LVET = 0.37*sqrt(T)
								% 5: Q=0 or max(Q) after min Q
MethodP_out		= 4;%[1:5];		% 0: P_out = P_out,reference
								% 1: P_out = 0
								% 2: P_out = P_dias/2
								% 3: P_out = 0.7 * DBP (Parragh 2014)
								% 4: Decay time (Jorge): t = LVET
								% 5: Decay time (Jorge): t = LVET + 1/3 t_dias
MethodR_T		= 1;			% 0: R_T = R_T,reference
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
MethodPWV		= 1;%[1:5];			% 0: PWV = PWV,reference
								% 1: TT: ascending-descending aorta Q
								% 2: TT: carotid-femoral P
								% 3: Least-squares: ascending-descending aorta Q
								% 4: Least-squares: carotid-femoral P
								% 5: Sum of squares
MethodZ_0		= 17;			% 0: Z_0 = Z_0,reference
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
switch REF_Dataset
	case {'Wk2','Wk3'}
		PATHS.Input			= [PATHS.Root,'Windkessel/',REF_Dataset,'_input/'];
	case 'pwdb_4374'
		PATHS.Input			= [PATHS.Root,'VirtualDB/pwdb/Nektar_outputFiles/Simulations/'];
end
PATHS.Output				= [PATHS.ParameterEstimation,'Output/',EST_Dataset,'/'];
PATHS.Figures				= [PATHS.ParameterEstimation,'Figures/'];

addpath(PATHS.Matlab)


%%	Load CV parameter estimation data
%	Loop through each estimation method
load([PATHS.Input,REF_Dataset,'_reference_v2.mat'],'REF')

%	Physiological filter based on the Normotensive and Hypertensive datasets
addpath([PATHS.Root,'Others/Physiological_BP_Filter/'])
switch REF_Dataset
	case {'Wk2','Wk3'}
		for jj = 1:length(REF)
			REF(jj).Pressure		= REF(jj).Pressure*133.322368;
			REF(jj).Pressure_CBP	= REF(jj).Pressure;
		end
% 		REF = Physiological_BP_Filter(PATHS,REF);
		
		for jj = 1:length(REF)
			REF_temp.LVET(jj) = REF(jj).LVET*10^3;
			REF_temp.P_out(jj) = REF(jj).P_out;
			REF_temp.R_T(jj) = REF(jj).R_T;
			REF_temp.C_T(jj) = REF(jj).C_T;
			REF_temp.Z_0(jj) = REF(jj).Z_0;
		end
		REF = REF_temp; clear REF_temp
		
		k = 1;
		for a1 = MethodLVET
			Methods = [num2str(a1),'000_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).LVET = [EST.LVET]*10^3;
			k = k+1;
		end
		k = 1;
		for a2 = MethodP_out
			Methods = ['0',num2str(a2),'00_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).P_out = [EST.P_out]/133.322368;
			k = k+1;
		end
		k = 1;
		for a3 = MethodR_T
			Methods = ['00',num2str(a3),'0_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).R_T = [EST.R_T]/133.322368/10^6;
			k = k+1;
		end
		k = 1;
		for a4 = MethodC_T
			Methods = ['000',num2str(a4),'_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).C_T = [EST.C_T]*(133.322368*10^6);
			k = k+1;
		end
		switch REF_Dataset
			case 'Wk3'
				k = 1;
				for a6 = MethodZ_0
					Methods = ['0000_0_',num2str(a6)];
					load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
					EST_temp(k).Z_0 = [EST.Z_0]/133.322368/10^6;
					k = k+1;
				end
		end
		
	case 'pwdb_4374'
		for jj = 1:length(REF)
			REF(jj).Pressure		= REF(jj).Asc_Ao.ONE_CYCLE.P;
			REF(jj).Pressure_CBP	= REF(jj).Asc_Ao.ONE_CYCLE.P;
		end
		REF = Physiological_BP_Filter(PATHS,REF);
		
		for jj = 1:length(REF)
			REF_temp.LVET(jj) = REF(jj).HAEMO.LVET*10^3;
			REF_temp.P_out(jj) = REF(jj).P_out/133.322368;
			REF_temp.R_T(jj) = REF(jj).R_T/133.322368/10^6;
			REF_temp.C_T(jj) = REF(jj).C_T*(133.322368*10^6);
			REF_temp.PWV(jj) = REF(jj).PWV;
			REF_temp.Z_0(jj) = REF(jj).Z_0/133.322368/10^6;
		end
		REF = REF_temp; clear REF_temp
		
		k = 1;
		for a1 = MethodLVET
			Methods = [num2str(a1),'000_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).LVET = [EST.LVET]*10^3;
			k = k+1;
		end
		k = 1;
		for a2 = MethodP_out
			Methods = ['0',num2str(a2),'00_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).P_out = [EST.P_out]/133.322368;
			k = k+1;
		end
		k = 1;
		for a3 = MethodR_T
			Methods = ['00',num2str(a3),'0_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).R_T = [EST.R_T]/133.322368/10^6;
			k = k+1;
		end
		k = 1;
		for a4 = MethodC_T
			Methods = ['000',num2str(a4),'_0_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).C_T = [EST.C_T]*(133.322368*10^6);
			k = k+1;
		end
		k = 1;
		for a5 = MethodPWV
			Methods = ['0000_',num2str(a5),'_0'];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).PWV = [EST.PWV];
			k = k+1;
		end
		k = 1;
		for a6 = MethodZ_0
			Methods = ['0000_0_',num2str(a6)];
			load([PATHS.Output,'EST_',EST_Dataset,'_',Methods],'EST')
			EST_temp(k).Z_0 = [EST.Z_0]/133.322368/10^6;
			k = k+1;
		end
end
EST = EST_temp; clear EST_temp

% PWV_outliers = find (([EST.PWV] - [REF.PWV]) > 3.1);


%%	PLOTS general format
addpath([PATHS.Root,'/Others/PlotBlandAltman'])
addpath([PATHS.Root,'/Others/PlotFormat'])
addpath([PATHS.Root,'/Others/PlotSave'])

PLOTS.Format				= 1;
PLOTS.BA					= 1;

%	Fontsize
PLOTS.FontSize				= 45;
PLOTS.LegendSize			= 50;
PLOTS.FontSizeAnno			= 40;
PLOTS.LineWidth				= 2;
PLOTS.MarkerSize			= 20;

%	Title, legend, axis labels, and annotations
PLOTS.Title				= 1;	% 1: Display title on each figure
PLOTS.Legend			= 0;	% 1: Display Legend_text on each figure
PLOTS.XLabel			= 1;	% 1: Display x-axis labels
PLOTS.YLabel			= 1;	% 1: Display y-axis labels
PLOTS.Annotation		= 2;	% 1: Display annotations

%	Axis properties	
PLOTS.XLim				= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.YLim				= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.Grid				= 1;	% 1: grid on

%	Error box
PLOTS.ErrorBox			= 0;	% 1: Display absolute error box
								% 2: Display relative error box

%	Text
PLOTS.Legend_text		= {''};
PLOTS.XLabel_text_all	= {'Reference [ms]','Reference [mmHg]',	'Reference [mmHg s/mL]',...
	'Reference [mL/mmHg]','Reference [m/s]','Reference [mmHg s/mL]'};
PLOTS.YLabel_text_all	= {'Est. - Ref. [ms]','Est. - Ref. [mmHg]','Est. - Ref. [mmHg s/mL]',...
	'Est. - Ref. [mL/mmHg]','Est. - Ref. [m/s]','Est. - Ref. [mmHg s/mL]'};

switch REF_Dataset
	case 'Wk2'
		CV						= {'LVET','P_out','R_T','C_T'};
		PLOTS.Title_text_all	= {'LVET','P_{out}','R_T','C_T'};
		
		XLim_values				= {[200 350],	[25 40],	[0.4 0.6],	[2.0 2.6]	};
		XTick					= {50,			5,			0.1,		0.2,		};
		YLim_values				= {[-2 2],		[-25 25],	[-0.1 0.1],	[-1.5 1.5]	};
		YTick					= {1,			10,			0.05,		1,			};
		
		PLOTS.BA_sf_all			= [1, 1, 3, 1];
		
	case 'Wk3'
		CV						= {'LVET','P_out','R_T','C_T','Z_0'};
		PLOTS.Title_text_all	= {'LVET','P_{out}','R_T','C_T','Z_0'};
		
		XLim_values				= {[200 350],	[25 40],	[0.4 0.6],	[2.0 2.6],		[0 0.1]};
		XTick					= {50,			5,			0.1,		0.2,			0.02};
		YLim_values				= {[-2 2],		[-25 25],	[-0.1 0.1],	[-0.02 0.02],	[-0.07 0.07]};
		YTick					= {1,			10,			0.05,		0.01,			0.025};
		
		PLOTS.BA_sf_all			= [1, 1, 3, 3, 3];

	case 'pwdb_4374'
% 		CV						= {'LVET','P_out','R_T','C_T','PWV','Z_0'};
% 		PLOTS.Title_text_all	= {'LVET','P_{out}','R_T','C_T','PWV','Z_0'};
		CV						= {'PWV'};
		PLOTS.Title_text_all	= {'PWV'};
		
% 		XLim_values				= {[250 400],	[32 34],	[0 2],		[0 3.5],	[0 20],		[0 0.07]};
% 		XTick					= {50,			1,			1,			1,			5,			0.02};
% 		YLim_values				= {[-10 10],	[-15 15],	[-0.2 0.2],	[-2 2],		[-6 6],		[-0.13 0.13]};
% 		YTick					= {5,			10,			0.1,		1,			3,			0.05};
		XLim_values				= {[0 20]};
		XTick					= {5};
		YLim_values				= {[-6 6]};
		YTick					= {3};

		
% 		PLOTS.BA_sf_all			= [1, 1, 3, 1, 1, 3];
		PLOTS.BA_sf_all			= [1];

end

Markers					= {'ko'};

for jj = 1:length(CV)
	
	PLOTS.Title_text	= PLOTS.Title_text_all{jj};
	PLOTS.XLabel_text	= PLOTS.XLabel_text_all{jj};
	PLOTS.YLabel_text	= PLOTS.YLabel_text_all{jj};
	PLOTS.XLim_values	= XLim_values{jj};
	PLOTS.YLim_values	= YLim_values{jj};
	PLOTS.XTick			= XTick{jj};
	PLOTS.YTick			= YTick{jj};
	PLOTS.CV			= CV{jj};
	PLOTS.BA_sf			= PLOTS.BA_sf_all(jj);
	
	PLOTS.Markers		= Markers;

	BlandAltman_CV_parameters(PATHS,PLOTS,REF,EST,EST_Dataset,Sc)
	
end

clear


%%	Function definitions
function BlandAltman_CV_parameters(PATHS,PLOTS,REF,EST,EST_Dataset,Sc)
PLOTS.PlotName	= [EST_Dataset,'_',PLOTS.CV];
eval(['PlotBlandAltman_v2([REF.',PLOTS.CV,'],[EST.',PLOTS.CV,'],PLOTS);']);
PlotSave(PATHS.Figures,[PLOTS.PlotName,'_B-A_',Sc])
end
