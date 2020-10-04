%%	Estimate CBP waveform using a transfer function to a peripheral pressure waveform
%	Input
%	-Peripheral pressure waveform
%
%	Output
%	-CBP waveform
% 
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (29/11/18)
%	v1.1 (20/12/18) - BP B-A plots
%
%==========================================================================
clear;
close all;


%%	Parameters
sig_type = 'carotABP';
options.do_plot = 0;
options.retain_onset_time = 1;


%%	Normotensive reference dataset
%	Load peripheral waveform data
% load('/Users/joh15/Clinical_Data/Sam_Normotensive/Normotensive_jmh','REF')
% load('/Users/joh15/Clinical_Data/Sam_Hypertensive/Hypertensive_jmh','REF')
load('/Users/joh15/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_12_12/55_art_reference','REF')

%	Format data
for jj = 1:length(REF)
% 	S(jj).v		= REF(jj).Carotid.ONE_CYCLE.P;			%In vivo
	S(jj).v		= REF(jj).L_CommonCarotid.ONE_CYCLE.P;	%In silico
	S(jj).fs	= round(1 / (REF(jj).Asc_Ao.ONE_CYCLE.t(2) - REF(jj).Asc_Ao.ONE_CYCLE.t(1)));
end


%%	TF
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/TransferFunction/')

for jj = 1:length(REF)
	%	TF estimate for a given patient
	CBP_from_Carotid(jj) = per_cen_TF_v2(S(jj),sig_type,options);
	
	if options.do_plot == 1
		%	Compare to SphygmoCor results (if available)
		plot(REF(jj).Asc_Ao.ONE_CYCLE.t,REF(jj).Asc_Ao.ONE_CYCLE.P,'--b','LineWidth',2)
		hold off
		axis tight
		legend('Carotid','TF','SphygmoCor')
	end	
end


%%	SBP errors
for jj = 1:length(REF)
	SBP_ref(jj)			= max(REF(jj).Asc_Ao.ONE_CYCLE.P);
	MBP_ref(jj)			= mean(REF(jj).Asc_Ao.ONE_CYCLE.P);
	DBP_ref(jj)			= min(REF(jj).Asc_Ao.ONE_CYCLE.P);
	PP_ref(jj)			= SBP_ref(jj) - DBP_ref(jj);
	SBP_est(jj)			= max(CBP_from_Carotid(jj).centrABP.v);
	MBP_est(jj)			= mean(CBP_from_Carotid(jj).centrABP.v);
	DBP_est(jj)			= min(CBP_from_Carotid(jj).centrABP.v);
	PP_est(jj)			= SBP_est(jj) - DBP_est(jj);
	
	ERR(jj).SBP			= SBP_est(jj) - SBP_ref(jj);
% 	ERR(jj).TF_Carotid	= max(CBP_from_Carotid(jj).centrABP.v) - max(REF(jj).L_CommonCarotid.ONE_CYCLE.P);
end

%	Mean and SD
FT_FT_mean		= mean([ERR.SBP]);
FT_FT_SD		= std([ERR.SBP]);
% FT_Carotid_mean	= mean([ERR.TF_Carotid]);
% FT_Carotid_SD	= std([ERR.TF_Carotid]);

fprintf('SBP error: TF to TF (mean +/- SD)	= %.2f +/- %.2f\n', FT_FT_mean, FT_FT_SD)
% fprintf('SBP error: TF to Carotid (mean +/- SD)	= %.2f +/- %.2f\n', FT_Carotid_mean, FT_Carotid_SD)


%%	PATHS
PATHS.Figures				= '/Users/joh15/Haemodynamic_Tools/Version6/Others/TransferFunction/Figures/';
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotBlandAltman')
addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotSave')


%%	PLOTS
PLOTS.Format				= 1;

%	Others
PLOTS.PUnits				= 133.322368;	% From mmHg to Pa

%	Fontsize
PLOTS.FontSize				= 4;	% 1: Use FontSizeS for the whole figure
									% 2: Use FontSizeM for the whole figure
									% 3: Use FontSizeL for the whole figure
									% 4: Use FontSizeXL for the whole figure
PLOTS.FontSizeXS			= 10;
PLOTS.FontSizeS				= 16;
PLOTS.FontSizeM				= 20;
PLOTS.FontSizeL				= 26;
PLOTS.FontSizeXL			= 36;
PLOTS.LineWidthS			= 1;
PLOTS.LineWidthM			= 2;
PLOTS.LineWidthL			= 4;

%	Title and other text
PLOTS.Title			= 1;	% 1: Display title on each figure
PLOTS.Legend		= 0;	% 1: Display Legend_text on each figure
PLOTS.XLabel		= 0;	% 1: Display x-axis labels
PLOTS.YLabel		= 0;	% 1: Display y-axis labels
PLOTS.Annotation	= 0;	% 1: Display annotations

%	Axis properties	
PLOTS.XLim			= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.YLim			= 1;	% 0: no limits, no axis
								% 1: 'custom'
								% 2: 'custom tight'
PLOTS.Grid			= 1;	% 1: grid on
PLOTS.ErrorBox		= 0;	% 1: Display absolute error box
								% 2: Display relative error box

PLOTS.Legend_text	= {'Physical','Optimised'};
PLOTS.XLabel_text	= {'Mean value(mmHg)'};
PLOTS.YLabel_text	= {'Estimation - Reference (mmHg)'};
PLOTS.XLim_values	= [0 250];
PLOTS.YLim_values	= [-10 10];
PLOTS.XTick			= 50;
PLOTS.YTick			= 5;
PLOTS.LineWidthS	= 1;
PLOTS.LineWidthM	= 2;
PLOTS.FontSizeL		= 16;
PLOTS.Markers		= 'o';

%	B-A SBP
PLOTS.Title_text	= 'Systolic CBP';
PLOTS.PlotName		= '55art_sCBP';	%[SBP, MBP, DBP, PP]

PlotBlandAltman(SBP_ref,SBP_est,PLOTS)
PlotSave(PATHS.Figures,PLOTS.PlotName)

%	B-A MBP
PLOTS.Title_text	= 'Mean CBP';
PLOTS.PlotName		= '55art_mCBP';	%[SBP, MBP, DBP, PP]

PlotBlandAltman(MBP_ref,MBP_est,PLOTS)
PlotSave(PATHS.Figures,PLOTS.PlotName)

%	B-A DBP
PLOTS.Title_text	= 'Diastolic CBP';
PLOTS.PlotName		= '55art_dCBP';	%[SBP, MBP, DBP, PP]

PlotBlandAltman(DBP_ref,DBP_est,PLOTS)
PlotSave(PATHS.Figures,PLOTS.PlotName)

%	B-A PP
PLOTS.Title_text	= 'Pulse CBP';
PLOTS.PlotName		= '55art_pCBP';	%[SBP, MBP, DBP, PP]

PlotBlandAltman(PP_ref,PP_est,PLOTS)
PlotSave(PATHS.Figures,PLOTS.PlotName)

%	Max SBP error
PLOTS.PlotName = '55art_Highest_SBP_error';

[max_err, max_loc] = max(abs([ERR.SBP]));
figure
hold on
plot(REF(max_loc).Asc_Ao.ONE_CYCLE.t,REF(max_loc).Asc_Ao.ONE_CYCLE.P,'k')
t_TF_max = 0:1/S(max_loc).fs:(length(CBP_from_Carotid(max_loc).centrABP.v) - 1)*1/S(max_loc).fs;
plot(t_TF_max,CBP_from_Carotid(max_loc).centrABP.v,'--b')
hold off
ylim([40 160])
xlabel('Time [s]')
ylabel('Pressure [mmHg]')
title('Highest SBP error')
set(gca,'FontSize',24)

PlotSave(PATHS.Figures,PLOTS.PlotName)

%	Min SBP error
PLOTS.PlotName = '55art_Lowest_SBP_error';

[min_err, min_loc] = min(abs([ERR.SBP]));
figure
hold on
plot(REF(min_loc).Asc_Ao.ONE_CYCLE.t,REF(min_loc).Asc_Ao.ONE_CYCLE.P,'k')
t_TF_min = 0:1/S(min_loc).fs:(length(CBP_from_Carotid(min_loc).centrABP.v) - 1)*1/S(min_loc).fs;
plot(t_TF_min,CBP_from_Carotid(min_loc).centrABP.v,'--b')
hold off
ylim([40 160])
xlabel('Time [s]')
ylabel('Pressure [mmHg]')
title('Lowest SBP error')
set(gca,'FontSize',24)

PlotSave(PATHS.Figures,PLOTS.PlotName)


%%
restoredefaultpath
