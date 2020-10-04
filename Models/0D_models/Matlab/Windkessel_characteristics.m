%%	Windkessel dataset characteristics
%	Input
%
%	Output
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (24/08/19) - based on Hypertensive_Characteristics
%
%==========================================================================
clear;
close all;


%%	.mat data
Wk = 'Wk3';
load(['~/Haemodynamic_Tools/Version6/Windkessel/',Wk,'_input/',Wk,'_reference_v2.mat'])


%%	Read data
for jj = 1:length(REF)
	SBP_c(jj) = max(REF(jj).Pressure);
	MBP_c(jj) = mean(REF(jj).Pressure);
	DBP_c(jj) = min(REF(jj).Pressure);
	PP_c(jj) = SBP_c(jj) - DBP_c(jj);
	
	REF(jj).t			= REF(jj).t_in;	%[s]
	dt					= REF(jj).t(2) - REF(jj).t(1);
	
	REF(jj).Q			= REF(jj).Q_in;
	
	SV(jj)		= REF(jj).SV;	%[mL]
	HR(jj)		= REF(jj).HR;	%bpm
	CO(jj)		= SV(jj)*HR(jj)/1000;	%[L/min]
end

%%	Calculate mean and standard deviation for entire population
%	Central BP mean/SD values
SBP_mean_c	= mean(SBP_c);
SBP_SD_c	= std(SBP_c);
MBP_mean_c	= mean(MBP_c);
MBP_SD_c	= std(MBP_c);
DBP_mean_c	= mean(DBP_c);
DBP_SD_c	= std(DBP_c);
PP_mean_c	= mean(PP_c);
PP_SD_c		= std(PP_c);

%	Cardiac values
SV_mean		= mean(SV);	% Stroke volume [mL]
SV_SD		= std(SV);
HR_mean		= mean(HR);	% Heart rate [bpm]
HR_SD		= std(HR);
CO_mean		= mean(CO);
CO_SD		= std(CO);

% 	Arterial
% PVR_mean	= mean([REF.R_T]);
% PVR_SD		= std([REF.R_T]);
% tau_mean	= mean([REF.tau]);
% tau_SD		= std([REF.tau]);

% PWV_mean	= mean(PWV);
% PWV_SD		= std(PWV);


%%	Calculate max and min BP values for entire population
%	Central BP values
SBP_max_c = max(SBP_c);
SBP_min_c = min(SBP_c);
MBP_max_c = max(MBP_c);
MBP_min_c = min(MBP_c);
DBP_max_c = max(DBP_c);
DBP_min_c = min(DBP_c);
PP_max_c = max(PP_c);
PP_min_c = min(PP_c);


%%	Generate .txt file with mean +/- SD
fileID = fopen([Wk,'_characteristics.txt'],'w');
fprintf(fileID,[Wk,' characteristics:\n']);
fprintf(fileID,'Population size: %i\n',length(REF));

fprintf(fileID,['Mean +/- standard deviation\n']);
fprintf(fileID,'DBP_c [mmHg] & %3.1f $\\pm$ %3.1f\n',DBP_mean_c, DBP_SD_c);
fprintf(fileID,'MBP_c [mmHg] & %3.1f $\\pm$ %3.1f\n',MBP_mean_c, MBP_SD_c);
fprintf(fileID,'SBP_c [mmHg] & %3.1f $\\pm$ %3.1f\n',SBP_mean_c, SBP_SD_c);
fprintf(fileID,'PP_c [mmHg] & %3.1f $\\pm$ %3.1f\n',PP_mean_c, PP_SD_c);

fprintf(fileID,'SV [mL] & %3.1f $\\pm$ %3.1f\n',SV_mean, SV_SD);
fprintf(fileID,'HR [bpm] & %3.1f $\\pm$ %3.1f\n',HR_mean, HR_SD);
fprintf(fileID,'CO [L/min] & %3.1f $\\pm$ %3.1f\n',CO_mean, CO_SD);

% fprintf(fileID,'PWV [m/s] & %3.1f $\\pm$ %3.1f\n',PWV_mean, PWV_SD);
%
fclose(fileID);


