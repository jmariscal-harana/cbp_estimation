%%	Physiological BP filter extracted from Normotensive/Hypertensive datasets
%	Input
%	-REF: reference data
%	-Scale_P: 
%
%	Output
%	-REF_filtered: filtered reference data
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (22/05/19)
%
%==========================================================================

function [REF_filtered] = Physiological_BP_Filter(PATHS,REF)
addpath([PATHS.Root,'Others/ExtractReferenceBP/'])
Scale_P = 133.32;

%1. Extract DBP, SBP, PP
for jj = 1:length(REF)
	REF_filt.Pressure	= REF(jj).Pressure_CBP;
	REF_temp(jj)		= ExtractReferenceBP('Waveform',REF_filt);
end

%2. Apply physiological BP filter extracted from Normotensive/Hypertensive datasets and save as REF_filtered
warning('Physiological BP filter: SI units required for input BP values.')
if mean([REF_temp.SBP]) < 10*Scale_P || mean([REF_temp.SBP]) > 500*Scale_P, error('The mean SBP for this dataset is %f',mean([REF_temp.SBP])), end

kk = 0;
for jj = 1:length(REF)
	if REF_temp(jj).SBP > 220*Scale_P || REF_temp(jj).DBP < 44*Scale_P || REF_temp(jj).PP < 18*Scale_P || REF_temp(jj).PP > 109*Scale_P
		%non-physiological virtual patient
	else
		%physiological virtual patient
		kk = kk + 1;
		REF(jj).ID_original = jj;
		REF_filtered(kk) = REF(jj);
	end
end

fprintf('Filtered (discarded) subjects: %i out of %i.\n',length(REF) - kk, length(REF)) 

end
