%%	Extract BP values from waveform or from DBP/SBP
%	Input
%	-Scenario: waveform or DBP/SBP values
%	-REF: in vivo/in silico data
%
%	Output
%	-REF: DBP, MBP, SBP, PP
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (22/05/19) - based on Wk_ExtractReferenceBP
%
%==========================================================================

function [REF] = ExtractReferenceBP(Scenario,REF)
switch Scenario
	case 'Waveform'	% Reference BP waveform is available
		REF.DBP	= min(REF.Pressure);	% DBP
		REF.MBP	= mean(REF.Pressure);	% MBP
		REF.SBP	= max(REF.Pressure);	% SBP
		
	case 'BP'	% Only reference BP values are available
		REF.DBP	= REF.DBP;				% DBP
		REF.SBP	= REF.SBP;				% SBP
		REF.MBP = REF.DBP + 0.4*(REF.SBP - REF.DBP);

end
REF.PP = REF.SBP - REF.DBP;
end