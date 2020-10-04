function [ERRORS] = CalculateBPErrors(REF,EST)
for jj = 1:length(REF)
	ERRORS(jj).SBP_abs	= EST(jj).SBP - REF(jj).SBP;
	ERRORS(jj).MBP_abs	= EST(jj).MBP - REF(jj).MBP;
	ERRORS(jj).DBP_abs	= EST(jj).DBP - REF(jj).DBP;
	ERRORS(jj).PP_abs	= EST(jj).PP - REF(jj).PP;
	
	ERRORS(jj).SBP_rel	= 100*ERRORS(jj).SBP_abs / REF(jj).SBP;
	ERRORS(jj).MBP_rel	= 100*ERRORS(jj).MBP_abs / REF(jj).MBP;
	ERRORS(jj).DBP_rel	= 100*ERRORS(jj).DBP_abs / REF(jj).DBP;
	ERRORS(jj).PP_rel	= 100*ERRORS(jj).PP_abs / REF(jj).PP;
end
end