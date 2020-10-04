function [REF] = Wk_VirtualParameters(PATHS)
%%	Virtual populations: baseline vascular parameters (mean and SD where possible)
REF.LIT.P_out_baseline		= 33.2;						% mmHg;	(1) capillary pressure results from http://www.ncbi.nlm.nih.gov/pubmed/8458818
REF.LIT.P_out_SD			= 1.5*sqrt(13);				% (1) capillary pressure results from http://www.ncbi.nlm.nih.gov/pubmed/8458818

%
REF.LIT.P_0_baseline		= REF.LIT.P_out_baseline;	% mmHg;

%
R_T							= [0.80, 1.42, 1.8, 0.93];	% [mmHg*s/mL]; (1) from Table 6 - 10.1007/s10439-015-1313-8
														%	WARNING: UNITS are mmHg*s*m2/mL!!! (2) from Table 1, 12 normal subjects (37+/-2 yo) - http://ajpheart.physiology.org/content/237/5/H550
														% (3) approximated from Table 1 - 10.1161/CIRCULATIONAHA.105.541078 - data for >45 years old
														% (4) from Figure 7 - 10.1016/0002-9149(85)90659-9 - calculated as average from data for 20-60 years old
														%	NOTE: convert from dyne*s*m2/cm5 to mmHg*s/mL using 1/(133.322368*10)

R_T_SD						 = [0.09];					%	WARNING: UNITS are mmHg*s*m2/mL!!! (1) from Table 1, 12 normal subjects (37+/-2 yo) - http://ajpheart.physiology.org/content/237/5/H550

[R_T(2),R_T_SD(1)] = R_T_correction_P_out(REF, R_T(2), R_T_SD(1));	% Correction for R_T(2) and R_T_SD(1) for body surface taking into account P_out

REF.LIT.R_T_baseline		= R_T(2);
REF.LIT.R_T_SD				= R_T_SD(1);					

%
C_T							= [1.7, 0.685, 1.26];		% (1) from Table 3 - DOI: 10.1007/s10439-015-1313-8, and this from Figure 1 - 10.1161/01.HYP.34.4.808; calculated as SV/PP_brachial, which may lead to underestimation of C_T,central
														% mL/mmHg; (2) average from Table 2 - DOI: 10.1042/cs0950669; calculated from P_sys,end - P_dias
														%	WARNING: UNITS are mL/(mmHg*m2)!!! (3) from Figure 5, 12 normal subjects (37+/-2 yo) - http://ajpheart.physiology.org/content/237/5/H550

C_T_SD						= [0, 0.25, 0.04];			% (1) not available
														% mL/mmHg; (2) average from Table 2 - DOI: 10.1042/cs0950669; calculated from P_sys,end - P_dias
														%	WARNING: UNITS are mL/(mmHg*m2)!!! (3) from Figure 5, 12 normal subjects (37+/-2 yo) - http://ajpheart.physiology.org/content/237/5/H550

[C_T(3), C_T_SD(3)] = C_T_correction(C_T(3), C_T_SD(3));			% Correction for C_T(3) and C_T_SD(3) for body surface 

REF.LIT.C_T_baseline	= C_T(3);
REF.LIT.C_T_SD			= C_T_SD(3);

%	Wk3 parameters
Z_0					= [];								% dyne*s/cm^5 (1) 10.1016/0002-9149(85)90659-9
														

PWV					= [4.6];					% m/s; (1) from abstract - DOI: 10.1016/j.jcmg.2010.09.016
PWV_SD				= [1.5];					% m/s; (1) from abstract - DOI: 10.1016/j.jcmg.2010.09.016
REF.LIT.PWV_baseline	= PWV(1);
REF.LIT.PWV_SD			= PWV_SD(1);

%
Diam				= [0.031];					% m; (1) systolic diameter, from abstract - DOI: 10.1016/j.jcmg.2010.09.016
Diam_SD				= [0.004];					% m; (1) systolic diameter, from abstract - DOI: 10.1016/j.jcmg.2010.09.016
REF.LIT.Diam_baseline	= Diam(1);
REF.LIT.Diam_SD			= Diam_SD(1);

REF.LIT.Density			= 1060;						% kg/m3; from Caro CG et al. The mechanics of the circulation (1978)

%	Baseline cardiac parameters (mean and SD where possible)
run([PATHS.Wk_study,'/Inflow/Weissler_1961_Data'])

%
SV						= SV_normal;				% mL; (1) 17 normal subjects, Weissler 1961
REF.LIT.SV_baseline			= mean(SV);
REF.LIT.SV_SD				= std(SV);

%
HR						= HR_normal;				% bpm; (1) 17 normal subjects, Weissler 1961
REF.LIT.HR_baseline			= mean(HR);
REF.LIT.HR_SD				= std(HR);

end


%%	Function definitions
function [R_T, R_T_SD] = R_T_correction_P_out(REF, R_T, R_T_SD)
%%	Correction for when R_T,clinical = P_mean,clinical / Q_mean,clinical (does not take into account P_out)
%
%	1: R_T,target = (1 - P_out/P_mean,clinical) * R_T,clinical
%	2: R_T,target = R_T,clinical - P_out/Q_mean,clinical
%
Correction = 2;

P_out	= REF.LIT.P_out_baseline;

%	Data from Table 1, http://ajpheart.physiology.org/content/237/5/H550:
SBP		= 126;	% mmHg 
DBP		= 71;	% mmHg
%	If P_mean,clinical (MBP) is not available: MBP = DBP + 0.4*(SBP - DBP)
%	where DBP, SBP are diastolic and systolic blood pressure, respectively
MBP		= DBP + 0.4*(SBP - DBP);	% mmHg

Body_surface	= 1.8;	% m2
Cardiac_index	= 3832;	% mL/(min*m2)
Q_mean			= Cardiac_index*Body_surface/60;	% mL/s
R_T_surface		= R_T;

R_T				= R_T/Body_surface;	% R_T,clinical given in units mmHg*s*m2/mL

R_T = (MBP) / Q_mean;	

%	All methods yield similar results
switch Correction
	case 1
		R_T	= (1 - P_out/MBP)*R_T;
	case 2
		R_T	= R_T - P_out/Q_mean;
	case 3
		R_T = (MBP - P_out) / Q_mean;	
end

R_T_SD = R_T_SD/Body_surface;

end

function [C_T, C_T_SD] = C_T_correction(C_T, C_T_SD)
%%	Correction for body surface for http://ajpheart.physiology.org/content/237/5/H550
Body_surface	= 1.8;	% m2

C_T			= C_T*Body_surface;	% C_T,clinical given in units  mL/(mmHg*m2)
C_T_SD		= C_T_SD*Body_surface;

end
