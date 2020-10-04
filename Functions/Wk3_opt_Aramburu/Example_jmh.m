load('Example_jmh')
Scale_Q = 10^6;
Scale_P = 133.32;

addpath('~/Haemodynamic_Tools/Version6/Others//Wk3_opt_Aramburu/')

%	Scale for numerical reasons within Jorge's algorithm
t = REF.Asc_Ao.ONE_CYCLE.t;			%[s]
Q = REF.Asc_Ao.ONE_CYCLE.Q*Scale_Q;	%[mL/s]
P = REF.Asc_Ao.ONE_CYCLE.P/Scale_P;	%[mmHg]
P_out = REF.P_out/Scale_P;			%[mmHg]
R_T	= REF.R_T/(Scale_P*Scale_Q);	%[mmHg*s/mL]
C_T	= REF.C_T*(Scale_P*Scale_Q);	%[mL/mmHg]

[Z_0,R_2,C_T]		= Optimise_Wk3(t,Q,P,P_out,R_T*0.90,C_T,1,1);	%[mm*s/mL]
Z_0					= Z_0*(Scale_P*Scale_Q);	%[Pa*s/m3]