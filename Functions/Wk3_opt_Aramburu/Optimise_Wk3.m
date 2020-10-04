%%	Calculate optimal Wk3 parameters for given flow and pressure
%	Input
%	-t, Q_in, P_in: time, flow and pressure at given location
%	-P_out: outflow pressure
%	-P_0: initial pressure (commonly DBP)
%	-tau: time constant
%	-R_2: initial R_2 guess
%	-C: initial C guess
%	-plots: =1 generate plots	
%
%	Output
%	-R_1, R_2, C: optimised Wk3 parameter values
%
%==========================================================================
%	Jorge Aramburu, King's College London (original author)
%	v1.0 (13/09/18)
%	Jorge Mariscal-Harana, King's College London
%	v2.0 (03/12/18) - Clean up code
%
%==========================================================================

function [R_1, R_2, C, P_est] = Optimise_Wk3(t,Q_in,P_in,P_out,R_2,C,plots,fig_num)

Est_old = [R_2; C];		%Initial guess for R2 y C
Lim_t = [0 R_2; 0 C];
Lim   = [0 R_2; 0 C];

[R_1, R_2, C, Est_hist_Wk, Cost, fig_num] = Estimate_R2_and_C(t, P_in, Q_in, P_out, Est_old, Lim_t, Lim, plots, fig_num);

if plots == 1
% 	[meanP_in1, meanP_est1, fig_num] = Plot_Estimated_Windkessel(t, P_in, Q_in, R_1, R_2, C, P_out, fig_num);
	[meanP_in2, meanP_est2, fig_num, P_est] = Plot_Estimated_WindkesselB(t, P_in, Q_in, R_1, R_2, C, P_out, fig_num);
end

end