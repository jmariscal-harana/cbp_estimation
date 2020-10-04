%%	Calculates 2-element Windkessel pressure
%	Input
%	-Flow, time
%	-R, C, tau values
%	-P_0, P_out
%
%	Output
%	-P_wave
%	-P_components
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (19/10/17) 
%	v1.1 (15/01/18) - using tau instead of C_T for the calculation of P_Wk2
%	v2.0 (18/04/19) - input format changed from INPUT. to individual 
%		parameters
%
%==========================================================================

function [P_Wk2,P_Wk2_exponential,P_Wk2_integral,P_Wk2_exp_int] = Wk2_PressureFunction(P_out,P_0,R_T,C_T,Q_in,t_in)

tau		= R_T*C_T;
t_in	= linspace(t_in(1),t_in(end),length(Q_in))';
t_0		= t_in(1);
dt		= t_in(2) - t_in(1); 

P_Wk2 = P_out + (P_0 - P_out) * exp((t_0 - t_in)/tau) + exp(-t_in/tau) / C_T .* cumtrapz(Q_in.*exp(t_in/tau))*dt;

% P_Wk2_exponential = (P_0 - P_out) * exp((t_0 - t_in)/tau);
% P_Wk2_integral = cumtrapz(Q_in.*exp(t_in/tau))*dt;
% P_Wk2_exp_int = exp(-t_in/tau) / C_T .* P_Wk2_integral;

end
