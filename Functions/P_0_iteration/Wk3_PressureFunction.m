%%	Calculates 3-element Windkessel pressure
%	Input:
%	-Flow, time
%	-Z, R2, C values
%	-P_0, P_out
%
%	Output:
%	-P_wave
%	-P_components
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (11/12/17) - based on Wk2_PressureFunction.m
%	v1.1 (15/01/18) - using tau instead of C_T for the calculation of P_Wk3
%
%==========================================================================

function [P_Wk3,P_Wk3_exponential,P_Wk3_integral,P_Wk3_exp_int] = Wk3_PressureFunction(P_out,P_0,R_T,R_1,C_T,Q_in,t_in)

tau		= (R_T - R_1)*C_T; 
t_in	= linspace(t_in(1),t_in(end),length(Q_in))';
t_0		= t_in(1);
dt		= t_in(2) - t_in(1);

P_Wk3 = P_out + (P_0 - P_out - R_1*Q_in(1)) * exp((t_0-t_in)/tau) + R_1*Q_in + exp(-t_in/tau) / C_T .* cumtrapz(Q_in.*exp(t_in/tau))*dt;

% P_Wk3_exponential = P_out + (P_0 - P_out - R_1*Q_in(1)) * exp((t_0-t_in)/tau);
% P_Wk3_integral = cumtrapz(Q_in.*exp(t_in/tau))*dt;
% P_Wk3_exp_int = exp(-t_in/tau) / C_T .* P_Wk3_integral;

% figure,hold on
% plot(t_in,P_Wk3,'--')
% plot(t_in,P_Wk3_exponential)
% plot(t_in,P_Wk3_integral)
% plot(t_in,P_Wk3_exp_int)
% plot(t_in,R_1*Q_in)
% hold off
% % 
% legend('P', '1st exp', 'integral', '2nd exp * integral', 'R_1 * Q_{in}')

end
