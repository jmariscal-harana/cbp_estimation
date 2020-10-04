%%	Estimate CV parameters
%	Input
%	-PATHS: paths to required folders
%	-REF: reference in silico data
%	-EST: estimated parameters
%
%	Output
%	-EST: parameter estimation data
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (22/05/19) - based on EstimationParameters
%	v1.1 (13/06/19) - REF case outside individual functions; added PATHS
%
%==========================================================================

function [EST] = ParameterEstimation(PATHS,REF,EST,Plots)
EST.Plots = Plots;


%%	EST inputs
%	Flow and time from REF data
EST.t_in = REF.t_in;
EST.Q_in = REF.Q_in;

%	Start (t_0) and final (t_f) time, and cardiac period (T)
EST.t_0	= EST.t_in(1);
EST.t_f	= EST.t_in(end);
EST.T	= EST.t_f - EST.t_0; if EST.T > 2.5 || EST.T < 0.25, error('Cardiac period is %f [s]', EST.T), end

%	P_0 = DBP (run multiple cycles to reach quasi-steady state)
EST.P_0 =  REF.DBP;


%%	Estimate LVET from pressure or flow waveform
if EST.MethodLVET == 0
	EST.LVET = REF.LVET;
else
	EST	= EstimationLVET(PATHS,REF,EST);
end


%%	Estimate P_out from whole waveform or individual BP values
if EST.MethodP_out == 0
	EST.P_out = REF.P_out;
else
	EST	= EstimationPout(PATHS,REF,EST);
end


%%	Estimate R_T from REF pressure and flow data
if EST.MethodR_T == 0
	EST.R_T = REF.R_T;
else
	EST = EstimationR_T(REF,EST);
end


%%	Estimate PWV from pressure or flow waveforms
if EST.MethodPWV == 0
	EST.PWV = REF.PWV;
else
	EST = EstimationPWV(PATHS,REF,EST);
end


%%	Estimate Z_0 from REF pressure and flow data
Estimate_Z_0 = 1;
if isfield(EST,'Windkessel')
	if strcmp(EST.Windkessel,'Wk2')
		EST.Z_0 = 0;
		Estimate_Z_0 = 0;
	end
end
if Estimate_Z_0 == 1
	if EST.MethodZ_0 == 0
		EST.Z_0 = REF.Z_0;
	else
		EST = EstimationZ_0(PATHS,REF,EST);
	end
end


%%	Estimate C_T from whole waveform or individual BP values
if EST.MethodC_T == 0
	EST.C_T = REF.C_T;
else
	EST = EstimationC_T(PATHS,REF,EST);
end


%%	Estimate tau from R_T, C_T, Z_0
EST = EstimationTau(EST);


%%	Ad-hoc plots
% EST = P_0_iteration(EST,0);
% figure
% hold on
% plot(REF.Pressure_CBP,'k')
% plot(REF.Pressure,'--k')
% plot(EST.Pressure,'--b')
% hold off
% legend('CBP','Brachial','Estimation')


end

