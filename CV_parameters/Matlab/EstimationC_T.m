%%	Estimate C_T from whole waveform or individual BP values
function [EST] = EstimationC_T(PATHS,REF,EST)

% error('Check equation A8 from Westerhof1991 to estimate C_T from Z_c')

switch EST.MethodC_T			
	case {1,2,3}	%Exponential fit methods
		addpath([PATHS.Root,'Others/P_out'])
		Pressure = REF.Pressure;
		T = EST.t_in(end);
		LVET = EST.LVET;
		switch EST.MethodC_T
			case 1
				[~,index_sys] = min(abs(EST.t_in - LVET));
				P_sys = Pressure(index_sys);
				tau = (T - LVET) / log((P_sys - EST.P_out) / (REF.DBP - EST.P_out));
			case 2
				t_fit = LVET;
				[tau,~,~,~] = ExponentialFit(Pressure,T,t_fit,EST.P_out,'P_out',EST.Plots);
			case 3
				t_fit = LVET + (T - LVET)/3;
				[tau,~,~,~] = ExponentialFit(Pressure,T,t_fit,EST.P_out,'P_out',EST.Plots);
		end
		EST.C_T = tau/(EST.R_T - EST.Z_0);
		
	case 4	%Area-method, Randall1976
		%"The portion of the waveform, from t1 (0.06-0.08 sec after the
		%peak of the dicrotic wave) to t2 (0.20-0.24 sec after t1) was
		%chosen to avoid aberrations known to occur in the brachial
		%waveform, when compared with the aortic waveform during systole
		%and the dicrotic wave"
		P = REF.Pressure;
		t = linspace(REF.t_in(1),REF.t_in(end),length(P));
		T = t(end);
		dt = t(2) - t(1);
		LVET = EST.LVET;
		
		[~,t_1] = min(abs(t - (2/3*LVET + 1/3*T)));	%ignore 1/3 of diastole
		[~,t_2]	= min(abs(t - 0.9*T));				%ignore last 10% of period
		
		P_12 = P(t_1:t_2);
		P_out = EST.P_out;
		P_1 = P_12(1);
		P_2 = P_12(end);
		
		EST.C_T = trapz(P_12 - P_out)*dt / (EST.R_T*(P_1 - P_2));
		
		if EST.Plots == 1
			figure, hold on, 
			plot(t,P,'k')
			area(t(t_1:t_2),P_12)
			legend('P','Area between P_1 and P_2')
		end
		
		%	Liu1986 (alternative implementation)
		%	Option 1: C = SV / [K (P_es - DBP)] where K = (A_s + A_d)/A_d, and
		%	A_s and A_d are the areas under the pressure curve
		%	during systole and diastole, respectively
		%	Option 2: C = A_d / [R (P_es - DBP)]
		
	case 5	%Two-area method, Self1994
		P = REF.Pressure;
		Q = REF.Q_in;
		t = REF.t_in;
		
		addpath([PATHS.Root,'Others/InterpolateSpline/'])
		[P,Q] = InterpolateSpline(t,P,Q);
		t = linspace(REF.t_in(1),REF.t_in(end),length(P));
		dt = t(2) - t(1);
		LVET = EST.LVET;
% 		[~,t_SBP] = max(P);
% 		t_SBP = t(t_SBP);
% 		t_mid = t_SBP;
		R_T = EST.R_T;

		t_0 = 1;
		[~,t_1] = min(abs(t - LVET));
		t_2 = length(P);
		
		Q_01 = Q(t_0:t_1);
		Q_12 = Q(t_1:t_2);
		P_01 = P(t_0:t_1);
		P_12 = P(t_1:t_2);
		P_out = EST.P_out;
		
		% int(Q) - 1/R_T int(P - P_out) = C_T (P_2 - P_1)		
		RHS_matrix	= [P(t_1) - P(t_0); P(t_2) - P(t_1)];
		LHS_matrix	= [trapz(Q_01)*dt - 1/R_T*trapz(P_01 - P_out)*dt; 
					trapz(Q_12)*dt - 1/R_T*trapz(P_12 - P_out)*dt];
		
		C_T = linsolve(RHS_matrix,LHS_matrix);
		
		EST.C_T = C_T;
		
	case {6,7,8}	%DBP, PP and SV/PP methods
		PP_ref = REF.SBP - REF.DBP;
		SV = REF.SV;
		EST.C_T = SV/PP_ref;
		
		switch EST.MethodC_T
			case 6
				EST	= DBP_C_T_method(PATHS,REF,EST);
			case 7
				EST = PP_C_T_method(PATHS,REF,EST);
		end		
		
	case 9	%Optimised Wk3 parameters (faster and more general)
		addpath('~/Haemodynamic_Tools/Version6/Others/ImpedanceAnalysis/')
		t = REF.t_in;	%[s]
		Q = REF.Q_in;	%[m3/s]
		P = REF.Pressure;	%[Pa]
		P_out = EST.P_out;	%[Pa]
		R_T	= EST.R_T;	%[Pa*s/m3]
		%	Initial C_T value using the SV/PP method
		EST_temp.MethodC_T = 8;
		EST_temp = EstimationC_T(PATHS,REF,EST_temp);
		C_T_temp = EST_temp.C_T;	%[m3/Pa]
		Z_0 = EST.Z_0;
		R_2 = R_T - Z_0;
		[~,~,C_T] = RCR_EstimationCBP_v2(P,Q,t,Z_0,R_2,C_T_temp,P_out,EST.Plots);
		EST.C_T	= C_T;	%[Pa*s/m3]		
		
	otherwise
		error('Choose a valid method for the estimation of C_T')
end

if EST.C_T <= 0
	error('Total compliance <= 0')
end

end


%%	Function definitions
% %	DBP method
function [EST] = DBP_C_T_method(PATHS,REF,EST)
Scale_P = 133.32;

EST = EstimationTau(EST);

if EST.Plots == 1
	figure
	hold on
	t = linspace(0,REF.t_in(end),length(REF.Pressure_CBP));
	plot(t,ones(length(REF.Pressure_CBP),1)*min(REF.Pressure_CBP)/Scale_P,'--k','LineWidth',2)
	t = linspace(0,REF.t_in(end),length(REF.Pressure));
	plot(t,REF.Pressure/Scale_P,'-k','LineWidth',2)
end

%	Initial quasi-steady state (QSS) estimation
EST = P_0_iteration(EST,0,'Wk3');
if EST.Plots == 1
	plot(EST.t_in,EST.Pressure/Scale_P,'--r','LineWidth',4)
end

% plot(EST.Pressure,'--')

DBP_error = EST.Pressure(end)/REF.DBP;
DBP_tol		= 1;		%1% error tolerance for the iterative calculation of C_T

% 		Reduce the error (DBP,estimation - DBP,reference) after	updating C_T
while abs(1 - DBP_error)*100 > DBP_tol	
	if DBP_error < sqrt(0.5)
		EST.C_T	= EST.C_T / 0.5;
	else
		EST.C_T	= EST.C_T / DBP_error^2;
	end
	EST = EstimationTau(EST);
	EST.P_0 = REF.DBP;
	
	EST = P_0_iteration(EST,0,'Wk3');
	
	%	Error function = DBP,estimation / DBP,reference
	DBP_error = EST.Pressure(end)/REF.DBP;
	
% 	if EST.Plots == 1
% 		plot(EST.t_in,EST.Pressure/Scale_P,'--','LineWidth',2)
% 	end
end


%%	Plots
if EST.Plots == 1
	plot(EST.t_in,EST.Pressure/Scale_P,'-.b','LineWidth',4)
	hold off
	
	T = EST.t_in(end);
	
	legend('DBP','pBP','Initial','Optimal')
	xlabel('Time [s]')
	ylabel('Pressure [mmHg]')
	xlim([0,T])
	ylim([70,125])
	set(gca,'XTick',[0:0.2:T+0.2])
	set(gca,'YTick',[80:20:120])
	box on
	legend boxoff
	set(gca,'FontSize',40)
	
	addpath([PATHS.Root,'Others/PlotSave/'])
	PlotSave(PATHS.Figures,'C_T_DBP_method')
end

end

%	Pulse pressure method
function [EST] = PP_C_T_method(PATHS,REF,EST)
EST = EstimationTau(EST);

if EST.Plots == 1
	figure
	hold on
	plot(REF.Pressure,'k')
end

PP_error	= 100;
PP_tol		= 1e-3;	% mmHg; error tolerance for the iterative calculation of C_T

%	1. Initial estimation
PP_ref = REF.PP;
C_T = EST.C_T;
SV = REF.SV;

jj = 1;
while abs(PP_error) > PP_tol
	%	2. Run simulation using C_T(jj) to calculate PP_est(jj)
	EST = EstimationTau(EST);
	EST = P_0_iteration(EST,0,'Wk3');
	
	%	3. Calculate PP_delta(jj) and C_T(jj+1) for the next iteration
	SBP_est(jj)	= max(EST.Pressure);
	DBP_est(jj) = min(EST.Pressure);
	PP_est(jj)	= SBP_est(jj) - DBP_est(jj);
	PP_delta(jj)= PP_ref - PP_est(jj);
	
	C_T_1(jj)	= C_T(jj) + (-SV/(PP_est(jj)^2)) * PP_delta(jj);	%original Taylor expansion derivation
	if C_T_1(jj) <= 0
		C_T_1(jj)	= C_T(jj) + (-SV/(PP_est(jj)^2)) * PP_delta(jj) / 2;	%divide increase by 2
	end
	if C_T_1(jj) <= 0 
		error('C_T <= 0 using PP method')
	end
	% 	C_T_2(jj)	= C_T(jj)*(2 - PP_ref/PP_est(jj));
	% 	C_T_3(jj)	= C_T(jj)/(1 + PP_delta(jj)/PP_est(jj));
	
	C_T(jj+1)	= C_T_1(jj);
	EST.C_T		= C_T(jj+1);
	
	%	4. Calculate current relative error
	PP_error = PP_delta(jj)/PP_ref;
	
	jj = jj + 1;

	if EST.Plots == 1
		plot(EST.Pressure,'--')
	end
end


% %%
% PP(1)	= PP_ref;
% PP_delta(1) = 0;
% C(1)	= SV / PP(1);
% EST.C_T = C(1);
% 
% jj = 1;
% while abs(PP_error) > PP_tol
% 	%	2. Run simulation using C_T(jj) to calculate PP_est(jj)
% 	EST = EstimationTau(EST);
% 	EST = P_0_iteration(EST,1);
% 	
% 	%	3. Calculate PP_delta(jj) and C_T(jj+1) for the next iteration
% 	SBP_est(jj)	= max(EST.Pressure);
% 	DBP_est(jj) = min(EST.Pressure);
% 	PP(jj+1)		= SBP_est(jj) - DBP_est(jj);
% 	PP_delta(jj+1)	= PP(jj+1) - PP_ref;
% 	
% 	C(jj+1)	= C(jj) + SV/(PP(jj)^2)*PP_delta(jj+1);	%original Taylor expansion derivation
% 	% 	C_T_2(jj)	= C_T(jj)*(2 - PP_ref/PP_est(jj));
% 	% 	C_T_3(jj)	= C_T(jj)/(1 + PP_delta(jj)/PP_est(jj));
% 	
% 	C_T(jj+1)	= C(jj);
% 	EST.C_T		= C_T(jj+1);
% 	
% 	%	4. Calculate current estimation error
% 	PP_error = PP_delta(jj+1);
% 	
% 	jj = jj + 1;
% 	
% 	plot(EST.Pressure,'--')
% end


%%	Plots
if EST.Plots == 1
	figure
	hold on
	plot(REF.Pressure,'k')
	plot(EST.Pressure,'--b')
	hold off
	
	title('C_{T}: pulse pressure (PP) method')
	xlabel('Time [s]')
	ylabel('Pressure [mmHg]')
	legend('Reference','Estimation')
	set(gca,'FontSize',24)
	
	SaveEstimationPlots([PATHS.Root,'Windkessel/Parameter_study/'],...
		['C_T_PP_method_',EST.Windkessel])

end

end