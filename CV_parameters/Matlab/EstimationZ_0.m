%%	Estimate Z_0 from reference data

function [EST] = EstimationZ_0(PATHS,REF,EST)
switch EST.MethodZ_0
	case 1	%10.1161/01.CIR.62.1.105
		R_T = REF.MBP / mean(REF.Q_in);%in the article, R_T was calculated as MBP/Q_mean
		EST.Z_0 = 1/10*R_T;
		
	case 2	%10.1161/01.CIR.62.1.105
		R_T = REF.MBP / mean(REF.Q_in);%in the article, R_T was calculated as MBP/Q_mean
		EST.Z_0 = 1/20*R_T;
		
	case 3
		EST.Z_0 = (REF.MBP - REF.DBP) / max(REF.Q_in);
		
	case 4
		rho = REF.rho;
		A_d = min(REF.Asc_Ao.A_in);		%[m^2]
		Z_0 = rho*EST.PWV/A_d;			%[Pa*s/m3]
		EST.Z_0 = Z_0;
		
	case {5,6,7,8,9,10,11,12,13}
		addpath([PATHS.Functions,'ImpedanceAnalysis/'])
		
		%	Ensure that P vector starts from DBP
		P_raw = REF.Pressure;
		while P_raw(2) - P_raw(1) <= 0
			P_raw = [P_raw(2:end);P_raw(1)];
			disp('P_{raw} must start from DBP')
		end
		
		switch EST.MethodZ_0
			case 5	%2-12th harmonics; Nichols1977 - 10.1016/S0021-9290(99)00180-3
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,2,12,EST.Plots);
			case 6	%6-10th harmonics; Segers2000 - 10.1016/S0021-9290(99)00180-3
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,6,10,EST.Plots);
			case 7	%1-8th harmonics (Clarke et al (1978))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,1,8,EST.Plots);
			case 8	%1-9th harmonics (Peluso et al (1978))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,1,9,EST.Plots);
			case 9	%2-10th harmonics (Mitchell et al (1994))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,2,10,EST.Plots);
			case 10	%3-10th harmonics (Hughes and Parker (2009))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,3,10,EST.Plots);
			case 11	%4-10th harmonics (Tabima et al (2012))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,4,10,EST.Plots);
			case 12	%6-8th harmonics (Abel (1971))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,6,8,EST.Plots);
			case 13	%4-8th harmonics (Qureshi et al (2018))
				EST.Z_0 = Z_0_FromHarmonics(PATHS,REF.t_in,REF.Q_in,P_raw,EST.P_out,4,8,EST.Plots);
		end
		
	case {14 15 16 17}
		%	Ensure that P vector starts from DBP
		P_raw = REF.Pressure;
		while P_raw(2) - P_raw(1) <= 0
			P_raw = [P_raw(2:end);P_raw(1)];
			disp('P_{raw} must start from DBP')
		end
		
		%	Same vector length required
		addpath([PATHS.Functions,'InterpolateSpline/'])
		[P_ref, Q_ref] = InterpolateSpline(REF.t_in,P_raw,REF.Q_in);

		%	Early systole (_es) defined as [t = 0, t = Q_peak]
		[~,Q_loc_max] = max(Q_ref);
		loc_max = Q_loc_max;

		%	Early systole (_es) defined as [t = 0, t = P_peak]
% 		[~,P_loc_max] = max(Q_ref);
% 		loc_max = P_loc_max;

		Q_es = Q_ref(1:loc_max);
		NonZeroQ = Q_es ~= 0;
		Q_es = Q_es(NonZeroQ);	%remove 0 Q-values
		P_es = P_ref(NonZeroQ);
		
		if EST.Plots == 1
			figure
			hold on
			plot(Q_ref,P_ref,'ko')
			plot(Q_es,P_es,'b--')
			xlabel('Flow [m^3/s]'), ylabel('Pressure [Pa]')
		end
		
		switch EST.MethodZ_0
			case 14	%mean of P_es to Q_es ratio
				Z = ((P_es - P_es(1))./(Q_es - Q_es(1)));
				EST.Z_0 = mean(Z(~isnan(Z)));
								
				if EST.Plots == 1
					x = [0, Q_es(end)];
					y = [P_es(1), EST.Z_0*Q_es(end) + P_es(1)];
					plot(x,y,'r-.')
					title('Z_0 (time-domain): mean [(P - DBP) / Q]')
					legend('Cycle','Early systole','Estimated Z_0','Location','southeast')
					set(gca,'FontSize',24)
					
					SaveEstimationPlots([PATHS.Root,'Windkessel/Parameter_study/',...
						'Z_0_time_mean_PQ'])
				end
				
			case 15	%least-squares fit for early-systole
				coeffs	= polyfit(Q_es,P_es,1);
				EST.Z_0 = coeffs(1);
				
				if EST.Plots == 1
					fittedX = linspace(min(Q_es), max(Q_es), length(Q_es));
					fittedY = polyval(coeffs, fittedX);
					plot(fittedX,fittedY,'r-.')
					title('Z_0 (time-domain): least-square fit')
					legend('Cycle','Early systole','Least-square fit','Location','southeast')
					set(gca,'FontSize',24)
					
					SaveEstimationPlots([PATHS.Root,'Windkessel/Parameter_study/',...
						'Z_0_time_least_square'])
				end
				
			case 16	%P_es to Q_es ratio at time of maximum dQ
				
				[~,dQ_loc_max] = max(diff(Q_es));
				% 	[~,dP_loc_max] = max(diff(P_es));
				
				loc_max = dQ_loc_max;
				
				%	Ensure that the max is not found too early
				kk = 0;
				while loc_max == 1 && kk < 5
					kk = kk+1;
					[~,dQ_loc_max] = max(diff(Q_es(kk+1:end)));
				% 	[~,dP_loc_max] = max(diff(P_es(kk+1:end)));
							
					loc_max = dQ_loc_max;
				end
				loc_max = kk + loc_max;
				
				Q = Q_es(loc_max);
				P = P_es(loc_max);
				EST.Z_0 = (P - P_es(1))/(Q - Q_es(1));
				
				if EST.Plots == 1
					x = [0, Q_es(end)];
					y = [P_es(1), EST.Z_0*Q_es(end) + P_es(1)];
					plot(x,y,'r-.')
					title('Z_0 (time-domain): P / Q at max dQ')
					legend('Cycle','Early systole','Estimated Z_0','Location','southeast')
					set(gca,'FontSize',24)
					
					SaveEstimationPlots([PATHS.Root,'Windkessel/Parameter_study/',...
						'Z_0_time_mean_PQ_max_dP'])
				end
				
			case 17	%P_es to Q_es ratio [t = 0, t = maximum dQ]
				[~,dQ_loc_max] = max(diff(Q_es));
				% 				[~,dP_loc_max] = max(diff(P_es));
				
				loc_max = dQ_loc_max;
				
				%	Ensure that the max is not found too early
				kk = 0;
				while loc_max == 1 && kk < 5
					kk = kk+1;
					[~,dQ_loc_max] = max(diff(Q_es(kk+1:end)));
				% 	[~,dP_loc_max] = max(diff(P_es(kk+1:end)));
							
					loc_max = dQ_loc_max;
				end
				loc_max = kk + loc_max;
				
				Q = Q_es(1:loc_max);
				P = P_es(1:loc_max);
				Z = ((P - P(1))./(Q - Q(1)));
				EST.Z_0 = mean(Z(~isnan(Z)));
				
				if EST.Plots == 1
					x = [Q_es(1), Q_es(end)];
					y = [P_es(1), EST.Z_0*Q_es(end) + P_es(1)];
					plot(x,y,'r-.')
					title('Z_0 (time-domain): mean (P / Q) until max dQ')
					legend('Cycle','Early systole','Estimated Z_0','Location','southeast')
					set(gca,'FontSize',24)
					
					SaveEstimationPlots([PATHS.Root,'Windkessel/Parameter_study/',...
						'Z_0_time_mean_PQ_0_to_max_dP'])
				end
				
		end
		
	case 18	%Optimised Wk3 parameters (faster and more general)
		addpath([PATHS.Functions,'ImpedanceAnalysis/'])
		t = REF.t_in;	%[s]
		Q = REF.Q_in;	%[m3/s]
		P = REF.Pressure;	%[Pa]
		P_out = EST.P_out;			%[Pa]
		R_T	= EST.R_T;				%[Pa*s/m3]
		%	Initial C_T value using the SV/PP method
		EST_temp.MethodC_T = 8;
		EST_temp = EstimationC_T(PATHS,REF,EST_temp);	%[m3/Pa]
		C_T_temp = EST_temp.C_T;		
		Z_0_temp = R_T*0.05;
		R_2 = R_T - Z_0_temp;
		[Z_0] = RCR_EstimationCBP_v2(PATHS,P,Q,t,Z_0_temp,R_2,C_T_temp,P_out,EST.Plots);
		EST.Z_0	= Z_0;	%[Pa*s/m3]
		
	otherwise
		error('Choose a valid method for the estimation of Z_0')
end

if EST.Z_0 > EST.R_T
	error('Z_0 > R_T')
elseif EST.Z_0 <= 0
	error('Z_0 <= 0')
elseif isempty(EST.Z_0)
	error('Z_0 is empty')
end

end

