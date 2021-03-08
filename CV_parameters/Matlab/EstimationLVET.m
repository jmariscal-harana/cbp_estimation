%%	Estimate LVET from pressure or flow waveform

function [EST] = EstimationLVET(PATHS,REF,EST)

switch EST.MethodLVET
	case 1
		addpath([PATHS.Functions,'LVET/'])
		[ind_dic,ind_dia] = valve_closure(PATHS,REF.Pressure',REF.t_in(end),EST.Plots);
 		LVET_index = ind_dic;
		LVET_index = ind_dia;
		temp = LVET_index/length(REF.Pressure);
		LVET_index = round(temp*length(REF.t_in));
		EST.LVET = REF.t_in(LVET_index);
		
	case 2	%10.1088/1361-6579/aabe6a
		addpath([PATHS.Functions,'PulseAnalyse/'])
		Wave.v = REF.Pressure;
		Wave.fs = round(1/(REF.t_in(2) - REF.t_in(1)));	%1/dt
		Options.do_plot = false;
		[~,fid_pts] = PulseAnalyse10(Wave,Options);
		LVET_index = fid_pts.dic;	%dicrotic notch: < LVET (i.e. valve closure) for P_{Asc Ao}
		LVET_index = fid_pts.dia;	%diastolic peak: > LVET for P_{Asc Ao}
		if isinteger(LVET_index) && LVET_index > 0
			EST.LVET = REF.t_in(LVET_index);
		else
			EST.MethodLVET = 4;
			EST = EstimationLVET(PATHS,REF,EST);
			EST.MethodLVET = 2;
		end
		
	case 3	%10.1186/s12938-017-0341-z
		P = REF.Pressure;
		t = REF.t_in;
		%	Require length(t) = length(P) and length(t_dP) = length (dP)
		t = linspace(t(1),t(end),length(P))';
		t_dP = linspace(t(1),t(end),length(P)-1)';

		HR = 60/t_dP(end);
		WF = (0.5 - abs(0.5 - HR*t_dP/60)).^2;
		dP = diff(P);
		[~,LVET_index] = min(dP.*WF);
		EST.LVET = t_dP(LVET_index);
		if EST.Plots == 1
			Plot_WF(t,t_dP,P,dP,WF,LVET_index)
		end
		
	case 4
		A_men = 0.37;
		A_women = 0.40;
% 		EST.LVET = mean([A_men,A_women])*sqrt(REF.t_in(end));	%empirical eq: 10.1111/j.1542-474X.1997.tb00325.x
		EST.LVET = A_men*sqrt(REF.t_in(end));	%empirical eq: 10.1111/j.1542-474X.1997.tb00325.x
		
	case 5
		addpath([PATHS.Functions,'LVET/'])
		EST.LVET = LVET_from_Q(PATHS,REF.Q_in,REF.t_in,EST.Plots);
		
	otherwise
		error('Choose a valid method for the estimation of LVET')
end

end


%%	Function definitions
function Plot_WF(t,t_dP,P,dP,WF,LVET_index)

figure
hold on
plot(t,(P-min(P))/(max(P)-min(P)),'k','LineWidth',2)
plot(t_dP,(dP-min(dP))/(max(dP)-min(dP)),'--b')
plot(t_dP,(dP.*WF-min(dP.*WF))/(max(dP.*WF)-min(dP.*WF)),'-.r')
plot(t(LVET_index),(P(LVET_index)-min(P))/(max(P)-min(P)),'ko','MarkerSize',20)
hold off

title('t_{sys} (pressure): weight function (WF)')
xlabel('Time [s]')
ylabel('Normalised magnitude')
legend('Pressure','dP/dt','(dP/dt)WF','Location','southeast')
set(gca,'FontSize',24)

% SaveEstimationPlots('~/Haemodynamic_Tools/Version6/Windkessel/Parameter_study/LVET_from_P_WF')

end

