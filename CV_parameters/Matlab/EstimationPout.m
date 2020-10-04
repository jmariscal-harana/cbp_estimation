%%	Estimate P_out from whole waveform or individual BP values
function [EST] = EstimationPout(PATHS,REF,EST)

switch EST.MethodP_out		
	case 1
		EST.P_out = 0;
		
	case 2
		EST.P_out = 0.5*REF.DBP;
		
	case 3
		EST.P_out = 0.7*REF.DBP;
		
	case {4,5}	%Decay time methods
		addpath([PATHS.Root,'Others/P_out'])
		Pressure = REF.Pressure;
		T = EST.t_in(end);
		LVET = EST.LVET;
		
		switch EST.MethodP_out
			case 4
				t_fit = LVET;
			case 5
				t_fit = LVET + (T - LVET)/3;
		end
		
		switch EST.MethodC_T
			case 0
				EST.tau = REF.tau;
				[~,EST.P_out] = ExponentialFit(Pressure,T,t_fit,EST.tau,'tau',EST.Plots);
			otherwise
				[EST.tau,EST.P_out] = ExponentialFit(Pressure,T,t_fit,0,'none',EST.Plots);
		end
		
	otherwise
		error('Choose a valid method for the estimation of P_out')
end

if EST.P_out > 0.75*REF.DBP
	EST.P_out = REF.DBP * 1/2;
	warning('Estimated P_out > 3/4 * DBP; P_out set to DBP/2')
end

end