%%	Estimate R_T from whole waveform or individual BP values

function [EST] = EstimationR_T(REF,EST)

switch EST.MethodR_T
	case 1	%Scenario 1
		EST.R_T = (REF.MBP - EST.P_out) / mean(REF.Q_in);
	case 2	%Scenario 2
		EST.R_T = (REF.MBP - EST.P_out) / mean(REF.Q_in);
		
	otherwise
		error('Choose a valid method for the estimation of R_T')
end

if EST.R_T <= 0
	error('Total resistance <= 0')
end

end