%%	Estimate tau from R_T and C_T

function [EST] = EstimationTau(EST)
	EST.tau		= (EST.R_T - EST.Z_0)*EST.C_T;
end