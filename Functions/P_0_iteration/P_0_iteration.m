%%	Calculate P_0 ensuring quasi-steady state (QSS), which implies P_Wk(end) = P_Wk(1)
%	Input
%	-DATA: reference parameters required to run a Wk simulation
%	-Wk_cycles: compute current cycle (0) or all cycles (1)
%
%	Output
%	-DATA.P_0: updated P_0 QSS value
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (09/05/18) - based on Parameter_study.m
%	v1.1 (11/09/19) - linspace(Q) for t_in
%
%==========================================================================

function [DATA] = P_0_iteration(DATA,Wk_cycles,Windkessel)

t_original	= DATA.t_in;
Q_original	= DATA.Q_in;

t_in	= linspace(t_original(1),t_original(end),length(Q_original))'; 
Q_in	= Q_original(:);

j = 0;
P_0_error = 100;

switch Wk_cycles
	case 0
		%	Method 1: compute current cycle only (faster)		
		DATA.t_in	= t_in;
		DATA.Q_in	= Q_in;
		
		while abs(P_0_error) > 0.1	% P (end) - P_0 (current cycle)
						
			switch Windkessel
				case 'Wk2'
					[Pressure] = Wk2_PressureFunction(DATA.P_out,DATA.P_0,DATA.R_T,DATA.C_T,DATA.Q_in,DATA.t_in);
				case 'Wk3'
					[Pressure] = Wk3_PressureFunction(DATA.P_out,DATA.P_0,DATA.R_T,DATA.Z_0,DATA.C_T,DATA.Q_in,DATA.t_in);
			end
						
			j = j + 1;	%iteration counter
			P_0_error	= (Pressure(end) - Pressure(1))/Pressure(1)*100;	%percentage error
			
			if j > 99
				error('100 iterations and still no convergence! Check input parameters for subject: %s',DATA.ID)
			end
			
			DATA.P_0	= Pressure(end);	% P_0 extracted as the end-point of the current cardiac cycle
			
		end	

	case 1
		%	Method 2: keep previous cycles
		DATA.t_in	= t_in(1);
		DATA.Q_in	= Q_in(1);
		
		while abs(P_0_error) > 0.1	% P (end) - P_0 (current cycle)
			
			DATA.t_in	= cat(1, DATA.t_in, t_in(2:end)+j*t_in(end));
			DATA.Q_in	= cat(1, DATA.Q_in, Q_in(2:end));
			
			switch Windkessel
				case 'Wk2'
					[P_Wk_cycles] = Wk2_PressureFunction(DATA.P_out,DATA.P_0,DATA.R_T,DATA.C_T,DATA.Q_in,DATA.t_in);
				case 'Wk3'
					[P_Wk_cycles] = Wk3_PressureFunction(DATA.P_out,DATA.P_0,DATA.R_T,DATA.Z_0,DATA.C_T,DATA.Q_in,DATA.t_in);
			end
			
			j = j + 1;	%iteration counter
			P_0_error = (P_Wk_cycles(end) - P_Wk_cycles(end - length(t_in) + 1))/P_Wk_cycles(end - length(t_in) + 1)*100;	%percentage error
			
			if j > 99
				error('100 iterations and still no convergence! Check input parameters for subject: %s',DATA.ID)
			end
			
		end
		
		DATA.P_0 = P_Wk_cycles(end);	% P_0 extracted as the end-point of the last cardiac cycle
		
end

%	Recover the original time and flow vectors
DATA.t_in		= t_original;
DATA.Q_in		= Q_original;

switch Windkessel
	case 'Wk2'
		DATA.Pressure	= Wk2_PressureFunction(DATA.P_out,DATA.P_0,DATA.R_T,DATA.C_T,DATA.Q_in,DATA.t_in);
	case 'Wk3'
		DATA.Pressure	= Wk3_PressureFunction(DATA.P_out,DATA.P_0,DATA.R_T,DATA.Z_0,DATA.C_T,DATA.Q_in,DATA.t_in);
end

end

