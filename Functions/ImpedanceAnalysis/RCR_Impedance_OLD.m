function [R_1, R_2, C_T] = RCR_Impedance(P_in,Q_in,P_out,T,Plots)
%%	INPUT
%	Input signals
P = P_in(1:end-1) - P_out;	%remove last sample shared with next cycle
Q = Q_in(1:end-1);			%remove last sample shared with next cycle

%	Same vector length required
if length(P) > length(Q)
	[Q] = SameVectorLength(P,Q);
elseif length(P) < length(Q)
	[P] = SameVectorLength(P,Q);
end


%%	Calculate impedance from P and Q using FFT
P_fft = fft(P);
Q_fft = fft(Q);
Z_fft = P_fft./Q_fft;

M = length(Z_fft);

R_1 = real(Z_fft(floor(M/2+1)));
R_2 = real(Z_fft(1))-R_1;
C_T = (-R_2-sqrt(R_2^2-4*imag(Z_fft(2))^2))/2/imag(Z_fft(2))/2/pi/R_2*T;

if Plots == 1
	INPUT.P_out = P_out;
	INPUT.P_0	= P(1);
	INPUT.R_T	= R_1 + R_2;
	INPUT.Z_0	= R_1;
	INPUT.C_T	= C_T;
	INPUT.Q_in	= Q_in;
	INPUT.t_in	= [0:1/round((length(Q_in)-1)/T):T]';
	
	addpath('~/Haemodynamic_Tools/Version6/Others/P_0_iteration/')	
	
	[INPUT] = P_0_iteration(INPUT,0,'Wk3');
		
	figure,	hold on
	plot(INPUT.t_in,P_in,'k')
	plot(INPUT.t_in,P_Wk3,'--b')	
end


end


%%	Function definitions
%	Length of 'ShortVector' becomes that of the other vector via spline interpolation
function [ShortVector] = SameVectorLength(Vector_1,Vector_2)

if length(Vector_1) < length(Vector_2)
	N_old = length(Vector_1);
	N_new = length(Vector_2);
	ShortVector = spline(linspace(0,1,N_old),Vector_1,linspace(0,1,N_new))';
elseif length(Vector_1) > length(Vector_2)
	N_old = length(Vector_2);
	N_new = length(Vector_1);
	ShortVector = spline(linspace(0,1,N_old),Vector_2,linspace(0,1,N_new))';
end

end