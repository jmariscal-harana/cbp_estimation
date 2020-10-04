function [Z_0] = Z_0_FromHarmonics(PATHS,t_in,Q_in,P_in,P_out,H_1,H_2,Plots,Z_0_ref)
%%	INPUT
%	Harmonic range for the calculation of Z_0
Harmonic_range		= [H_1:H_2];			%chosen harmonic range for Z_0 estimation 
Harmonic_range_id	= Harmonic_range + 1;	%exclude DC component

%	Input signals
%	Same vector length required
addpath([PATHS.Root,'Others/InterpolateSpline/'])
[Q_in, P_in] = InterpolateSpline(t_in,Q_in,P_in);

P = P_in(1:end-1) - P_out;	%remove last sample shared with next cycle
Q = Q_in(1:end-1);			%remove last sample shared with next cycle


%%	Calculate impedance from P and Q using FFT
P_fft = fft(P);
Q_fft = fft(Q);

P_fft = P_fft(1:Harmonic_range_id(end));	%harmonic range only
Q_fft = Q_fft(1:Harmonic_range_id(end));	%harmonic range only

Z_fft = P_fft./Q_fft;
Z_abs = abs(Z_fft);

x_Z			= 0:Harmonic_range(end);
x_Z_id		= 1:Harmonic_range_id(end);
x_Z_0		= Harmonic_range;
x_Z_0_id	= Harmonic_range_id;

Z_0 = mean(Z_abs(x_Z_0_id));

if Plots == 1
	figure
	hold on
	if nargin == 9
		plot([1 x_Z_id(end)]-1,[Z_0_ref Z_0_ref],'--r')
	end
	plot(x_Z,Z_abs(x_Z_id),'-ok','MarkerFaceColor','k')
	plot(x_Z_0,Z_abs(x_Z_0_id),'-ob','MarkerFaceColor','b')
	hold off
	if nargin == 9
		legend('Reference R_1','Impedance','Impedance (Z_0 calculation)')
	else
		legend('Impedance','Impedance (Z_0 calculation)')
	end
	axis([-1 Harmonic_range(end) 0 1.1*max(Z_abs)])
	set(gca,'XTick',[x_Z])
	set(gca,'fontsize',18)
	xlabel('Harmonic number')
	ylabel('Modulus Z [mmHg s/mL]')
	box off
end

end

