function [R_1, R_2, C_T, P_out] = RCR_EstimationCBP_v2(P_in,Q_in,t,R_1,R_2,C_T,P_out,Plots)
P_in = P_in(:)';	
Q_in = Q_in(:)';

%	Same vector length required
addpath('~/Haemodynamic_Tools/Version6/Others/InterpolateSpline/')
[P_in,Q_in,t] = InterpolateSpline(t,P_in,Q_in);

P_in = P_in(:)';	
Q_in = Q_in(:)';
t = t(:)';
MP = length(P_in);

P_in(end) = [];
Q_in(end) = [];
t(end) = [];

dt=t(2)-t(1);       % resolución temporal (paso en el tiempo)
M=MP-1;   % nº de muestras tomadas (nº natural)
t=[0:dt:(M-1)*dt];
fm=1/dt;       % frecuencia de muestreo
T=M*dt;       % periodo de observación (ciclo natural del latido)
fo=1/T;       % resolución espectral (paso en la frecuencia)
f=0:fo:(M-1)*fo;  % vector de frecuencias de la fft (solo la mitad aprox. son positivas)

Q_fft=fft(Q_in);      % componentes de frecuencia de la corriente

r1=mean(P_in); % valor medio de p
r2=min(P_in); % valor mínimo de p
r3=max(P_in); % valor máximo de p

w1=1;
w2=1;
w3=1;

s=0.01; % sensibilidad
iter=10000; % nº de iteraciones
itern=1;

% Obtención de psim_tot y del error inicial
Z=R_1+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2);
Psim=Q_fft(1:floor(M/2+1)).*Z;
Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
psim=real(ifft(Psim));
psim_tot=psim+P_out;
% ferr=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;	%error function looking at SBP, MBP, DBP errors
ferr(itern)=mean(sqrt((psim_tot-P_in).^2));		%error function looking at RMS errors

while (iter>0)
	% Obtención de psim_tot y del error con R_1+
	Z=R_1*(1+s)+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferrR_1p=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrR_1p=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con R_1-
	Z=R_1*(1-s)+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferrR_1m=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrR_1m=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con R_2+
	Z=R_1+R_2*(1+s)./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2*(1+s));
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferrR_2p=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrR_2p=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con R_2-
	Z=R_1+R_2*(1-s)./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2*(1-s));
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferrR_2m=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrR_2m=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con C+
	Z=R_1+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*(1+s)*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferrCp=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrCp=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con C-
	Z=R_1+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*(1-s)*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferrCm=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrCm=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con pout+
	Z=R_1+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out*(1+s);
% 	ferrpoutp=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrpoutp=mean(sqrt((psim_tot-P_in).^2));
	
	% Obtención de psim_tot y del error con pout-
	Z=R_1+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out*(1-s);
% 	ferrpoutm=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferrpoutm=mean(sqrt((psim_tot-P_in).^2));
	
	% Determinación del mínimo error
	[merr,mmerr]=min([ferr(itern) ferrR_1p ferrR_1m ferrR_2p ferrR_2m ferrCp ferrCm ferrpoutp ferrpoutm]);
	
	if (mmerr==1)
		break
	elseif (mmerr==2)
		R_1=R_1*(1+s);
	elseif (mmerr==3)
		R_1=R_1*(1-s);
	elseif (mmerr==4)
		R_2=R_2*(1+s);
	elseif (mmerr==5)
		R_2=R_2*(1-s);
	elseif (mmerr==6)
		C_T=C_T*(1+s);
	elseif (mmerr==7)
		C_T=C_T*(1-s);
	elseif (mmerr==8)
		P_out=P_out*(1+s);
	else
		P_out=P_out*(1-s);
	end
	
	% Obtención de psim_tot y del error final
	Z=R_1+R_2./(1+1i*2*pi*f(1:floor(M/2+1))*C_T*R_2);
	Psim=Q_fft(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+P_out;
% 	ferr=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	ferr(itern+1)=mean(sqrt((psim_tot-P_in).^2));
	
	%	Stop iterations if the relative difference between errors is smaller than 0.1%
	if abs(ferr(itern+1) - ferr(itern))/ferr(itern)*100 < 0.1
		itern=itern+1;
		break
	end
	
	% Decremantamos el nº de iteraciones
	iter=iter-1;
	itern=itern+1;
end

%	Condition: relative MBP error < 5%
if (mean(psim_tot) - mean(P_in)) / mean(P_in) * 100 >= 5
	warning('Relative MBP error >= %f',(mean(psim_tot) - mean(P_in)) / mean(P_in) * 100)
end

if Plots == 1
	figure
	plot(t,P_in,'b',t,psim_tot,'r'); grid on
	disp(['R_1 = ' num2str(R_1)])
	disp(['R_2 = ' num2str(R_2)])
	disp(['C_T = ' num2str(C_T)])
	disp(['P_{out} = ' num2str(P_out)])
	disp(['Iterations = ' num2str(itern)])
	disp(['Relative MBP error (EST - REF) [%] = ', num2str((mean(psim_tot) - mean(P_in)) / mean(P_in) * 100)])
end
end

