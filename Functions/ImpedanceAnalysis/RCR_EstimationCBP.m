close all
clear
clc

load Sujeto_1D
p=REF(1).P_central';
q=REF(1).Q_central';
t=REF(1).t_central';

M=length(t);   % nº de muestras tomadas (nº natural)
tm=1e-3;       % resolución temporal (paso en el tiempo)
fm=1/tm;       % frecuencia de muestreo
to=M*tm;       % periodo de observación (ciclo natural del latido)
fo=1/to;       % resolución espectral (paso en la frecuencia)
f=0:fo:(M-1)*fo;  % vector de frecuencias de la fft (solo la mitad aprox. son positivas)

Q=fft(q);      % componentes de frecuencia de la corriente

r1=mean(p); % valor medio de p
r2=min(p); % valor mínimo de p
r3=max(p); % valor máximo de p

w1=1;
w2=1;
w3=1;

s=0.001; % sensibilidad

R1=0.05; %abs(randn); % valor inicial de R1
R2=1.0; %abs(randn); % valor inicial de R2
C=2.0; %abs(randn); % valor inicial de C
pout=30; %abs(randn); % valor inicial de pout

% Obtención de psim_tot y del error inicial
Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2);
Psim=Q(1:floor(M/2+1)).*Z;
Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
psim=real(ifft(Psim));
psim_tot=psim+pout;
ferr=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;

iter=1000; % nº de iteraciones
itern=1;

while (iter>0)
	% Obtención de psim_tot y del error con R1+
	Z=R1*(1+s)+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferrR1p=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con R1-
	Z=R1*(1-s)+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferrR1m=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con R2+
	Z=R1+R2*(1+s)./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2*(1+s));
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferrR2p=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con R2-
	Z=R1+R2*(1-s)./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2*(1-s));
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferrR2m=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con C+
	Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*(1+s)*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferrCp=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con C-
	Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*(1-s)*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferrCm=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con pout+
	Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout*(1+s);
	ferrpoutp=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Obtención de psim_tot y del error con pout-
	Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout*(1-s);
	ferrpoutm=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Determinación del mínimo error
	[merr,mmerr]=min([ferr ferrR1p ferrR1m ferrR2p ferrR2m ferrCp ferrCm ferrpoutp ferrpoutm]);
	
	if (mmerr==1)
		break
	elseif (mmerr==2)
		R1=R1*(1+s);
	elseif (mmerr==3)
		R1=R1*(1-s);
	elseif (mmerr==4)
		R2=R2*(1+s);
	elseif (mmerr==5)
		R2=R2*(1-s);
	elseif (mmerr==6)
		C=C*(1+s);
	elseif (mmerr==7)
		C=C*(1-s);
	elseif (mmerr==8)
		pout=pout*(1+s);
	else
		pout=pout*(1-s);
	end
	
	% Obtención de psim_tot y del error final
	Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2);
	Psim=Q(1:floor(M/2+1)).*Z;
	Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];
	psim=real(ifft(Psim));
	psim_tot=psim+pout;
	ferr=w1*(mean(psim_tot)-r1)^2+w2*(min(psim_tot)-r2)^2+w3*(max(psim_tot)-r3)^2;
	
	% Decremantamos el nº de iteraciones
	iter=iter-1;
	itern=itern+1;
end
plot(t,p,'b',t,psim_tot,'r'); grid on
disp(['R1= ' num2str(R1)])
disp(['R2= ' num2str(R2)])
disp(['C= ' num2str(C)])
disp(['pout= ' num2str(pout)])
disp(['nº de iter.= ' num2str(itern)])