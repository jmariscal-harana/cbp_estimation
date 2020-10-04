%	Optimised Wk3 calculation to get Z_0
clear
clc

load Sujeto_1D
p=REF(1).P_central;
q=REF(1).Q_central;
t=REF(1).t_central;

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

addpath('~/Haemodynamic_Tools/Version6/Others//Wk3_opt_Aramburu/')
[Z_0,R_2,C_T]		= Optimise_Wk3(t,q,p,pout,R2,C,1,1);	%[mm*s/mL]
