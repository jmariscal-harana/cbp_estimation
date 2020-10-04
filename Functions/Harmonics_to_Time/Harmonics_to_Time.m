function [tfunc, funcf] = Harmonics_to_Time(HarmonicsFile)
%	Read .BCS file
Harmonics = textread(HarmonicsFile);

%	Extract number of harmonics, period, and DC component
NHarm = Harmonics(1,1);
T = Harmonics(1,2);
FT = Harmonics(1,3);

%	Extract amplitude and phase for each harmonic	
H_amp = Harmonics(2:end,1);
H_ph = Harmonics(2:end,2);

%	Calculate sampling rate
if mod(NHarm,2) == 0	%odd NHarm
	SR = 2*NHarm/T;
else
	SR = 2*(NHarm-1)/T;
end

%	Calculate time and function values
tfunc = 0:1/SR:T;
funcf = FT(1);
for j=1:NHarm
    funcf = funcf + H_amp(j)*sin(2*j*pi*tfunc/T + H_ph(j));
end

end