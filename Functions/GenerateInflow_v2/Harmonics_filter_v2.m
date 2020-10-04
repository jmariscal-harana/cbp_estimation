%:::::::::::::::::::::::::::::::::%
% FILTER WAVEFORM USING HARMONICS %
%:::::::::::::::::::::::::::::::::%
%
% Jordi Alastruey
% King's College London
% March 2014

%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (23/01/19)
%
%==========================================================================

function [tfunc,funcf,T] = Harmonics_filter_v2(tfunc,func,NHarm,SR,fmt,OutputFile)

% Cardiac period in s
T = tfunc(end) - tfunc(1);

if nargin == 6
	[FT,H_amp,H_ph] = Harmonics(func(1:end),NHarm,fmt,OutputFile,T);
else
    [FT,H_amp,H_ph] = Harmonics(func(1:end),NHarm,fmt);
end

tfunc = 0:1/SR:T;
funcf = FT(1);
for j=1:NHarm
    funcf = funcf + H_amp(j)*sin(2*j*pi*tfunc/T + H_ph(j));
end

T = tfunc(end) - tfunc(1);

end