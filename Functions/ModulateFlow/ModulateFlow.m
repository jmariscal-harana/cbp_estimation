%% Remove noise from VTI/flow waveforms during diastolic downslope
%
% Samuel Vennin
% King's College London
% 11 October 2018
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (31/03/19) - SG-filtering OR spline interpolation; changed
%		conditions for anlaysis of downslope
%
%==========================================================================


%%
function [U_final] = ModulateFlow(U_start,Plots)
l = length(U_start);

if l < 100
% 	FRAMELEN = round(l/10)*2 + 1;
% 	ORDER = round(FRAMELEN/4);
% 	U_filter = sgolayfilt(U_start, ORDER, FRAMELEN);	%adaptative Savitzky-Golay filtering
	U_filter = spline(1:l,U_start,1:(l-1)/99:l);		%spline interpolation with 100 points
else
	U_filter = U_start;
end

l = length(U_filter);
[Umax ind_max] = max(U_filter);

%-Get systolic downslope right (a.k.a delete speckle contribution to US scan)
U_temp = U_filter(ind_max:end);
Noise  = find(U_temp < 0.2*max(U_temp));	%find U < 0.1*U_max to bypass noise
Noise2  = find(U_temp > 0.5*max(U_temp));	%find U > 0.1*U_max to exclude from downslope projection

U_temp2 = U_temp(Noise2(end):Noise(1));	
dU_temp = diff(U_temp2);
[slope ind_min] = min(dU_temp);	%find largest negative slope
ind_min = ind_min + Noise2(end);
% temp = U_filter(ind_max+ind_min) + slope*[1:1000]; if sum(temp < 0) == 0, error('x-intercept required'), end
temp = U_temp(ind_min) + slope*[1:1000]; if sum(temp < 0) == 0, error('x-intercept required'), end
temp = temp(temp>=0);
U_temp = [U_filter(1:ind_max+ind_min-1) temp];
U_final = [U_temp zeros(1,l-length(U_temp))];

if Plots ~= 0
	figure, hold on
	plot(U_start,'-k')
	plot(U_filter,'--b')
	plot(U_final,'--r')
	ylabel('VTI [m/s]')
	xlabel('Datapoints')
	legend('Input flow','Cubic spline interpolation','Downslope projection')
	set(gca,'fontsize',20)
end

% %- Get systolic upstroke right
% U_temp = U_final(1:ind_max);
% dU_temp = diff(U_temp);
% [slope ind_max_diff] = max(dU_temp);
% U_temp = U_final(ind_max_diff)/(ind_max_diff-1)*[0:ind_max_diff-1];
% U_final = [U_temp U_final(ind_max_diff+1:end)];

end