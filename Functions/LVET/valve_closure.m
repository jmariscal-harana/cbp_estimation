%% valve_closure
%
% Samuel Vennin
% King's College London
% 11 February 2016

%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v1.0 2018
%
%==========================================================================

function [ind_dic,ind_dia] = valve_closure(PATHS,P,T,Plots)
Scale_P = 133.32;

time = linspace(0,T,length(P));
dP = diff(P);
[~, imax] = max(P);
DBP = min(P);
P_len = length(P);
P_temp = P(imax+1:round(0.8*P_len));	%pressure from global maximum to 80% of T
dP_temp = diff(P_temp);
%    ddP_temp = diff(dP_temp);
i_intersec = intersection(dP_temp);	%find points before dP 0-crossings

% i_intersec = [];

if Plots == 1
	figure,	hold  on
	plot(time,P/Scale_P,'-k','LineWidth',2)
	plot(time([imax+1:round(0.8*P_len)]),P_temp/Scale_P,'--b','LineWidth',4)
	if isempty(i_intersec)==1
		lgd = legend('P','P (SBP to 0.8 T)');
	else
		lgd = legend('P','P (SBP to 0.8 T)');
% 		plot(time([imax+1+i_intersec]),P(imax+1+i_intersec),'ro')
% 		lgd = legend('P','P (SBP to 0.8 T)','0-crossings for dP/dt');
	end
end

if isempty(i_intersec)==1
	%        [c imaxddP] = max(ddP_temp);
	%        i_vc = imax + imaxddP + 1;
	[c, imin] = min(dP_temp);
	x = [1:P_len];
	P_down = P(imax+imin);
	down = P_down + c*x;
	P_temp_2 = [P(1:imax+imin) down];
	P_temp_2 = P_temp_2(P_temp_2>=DBP);
	i_vc_max = length(P_temp_2);	%point where downslope projection meets DBP
	P_temp_2(end:P_len) = DBP;
	if Plots == 1
		plot(time([imax+imin:P_len]),P_temp_2([imax+imin:P_len])/Scale_P,'-.r','LineWidth',3)
		lgd.String{length(lgd.String)} = 'Downslope projection';
	end
	imin = imin+imax;
	i_vc_1 = round((imin+i_vc_max)/2);	%arbitrarily halfway between P_down and i_vc_max
	i_vc_2 = i_vc_max;					%assuming similarity with flow waveform
	
else
	i_vc_2 = imax + i_intersec(end) + 1;	%last dP 0-crossing
	[~, i_vc_1] = min(P(imax+1:i_vc_2));		%minimum P between SBP and last dP 0-crossing 
	i_vc_1 = imax + i_vc_1;
end
ind_dic = i_vc_1;
ind_dia = i_vc_2;

if Plots == 1
	plot(time(ind_dic),P(ind_dic)/Scale_P,'vr','MarkerSize',10,'MarkerFaceColor','r')
	plot(time(ind_dia),P(ind_dia)/Scale_P,'^r','MarkerSize',10,'MarkerFaceColor','r')
	lgd.String{length(lgd.String)-1} = 'Dicrotic notch';
	lgd.String{length(lgd.String)} = 'LVET';
	xlabel('Time [s]')
	ylabel('Pressure [mmHg]')
	xlim([0,T])
% 	ylim([70,130])
	set(gca,'XTick',[0:0.2:T+0.2])
% 	set(gca,'YTick',[])
	box on
	legend boxoff
	set(gca,'FontSize',40)
	addpath([PATHS.Root,'Others/PlotSave/'])
	PlotSave(PATHS.Figures,'LVET_dP_analysis_P')
	
% 	figure,	hold  on
% 	plot(time(1:end-1),dP/Scale_P,'-k','LineWidth',2)
% 	plot(time([imax+1:round(0.8*P_len)-1]),dP_temp/Scale_P,'--b','LineWidth',4)
% 	plot(time(ind_dia),dP(ind_dia)/Scale_P,'^r','MarkerSize',10,'MarkerFaceColor','r')
% 	legend('dP/dt','dP/dt (SBP to 0.8 T)','LVET');
% 	xlabel('Time [s]')
% 	ylabel('dP/dt [mmHg/s]')
% 	xlim([0,T])
% % 	ylim([70,130])
% 	set(gca,'XTick',[0:0.2:T+0.2])
% % 	set(gca,'YTick',[])
% 	box on
% 	legend boxoff
% 	set(gca,'FontSize',40)
% 	PlotSave(PATHS.Figures,'LVET_dP_analysis_dP')
end

end

function [intersections] = intersection (signal)
i=1;
intersections = [];
for k=1:length(signal)-1
	if signal(k)*signal(k+1)<0
		intersections(i) = k;
		i = i+1;
	end
end
end