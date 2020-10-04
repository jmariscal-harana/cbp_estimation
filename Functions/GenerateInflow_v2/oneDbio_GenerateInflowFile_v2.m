%%	Generate the .bcs inflow file
%	Input
%	-Path:		folder where time-flow data is stored
%	-File_Name: self explanatory
%	-Plots:		plots on/off (1/0)
%
%	Output
%	-t_H:	time vector from harmonics
%	-Q_H:		flow vector from harmonics
%	-T:		cardiac period extracted from 't_H'
%	-HR:	heart rate extracted from 't_H'
%	-SV:	stroke volume extracted from 'Q_H'
%	
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (23/01/19)
%
%==========================================================================

function [t_H,Q_H,T,HR,SV] = oneDbio_GenerateInflowFile_v2(Path,File_Name,Plots,t,Q)
%%	Import flow rate Q and time t
if nargin == 5
else
	data	= dlmread([Path,File_Name]);
	
	t		= data(:,1); %[s]
	Q		= data(:,2); %[mL/s]
end

%%	Smooth using a piecewise cubic spline interpolation
t_fine = linspace(t(1),t(end),length(t)*10);
Q_fine = interp1(t,Q,t_fine,'spline');

if Plots == 1
	figure, hold on
	plot(t,Q,'b')
	plot(t_fine,Q_fine,'r')
	hold off;
end


%% Calculate the approximation to the inflow given by N_Harm
N_Harm	= floor((length(t_fine)-1)/2);
dt		= t_fine(2) - t_fine(1);
SR		= round(1/dt);

[t_H,Q_H,T] = Harmonics_filter_v2(t_fine,Q_fine,N_Harm,SR,Plots,[Path,File_Name,'_IN_1.bcs']);
fprintf('Nektar .bcs inflow file generated\n');

if Plots == 1
	figure, hold on
	plot(t_fine,Q_fine)
	plot(t_H,Q_H)
	hold off
end


%%	Calculate SV and HR
HR	= 60/T;	%[bpm]
dt	= t_H(2) - t_H(1);
SV	= trapz(Q)*dt;	%


%% Plots
% if Plots == 1
% 	
% 	FontSize = 24;
% 	
% 	figure;
% 	hs1 = subplot(1,1,1);
% 	hold on;
% 	hp2 = plot(t_fine,Q_fine,'r');
% 	hp3 = plot(t_H,Q_H,'b');
% 	hp1 = plot(t,Q,'k .');
% 	hp4 = plot([DATA.Q.tpeak DATA.Q.tpeak], [0 1.1*max(Q_fine)], 'k --');
% 	
% 	set(hp1,'MarkerFaceColor','k'); set(hp1,'MarkerEdgeColor','k'); set(hp1,'MarkerSize',22.0);
% 	set(hp2,'LineStyle','-'); set(hp2,'LineWidth',2.0);
% 	set(hp3,'LineStyle','-'); set(hp3,'LineWidth',2.0);
% 	
% 	set(hs1,'FontSize',FontSize);
% 	hx = xlabel('t (s)','FontSize',FontSize);
% 	hy = ylabel('Q (ml/s)','FontSize',FontSize);
% 	
% 	grid on;
% 	hold off
% 	axis tight
% 	
% 	%set(hs1,'XTick',[0 0.2 0.4 0.6 0.8 1])
% 	%set(hs1,'YTick',[0; 15; 30; 45])
% 	
% 	set(hs1,'Box','off')
% 	
% 	hl = legend('interpolated', 'harmonics', 'raw', 'tpeak');
% 	set(hl,'FontSize',FontSize);
% 	
% 	%	Plot inflow waveform
% 	figure
% 	set(gcf,'name','In vivo inflow')
% 	
% 	hold on
% 	inflow1 = plot(t,Q,'k o');
% 	inflow2 = plot(DATA.Q.tinflow,DATA.Q.Qinflow,'k');
% 	hold off
% 	
% 	set(inflow1,'MarkerFaceColor','r'); set(inflow1,'MarkerEdgeColor','k'); %set(hp5,'MarkerSize',15.0);
% 	set(inflow2,'LineWidth',1.5);
% 	
% 	grid on
% 	box off
% 	
% 	axis([0,ceil(T),-50,500])
% 	set(gca,'XTick', 0:0.2:ceil(T), 'YTick', 0:100:500)
% 	set(gca,'FontSize',FontSize);
% 	hx = xlabel('t (s)','FontSize',FontSize);
% 	hy = ylabel('Q (ml/s)','FontSize',FontSize);
% 	hl = legend('Raw data', 'Interpolated');
% 	legend('boxoff')
% 	set(hl,'FontSize',18);
% 	
% 	%	Save images and figures
% 	print('In_vivo_inflow','-depsc')
% 	savefig('In_vivo_inflow')
% 
% end

