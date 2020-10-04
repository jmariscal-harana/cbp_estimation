close all; clear; clc;

time = (0:999)/999;
width = 199;

figure
filename = ['PulseWave','.gif'];

j = 1;

DelayTime = 1/20;	% the smaller, the quicker it plays on a presentation

for jj = 1:length(time)+width
	
	wave(1:length(time))	= 0;
	
	if jj <= width
		wave(1:jj)	= sin(pi*(5*time(width-jj+1:width)));
	elseif 	jj > length(time)
		wave(jj-width:length(time))	= sin(pi*(5*time(1:width-jj+length(time)+1)));
	else
		wave(jj-width:jj)			= sin(pi*(5*time(1:200)));
	end
	pause(0.005)
	plot(time,wave,'r','LineWidth',3)
	set(gcf, 'color', [1 1 1])
	
% 	xlabel('axial coordinate','FontSize',24)
	set(gca,'fontsize',24)
	ylim([-0.2,1.2]);
	ax = gca;
	ax.XTick = [];
	ax.YTick = [];
	
	drawnow
	frame = getframe(1);
	im = frame2im(frame);
	[imind,cm] = rgb2ind(im,256);
	
	if j == 1 && jj == 1
		imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',DelayTime);
	elseif mod(jj,5) == 0
		imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DelayTime);
	end
	
end