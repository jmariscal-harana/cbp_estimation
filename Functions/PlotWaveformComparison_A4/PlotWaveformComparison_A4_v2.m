function PlotWaveformComparison_A4_v2(PATHS,PLOTS,REF,EST,REF_Dataset,Sc)

Plots_page = PLOTS.Horizontal*PLOTS.Vertical;

figure
kk = 0;

for jj = 1:length(REF)
	kk = kk + 1;
	
	subplot(PLOTS.Vertical,PLOTS.Horizontal,kk);
	hold on
	plot(REF(jj).t,REF(jj).Pressure_CBP,'k','LineWidth',2)
	for ll = 1:length(EST)
		plot(EST{ll}(jj).t,EST{ll}(jj).Pressure, PLOTS.Markers{ll},'LineWidth',1)
	end
	hold off
	
	% 	YMin = min([min(REF(jj).Asc_Ao.ONE_CYCLE.P/133.25)
	% 		min(EST_1D(jj).Asc_Ao.ONE_CYCLE.P);
	% % 		min(EST_Wk2(jj).Pressure)
	% 		min(EST_Wk3(jj).Pressure)]);
	% 	YMax = max([max(REF(jj).Asc_Ao.ONE_CYCLE.P/133.25)
	% 		max(EST_1D(jj).Asc_Ao.ONE_CYCLE.P);
	% % 		max(EST_Wk2(jj).Pressure)
	% 		max(EST_Wk3(jj).Pressure)]);
	
	XMin = 0;
	XMax = round(max(REF(jj).t)+0.049,1);
	% 	XMax = 1.2;	
	xlim([XMin, XMax])
	ylim(PLOTS.YLim_values)
	YMin = PLOTS.YLim_values(1);
	YMax = PLOTS.YLim_values(2);
	
	XTick = [XMin, XMax];
	% 	XTick = [XMin:0.5:XMax];
	YTick = PLOTS.YTick;
	set(gca,'XTick',XTick)
	set(gca,'YTick',[ceil(YMin/YTick)*YTick:YTick:floor(YMax/YTick)*YTick])
	set(gca,'FontSize',PLOTS.FontSize)
	
	if mod(jj,PLOTS.Horizontal) == 1
		ylabel('P [mmHg]')
	end
	if jj > length(REF) - PLOTS.Horizontal || mod(jj-1,Plots_page) >= Plots_page - PLOTS.Horizontal
		xlabel('t [s]')
	end
	
	if jj == length(REF)
		print('-fillpage',[PATHS.Figures,'/',REF_Dataset,'_waveform_',Sc,'_',num2str(ceil(jj/Plots_page))],'-dpdf')
	elseif isreal(jj/Plots_page) && rem(jj/Plots_page,1) == 0
		print('-fillpage',[PATHS.Figures,'/',REF_Dataset,'_waveform_',Sc,'_',num2str(jj/Plots_page)],'-dpdf')
		figure
		kk = 0;
	end
end

end
