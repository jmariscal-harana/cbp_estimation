%%	EXAMPLE:
% 	PLOTS.Format = 1;
% 	PLOTS.FontSize = 20;
% 	PLOTS.Title = 1;
% 	PLOTS.Title_text = 'Title text';
% 	PLOTS.Legend = 1;
% 	PLOTS.Legend_box = 0;
% 	PLOTS.Legend_text = 'Legend text';
% 	PLOTS.LegendSize = 20;
% 	PLOTS.XLabel = 1;
% 	PLOTS.XLabel_text = 'X-label text';
% 	PLOTS.YLabel = 1;
% 	PLOTS.YLabel_text = 'Y-label text';
% 	PLOTS.Annotation = 1;
% 	PLOTS.FontSizeAnno = 20;
% 	PLOTS.LineWidth = 2;
% 	PLOTS.BA = 1;
% 	PLOTS.XLim = 1;
% 	PLOTS.XLim_values = [0 1];
% 	PLOTS.XTick = 0.2;
% 	PLOTS.YLim = 1;
% 	PLOTS.YLim_values = [-1 1];
% 	PLOTS.YTick = 0.5;
% 	PLOTS.Grid = 1;
% 	PLOTS.Markers = {'ko', 'r^'}
% 	PLOTS.MarkerSize = 10;
%
%	FormatPlots(PLOTS)


function [] = FormatPlots(PLOTS)

set(gca,'FontSize',PLOTS.FontSize)

%	TITLE
if PLOTS.Title == 1
	title(PLOTS.Title_text)
end

%	LEGEND
if PLOTS.Legend	== 1
	legend(PLOTS.Legend_text,'location','northeast','FontSize',PLOTS.LegendSize)
	if isfield(PLOTS,'Legend_box')
		if PLOTS.Legend_box == 0
			legend('boxoff')
		end
	end
end

%	X, Y LABELS
if PLOTS.XLabel == 1
	xlabel(PLOTS.XLabel_text)
end
if PLOTS.YLabel == 1
	ylabel(PLOTS.YLabel_text)
end

%	AXIS LIMITS
switch PLOTS.XLim
	case 0
		set(gca,'XTick',[])
	case 1
		xlim(PLOTS.XLim_values)
		set(gca,'XTick',PLOTS.XLim_values(1):PLOTS.XTick:PLOTS.XLim_values(2))
	case 2
		xlim('auto')
end

switch PLOTS.YLim
	case 0
		set(gca,'YTick',[])
	case 1
		ylim(PLOTS.YLim_values)
		set(gca,'YTick',ceil(PLOTS.YLim_values(1)/PLOTS.YTick)*PLOTS.YTick:PLOTS.YTick:floor(PLOTS.YLim_values(2)/PLOTS.YTick)*PLOTS.YTick)
	case 2
		ylim('auto')
end

%	GRID
if PLOTS.Grid == 1
	grid on
end

%	ANNOTATIONS
if PLOTS.Annotation == 1
	if PLOTS.BA == 1
		Text_BA{1} = num2str(round(PLOTS.md(1)+2*PLOTS.sd(1),1));
		Text_BA{2} = num2str(round(PLOTS.md(1),1));
		Text_BA{3} = num2str(round(PLOTS.md(1)-2*PLOTS.sd(1),1));
		
		Y_range = 4*PLOTS.sd(1);
		xt = [PLOTS.XLim_values(2)*0.9, PLOTS.XLim_values(2)*0.9, PLOTS.XLim_values(2)*0.9];
		yt = [PLOTS.md(1)+2*PLOTS.sd(1) + 0.2*Y_range, PLOTS.md(1) + 0.2*Y_range, PLOTS.md(1)-2*PLOTS.sd(1) + 0.2*Y_range];
		for jj = 1:length(Text_BA)
			t = text(xt(jj),yt(jj),Text_BA{jj});
			t.FontSize = PLOTS.FontSizeAnno;
		end
	end
	
elseif PLOTS.Annotation == 2
	if PLOTS.BA == 1
		Text_BA{1} = ['+ LoA: ',num2str(round(PLOTS.md(1)+2*PLOTS.sd(1),PLOTS.BA_sf))];
		Text_BA{2} = ['Bias:  ',num2str(round(PLOTS.md(1),PLOTS.BA_sf))];
		Text_BA{3} = ['- LoA: ',num2str(round(PLOTS.md(1)-2*PLOTS.sd(1),PLOTS.BA_sf))];
		
		dim = [0.65 0.85 0.1 0.1];	%upper right
% 		dim = [0.65 0.40 0.1 0.1];	%lower right
		str = {Text_BA{1},Text_BA{2},Text_BA{3}};
		t = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','white');
		t.FontSize = PLOTS.FontSizeAnno;

% 		Y_range = 4*PLOTS.sd(1);
% 		xt = PLOTS.XLim_values(2)*0.9;
% 		yt = PLOTS.md(1)+2*PLOTS.sd(1) + Y_range;
% 		xt = 0.9;
% 		yt = 0.9;
% 		t = text(xt,yt,Text_BA);
end

end