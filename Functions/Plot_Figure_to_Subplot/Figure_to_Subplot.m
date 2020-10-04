function Figure_to_Subplot(rows,columns,subplot_ID)
% Create destination graph
h =  findobj('type','figure');
figure_ID = length(h) + 1;
figure(figure_ID)
subplot_num = length(subplot_ID);

ax = zeros(subplot_num,1);
for i = 1:subplot_num
    ax(i)=subplot(rows,columns,i);
end
% Now copy contents of each figure over to destination figure
% Modify position of each axes as it is transferred
for i = 1:subplot_num
	figure(subplot_ID(i))
	h = get(gcf,'Children');
	newh = copyobj(h,figure_ID);
	possub  = get(ax(i),'Position');
	set(newh,'Position',...
		[possub(1) possub(2) possub(3) possub(4)])
	delete(ax(i));
end
figure(figure_ID)
