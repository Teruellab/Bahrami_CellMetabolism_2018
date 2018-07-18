function dscatter_diffcells(xdim,ydim,thresholdVal,binrange)

% xmin = prctile(xdim,.1); xmax = prctile(xdim,99.9);
% ymin = prctile(ydim,.1); ymax = prctile(ydim,99.9);

xmin = binrange(1); xmax = binrange(end);
ymin = binrange(1); ymax = binrange(end);

redcolor = [231,76,60]./255;
bluecolor = [8,81,156]./255;

subtightplot(3,3,[4,8],[0.04 0.02],[0.07 .05], [0.05 0.05])
axh = dscatter(xdim(ydim>thresholdVal),ydim(ydim>thresholdVal),'msize',15);
colorscale1 = axh.Children.CData;
colorscale1 = colorscale1.*0.7 + 0.3;
redscale = [colorscale1.*redcolor(1),colorscale1.*redcolor(2),colorscale1.*redcolor(3)];

hold on
axh = dscatter_color(xdim(ydim<thresholdVal),ydim(ydim<thresholdVal),'msize',15);
colorscale2 = axh.Children(1).CData;
colorscale2 = colorscale2.*0.7 + 0.3;
bluescale = [colorscale2.*bluecolor(1),colorscale2.*bluecolor(2),colorscale2.*bluecolor(3)];

axh.Children(1).CData = bluescale;
axh.Children(2).CData = redscale;

bh = colorbar;
bh.TickLabels = '';


line([7 7],get(gca,'YLim'),'Color','k','LineStyle','--','LineWidth',1.5);
line(get(gca,'XLim'),[thresholdVal thresholdVal],'Color','k','LineStyle','--','LineWidth',1.5);

% [xs,ys ]= bin_scatter(xdim,ydim);
% plot(xs,ys,'LineWidth',2,'Color','k')
% scatter(xdim,ydim,'.')
set(gca,'FontName','Arial','FontSize',24,'XLim',[xmin xmax],'YLim',[ymin ymax])
% set(gca,'FontSize',24,'XLim',[4 10],'YLim',[4 10])

subtightplot(3,3,[1,2],[0.01 0.02],[0.07 .05], [0.05 0.05])
histogram(xdim(ydim>thresholdVal),binrange,'FaceColor',[231,76,60]./255);
hold on
histogram(xdim(ydim<thresholdVal),binrange,'FaceColor',[8,81,156]./255);

[bins1,edges1] = histcounts(xdim,binrange);
edges1 = edges1(1:end-1) + diff(edges1)./2;
plot(edges1,bins1,'Color','k','LineWidth',2,'LineStyle','--');

set(gca,'XTickLabel','','XTick',[],'YTick',get(gca,'YLim'),'XLim',[xmin xmax]);
hold off

subtightplot(3,3,[6,9],[0.01 0.02],[0.07 .05], [0.05 0.05])
histogram(ydim(ydim>thresholdVal),binrange,'FaceColor',[231,76,60]./255);
hold on
histogram(ydim(ydim<thresholdVal),binrange,'FaceColor',[8,81,156]./255);

[bins2,edges2] = histcounts(ydim,binrange);
edges2 = edges2(1:end-1) + diff(edges2)./2;
plot(edges2,bins2,'Color','k','LineWidth',2,'LineStyle','--');
view(90,-90)
set(gca,'XTickLabel','','XTick',[],'YTick',get(gca,'YLim'),'XLim',[ymin ymax]);
hold off
set(gcf,'color','white')

end
