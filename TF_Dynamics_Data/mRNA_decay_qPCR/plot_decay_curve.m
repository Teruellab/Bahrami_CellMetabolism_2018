
function plot_decay_curve(name,data,timerange,modelfit)

timevec = data(:,1);
chxvec = data(:,2);

confidence = confint(modelfit);
lowerint = confidence(1,:);
upperint = confidence(2,:);
lb  = lowerint(1).*exp(lowerint(2).*timerange);
ub = upperint(1).*exp(upperint(2).*timerange);
halflife = log(2)/-modelfit.b;
llife = log(2)/-lowerint(2);
ulife = log(2)/-upperint(2);
% fitintervals = predint(modelfit,timerange,.95,'functional','on');
% lb = fitintervals(:,1)';
% ub = fitintervals(:,2)';
X = [timerange, fliplr(timerange)];
F = [ub, fliplr(lb)];
% figure
patch(X,F,'r','FaceAlpha',0.2,'EdgeColor','none')
hold on 
modelhandle = plot(modelfit);
plot(timevec,chxvec,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
modelhandle.LineWidth = 1.5;
text(2.5, 0.75, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'])
legend off
title(name)
ylabel 'Normalized Value'
xlabel 'Time (hours)'
set(gca,'YLim',[0 2],'YTick',[0:0.5:2],'XLim',[-0.1 6.1],'XTick',[0:1:6],'FontName','Arial','FontSize',20)