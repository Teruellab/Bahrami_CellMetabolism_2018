datadir = 'D:\Experiments\20171212-5EU\Data\';

load([datadir,'firstreplicate.mat'])
firstrep = datamat2;

load([datadir,'secondreplicate.mat'])
secondrep = datamat;
datacols = [2,3];
counter = 0 ;
labelmat = {'FABP4 - replicate 1','CEBP\alpha - replicate 1','FABP4 - replicate 2','CEBP\alpha - replicate 2'};

for rep = 1:2
    if rep == 1
        datamat = firstrep;
    else
        datamat = secondrep;
    end
    for col = 1:2
        counter = counter+1;
        c = datacols(col);
        timerange = datamat(1,1):.1:datamat(end,1);

        modelfit = fit(datamat(:,1),datamat(:,c),'exp1');    
        confidence = confint(modelfit);
        lowerint = confidence(1,:);
        upperint = confidence(2,:);
        lb  = lowerint(1).*exp(lowerint(2).*timerange);
        ub = upperint(1).*exp(upperint(2).*timerange);
        halflife = log(2)/-modelfit.b;
        llife = log(2)/-lowerint(2);
        ulife = log(2)/-upperint(2);
        X = [timerange, fliplr(timerange)];
        F = [ub, fliplr(lb)];
        figure
        modelhandle = plot(modelfit);
        hold on 
        patch(X,F,'r','FaceAlpha',0.2,'EdgeColor','none')
        plot(datamat(:,1),datamat(:,c),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
        modelhandle.LineWidth = 1.5;
        set(gca,'YLim',[0 1.1])
        text(0.3, 0.2, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'])
        xlabel 'Time (hrs')
        ylabel 'Normalized Signal'
        title(labelmat{counter})
    end
end