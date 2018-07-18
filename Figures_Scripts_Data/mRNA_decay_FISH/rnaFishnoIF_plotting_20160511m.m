imagedir = 'Y:\michael\Fixed Cell Experiments\';
experimentdir = [imagedir,'2016\2016 From Asus Lap top\20160511-OP9-RNA-Fish\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 2:5; numrows = numel(rows);
cols = [4,5]; numcols = numel(cols);
sites = 1:9; numsites = numel(sites);
plates = [1:7]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

figure(1)
% figure(2)
% fishcell = cell(numplates);
for plate = 1:numplates
% figure
if plates(plate)==1
    cols = [5,6];
else
    cols =[4,5];
end
    for col = 1:numcols


        for row = 1:numrows

            fish = [];

            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fishdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                fish = [fish; IFdata(:,11)];

            end
            qcmat(row,col) = sum(fish>0)*100/numel(fish);
            
            figure(1)
            hold on
            subtightplot(numrows,numcols,posmat(row,col),[0.06,0.01])
            hold on
            histogram(fish,'EdgeColor','none','FaceAlpha',0.3,'normalization','probability')
%             set(gca,'XLim',[0 10])
            hold off 
            
            fishmat(row,col,plate) = nanmean(fish);


        end

    end

end
%% initial plotting code

fishmatcombined = squeeze(mean(fishmat,2));
time = [0,.5,1,1.5,3,6,9];
condList = {'FABP4','CEBPB','CEBPA','PPARG'};
unstimvals = [31, 22, 5, 18];
maxlim = [140, 60, 120, 75];
for c = 2
   figure
%    subtightplot(1,2,1,[0.02 0.1])
   plot(time,fishmatcombined(c,[1:7]),'Marker','o','Color','k','LineWidth',2,'LineStyle','-')
   hold on
   plot(time,squeeze(fishmat(c,1,[1:7])),'Marker','o','Color','k','LineWidth',1,'LineStyle','--')
   plot(time,squeeze(fishmat(c,2,[1:7])),'Marker','o','Color',[.5 .5 .5],'LineWidth',1,'LineStyle','--')
%    axis square
   legend ({'Average','Rep. 1','Rep. 2'})
%    set(gca,'YLim',[0 200])
%    xlabel 'Time (minutes)'
%    ylabel 'mRNA puncta'
   title(condList{c})
%    set(gca,'FontName','Arial','FontSize',12)
%    hold off
%    modelfit = fit(time',fishmatcombined(c,[1,3:6])','exp1');
%    halflife = log(2)/-modelfit.b;
%    subtightplot(1,2,2,[0.02 0.02])
%    plot(time,fishmatcombined(c,[1,3:6]),'LineStyle','none','Marker','o','Color','k','LineWidth',1)
%    text(time(4),.25*maxlim(c),['Half-life: ',num2str(round(10*halflife)/10),' hours'])
%    axis square
%    title(condList{c})
%    set(gca,'FontName','Arial','FontSize',12)
%    set(gca,'YLim',[0 maxlim(c)])
%    hold on
%    handle = plot(modelfit);
%    handle.LineWidth = 1.5;
%    axis square
%    line(get(gca,'XLim'),[unstimvals(c) unstimvals(c)],'LineStyle','--','LineWidth',1,'Color','k')
%    title(condList{c})
%    xlabel 'Time (hours)'
%    ylabel 'Average mRNA puncta per cell'
%    legend ({'Data','Fit Line','Unstimulated OP9'})
   hold off
end
%% plotting code with confidence bounds 2016-11-21
time = [0 .5 1 1.5 3];
interval = 1:5;
timemat = repmat(time,2,1);
timevec = timemat(:);
timerange = 0:.01:3;
gene = {'FABP4','CEBPB','CEBPA','PPARG'};
for c = 2
    fishvals = squeeze(fishmat(c,:,interval));
    fishvec = fishvals(:);
    modelfit = fit(timevec,fishvec,'exp1');
    confidence = confint(modelfit);
    lowerint = confidence(1,:);
    upperint = confidence(2,:);
    lb  = lowerint(1).*exp(lowerint(2).*timerange);
    ub = upperint(1).*exp(upperint(2).*timerange);
    halflife = round(10*log(2)/-modelfit.b)/10;
    llife = round(10*log(2)/-lowerint(2))/10;
    ulife = round(10*log(2)/-upperint(2))/10;
    X = [timerange, fliplr(timerange)];
    F = [ub, fliplr(lb)];
    figure
    modelhandle = plot(modelfit);
    modelhandle.LineWidth = 1.5;
    hold on
    patch(X,F,'r','FaceAlpha',0.2,'EdgeColor','none')
    plot(timevec,fishvec,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','none')
    title(gene{c})
    xlabel 'Time (hours)'
    ylabel 'mRNA puncta'
    set(gca,'XLim',[-0.1 3.1],'XTick',[0:.5:3],'FontName','Arial','FontSize',20)
    text(0,5, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'],'FontName','Arial','FontSize',16)
    
end

