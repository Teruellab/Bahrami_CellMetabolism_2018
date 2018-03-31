imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure S2\FigS2A-S2B\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 2:7; numrows = numel(rows);
cols = [3,4,5,8,9,10]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = 1:5; numplates = numel(plates);

posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';
ppargcondition = cell(numcols,1);
ppargplate= cell(numrows,numcols);
cellcount = zeros(numrows,numcols);
percdiffplate  = zeros(numrows,numcols,numplates);
binranges = 0:5:300;
titlevec = {'Lowest Titr Range','Low Titr Range','Mid Titr Range','High Titr Range','Highest Titr Range'};
for plate = 1:numplates
% figure;
% set(gcf,'color','w');

for col = 1:numcols
    tempcondition = [];
    for row = 1:numrows
        ppargwell = [];
        dnawell = [];
        ppargThresholdGuess = 50;
        for site = 1:numsites
            
            shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
            load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
            pparg = fixdata(:,6);
            dnacontent = fixdata(:,4);
            ppargwell = [ppargwell; pparg];
            dnawell = [dnawell; dnacontent];
            
        end
        cellcount(row,col) = numel(dnawell);
%         badDna = dnawell> 600000 | dnawell<100000;
        badPparg = ppargwell > 1000;
        badcells =  badPparg;
        dnawell(badcells) = [];
        ppargwell(badcells) = [];
        tempcondition = [tempcondition; ppargwell];
        percdiffplate(row,col,plate) = (sum(ppargwell>ppargThresholdGuess)*100)/numel(ppargwell);
        ppargplate{row,col} = ppargwell;
%         figure(gcf)
%         subtightplot(numrows,numcols,posmat(row,col),[0.03 0.03])
%         histogram(ppargwell,binranges,'EdgeColor','none'); 
%         set(gca,'XLim',[0 600],'YLim',[0 1000]);
%         hold on
%         line([ppargThresholdGuess ppargThresholdGuess],get(gca,'YLim'),'color','k','LineWidth',2)
%         hold off
%         title(num2str(percdiffplate(row,col)));
        
%         figure(2)
%         subtightplot(numrows,numcols,posmat(row,col),[0.03 0.03])
%         histogram(dnawell);
        
    end

end

controlmean = mean(percdiffplate(1,:,plate),2);
controlstd = std(percdiffplate(1,:,plate),0,2);

nocomp = percdiffplate(2:end,1:3,plate);
meanNoComp = mean(nocomp,2);
stdNoComp = std(nocomp,0,2);

comp = percdiffplate(2:end,4:6,plate);
meanComp = mean(comp,2);
stdComp = std(comp,0,2);
subtightplot(2,3,plate,[0.07 0.07],[0.05 0.02],.05)
barhandle = barwitherr([stdNoComp,stdComp],[meanNoComp,meanComp]);
barhandle(2).FaceColor = 'r';
set(gca,'YLim',[0 100],'XTickLabel',{'8hr','12h','16h','24h','48h'})
set(gca,'FontSize',20,'XLim',[0 6])
hold on
% line(get(gca,'XLim'),[controlmean+controlstd, controlmean+controlstd],'LineStyle','--','Color','k','LineWidth',1);
% line(get(gca,'XLim'),[controlmean-controlstd, controlmean-controlstd],'LineStyle','--','Color','k','LineWidth',1);
line(get(gca,'XLim'),[controlmean, controlmean],'LineStyle','--','Color','k','LineWidth',2);
hold off
set(gcf,'color','white')
title(titlevec{plates(plate)},'Units','normalized','Position',[.5 .9 1])
% figure,
% notBoxPlot([nocomp;comp]')

end


%% Plotting by time point
r = [231 76 60]./255; b = [8 81 156]./255;
averagemat = squeeze(percdiffplate(1,:,:));
bgdiff = mean(averagemat(:));
titlevec = {'8hr','12hr','16hr','24hr','48hr'};
for ind = 1:5
    test = squeeze(percdiffplate(ind+1,:,:))';
    nocomp = test(:,1:3);
    comp = test(:,4:6);
    figure
    handle = barwitherr([std(nocomp,0,2),std(comp,0,2)],[mean(nocomp,2),mean(comp,2)]);
    handle(1).FaceColor = r;
    handle(2).FaceColor = b;
    set(gca,'YLim',[0 100],'XTickLabel',{'0.083x','0.17x','0.33x','0.67x','1x'},'FontName','Arial','FontSize',14)
    xlabel 'Stimulus Dilution'
    ylabel '% Differentiated Cells'
    hold on
    line(get(gca,'XLim'),[bgdiff bgdiff],'Color','k','LineStyle','--','LineWidth',1)
    if ind==5
        legend({'Not Compensated','Compensated'})
    end
    title(titlevec{ind});
end


