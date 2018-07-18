
imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure S2\FigS2D\'];
datadir = [experimentdir,'Data\'];

rows = 2:4;
cols = [2,4:11];
sites = 1:4;

numrows = numel(rows);
numcols = numel(cols);
numsites = numel(sites);

posmat = reshape(1:numrows*numcols,numcols,numrows)';

ppargthreshold = 60;
percdiff = zeros(numrows,numcols);

plate = 1;
for row = 1:numrows
    for col = 1:numcols
        ppargwell = [];
        dnawell = [];
        for site = 1:numsites
            
            shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
            datafile = [datadir,'fixdata_Plate_1_',shot,'.mat'];
            
            try

                load([datadir,'fixdata_Plate_1_',shot,'.mat'])
            catch
                continue
            end
            
            pparg = fixdata(:,6);
            dna = fixdata(:,4);
            ppargwell = [ppargwell;pparg];
            dnawell = [dnawell;dna];
            
            
            
        end
        figure(1)
        subtightplot(numrows,numcols,posmat(row,col))
        histogram(ppargwell,0:5:500)
        hold on
        line([ppargthreshold ppargthreshold],get(gca,'YLim'),'Color','k','LineStyle','--','LineWidth',1.5)
        
        percdiff(row,col) = 100*(sum(ppargwell>ppargthreshold))/numel(ppargwell);
        
    end
end
r = [[231,76,60]./255];
meanDiff = mean(percdiff(1:3,:));
meanStd = std(percdiff(1:3,:),0,1);
figure
barwitherr(meanStd(2:end),meanDiff(2:end),'FaceColor',r)
set(gca,'YLim',[0 40],'FontSize',16,'FontName','Arial','XLim',[0 9],...
'XTick', 1:8,'XTickLabel',{'2h','4h','8h','12h','18h','24h','36h','48h'})
line(get(gca,'XLim'),[meanDiff(1) meanDiff(1)],...
    'Color','k','LineStyle','--','LineWidth',1.5)
ylabel '% Differentiated Cells'
xlabel 'Stimulus Duration (hr)'
