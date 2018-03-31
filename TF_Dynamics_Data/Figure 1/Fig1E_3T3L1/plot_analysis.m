datadir = ['z:\michael\Fixed Cell Experiments\2015\20151202-3T3-Pulsing\10x\Data\'];

rows = 2:7;
cols = 2:11;
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
%             load([datadir,'fixdata_',shot,'.mat'])
            if ~exist([datadir,'fixdata_',shot,'.mat'],'file')
                continue
            end
            load([datadir,'fixdata_',shot,'.mat'])
            
            %%%%%%%Structure of fixdata.mat%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fixdata =
            % [nucx,nucy,nuc_area,nuc_totalIntensity,nuc_medianIntensity,
            % pparg_medianIntensity];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            
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
        
        figure(2)
        subtightplot(numrows,numcols,posmat(row,col))
        histogram(dnawell)
        
        percdiff(row,col) = 100*(sum(ppargwell>ppargthreshold))/numel(ppargwell);
        
    end
end

ctrl  = median(percdiff(:,2));

pvc = [percdiff(:,3),percdiff(:,7)];
figure,bar(median(pvc(:,1)),'FaceColor','red')
set(gca,'YLim',[0 100])
hold on
h1 = bar(median(pvc(:,2)));
h1.XData = 2;
set(gca,'XTick',1:2,'XTickLabel',{'1x, 1 pulse\newline48 hours','1x, 4 pulses\newline 12 hours each'})
line(get(gca,'XLim'),[ctrl ctrl],'LineStyle','--','Color','k','LineWidth',2)
set(gca,'FontSize',24)
hold off

cont = [percdiff(4:6,4),percdiff(4:6,5),percdiff(4:6,6),percdiff(4:6,10)];
pulse = [percdiff(4:6,7),percdiff(4:6,8),percdiff(4:6,9),percdiff(4:6,10)];
figure,barwitherr(fliplr(std(cont,0,1)),fliplr(mean(cont,1)),'FaceColor','Red')
set(gca,'YLim',[0 80],'XLim',[0 5])
hold on
line(get(gca,'XLim'),[ctrl ctrl],'LineStyle','--','Color','k','LineWidth',2)
set(gca,'XTickLabel',{'12h','24h','36h','48h'},'FontSize',24)
figure, notBoxPlot(cont)

figure,barwitherr(fliplr(std(pulse,0,1)),fliplr(mean(pulse,1)))
set(gca,'YLim',[0 80],'XLim',[0 5])
hold on
line(get(gca,'XLim'),[ctrl ctrl],'LineStyle','--','Color','k','LineWidth',2)
set(gca,'XTickLabel',{'12h','24h','36h','48h'},'FontSize',24)
figure,notBoxPlot(pulse)
