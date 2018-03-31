imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure S2\FigS2C\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 2:7; numrows = numel(rows);
cols = [2:11]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1:4]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';
ppargcondition = cell(numcols,1);
ppargplate= cell(numrows,numcols);
percdiffplate = zeros(numrows,numcols);
cellcount = percdiffplate;
binranges = 0:12:600;
platenames = {'Cort+IBMX','Cort, IBMX Constant','Dex+IBMX','Dex, IBMX Constant'};
for plate = 1:numplates
for col = 1:numcols
    tempcondition = [];
    ppargThresholdGuess = 140;
    

    for row = 1:numrows
        
  
        ppargwell = [];
        dnawell = [];
        
        for site = 1:numsites
            
            shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
            datafile = [datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat'];
            if ~exist(datafile)
                continue
            end
            try
                load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
            catch
                continue
            end
            pparg = fixdata(:,6);
            dnacontent = fixdata(:,4);
            ppargwell = [ppargwell; pparg];
            dnawell = [dnawell; dnacontent];
            
        end
        cellcount(row,col) = numel(dnawell);
        badDna = dnawell> 500000 | dnawell<100000;
        badPparg = ppargwell > 1000;
        badcells =  badDna;
        dnawell(badcells) = [];
        ppargwell(badcells) = [];
        percdiffplate(row,col) = (sum(ppargwell>ppargThresholdGuess)*100)/numel(ppargwell);
        ppargplate{row,col} = ppargwell;

%         
    end

end
meanctrl = mean(percdiffplate(:,1)); stdctrl = std(percdiffplate(:,1));
finaldiff = [percdiffplate(1:3,2:end),percdiffplate(4:6,2:3),percdiffplate(4:6,5)];

dmimean = mean(finaldiff(1:3,:),1);
dmistd = std(finaldiff(1:3,:),0,1);
% 
figure
set(gcf,'color','white')
barwitherr(dmistd(1:end),dmimean(1:end),'FaceColor','blue')
hold on
line(get(gca,'XLim'),[meanctrl meanctrl],'LineStyle','--','Color','k','LineWidth',2)
ylabel('% Differentiated')
title(platenames{plate})
set(gca,'XTickLabel',{'2h','4h','6h','8h','10h','12h','14h','16h','18h','24h','36h','48h'})
set(gca,'FontSize',20)
hold off

end

