imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure 1\Fig1C_3T3L1\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 5:7; numrows = numel(rows);
cols = [2,3,7,10,11]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = 1; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';
ppargcondition = cell(numcols,1);
ppargplate= cell(numrows,numcols);
cellcount = zeros(numrows,numcols);
percdiffplate  = cellcount;

binranges = 0:6:300;
platethresholdmat = [125 125 125 125 125 125 125 125 125 100];
% PlateNames = {'HalfDose','1xDose','2xDose','4xDose','6xDose'};
for plate = 1:numplates
figure;
set(gcf,'color','w');
for col = 1:numcols
    tempcondition = [];
    ppargThresholdGuess = 50;
    

    for row = 1:numrows
        

        ppargwell = [];
        dnawell = [];
        
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
%         dnawell(badcells) = [];
%         ppargwell(badcells) = [];
        tempcondition = [tempcondition; ppargwell];
        percdiffplate(row,col) = (sum(ppargwell>ppargThresholdGuess)*100)/numel(ppargwell);
        ppargplate{row,col} = ppargwell;
        figure(gcf)
        subtightplot(numrows,numcols,posmat(row,col),[0.03 0.03])
        histogram(ppargwell,binranges,'EdgeColor','none'); 
%         set(gca,'XLim',[0 600],'YLim',[0 1000]);

%         [f,x] = ksdensity(ppargwell);
%         plot(x,f,'LineWidth',2)
        hold on
        line([ppargThresholdGuess ppargThresholdGuess],get(gca,'YLim'),'color','red','LineWidth',.5)
        hold off
        title(num2str(percdiffplate(row,col)));
        

    end

end


end

