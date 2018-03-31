imagedir = 'D:\Experiments\';
experimentdir = [imagedir,'20171113-stimduration\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 2:7; numrows = numel(rows);
cols = [2:11]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = 1:2; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

percdiffplate  = zeros(numrows,numcols,numplates);

for plate = 1:numplates
figure;
set(gcf,'color','w');
for col = 1:numcols
    tempcondition = [];
    ppargThresholdGuess = 1e3;
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
        badDna = dnawell< 6e5| dnawell>2.3e6;
        badPparg = ppargwell > 1e4;
        badcells =  badPparg | badDna;
        dnawell(badcells) = [];
        ppargwell(badcells) = [];
        percdiffplate(row,col,plate) = (sum(ppargwell>ppargThresholdGuess)*100)/numel(ppargwell);
        figure(gcf)
        subtightplot(numrows,numcols,posmat(row,col),[0.03 0.03])
        histogram(ppargwell); 
%         set(gca,'XLim',[0 300],'YLim',[0 1000]);
%         hold on
        line([ppargThresholdGuess ppargThresholdGuess],get(gca,'YLim'),'color','red','LineWidth',.5)

        
    end

end

end

dmivec = [percdiffplate(1:3,:,1),percdiffplate(4:6,1:3)];
rosivec =  [percdiffplate(1:3,:,2),percdiffplate(4:6,1:3,2)];

dmimean = mean(dmivec);
dmistd = std(dmivec);
rosimean = mean(rosivec);
rosistd = std(rosivec);
figure
barwitherr([dmistd./sqrt(3);rosistd./sqrt(3)]',[dmimean;rosimean]')
ylabel '% Differentiated Cells'
xlabel 'Stimulus Duration (hours)'
set(gca,'XTickLabel',{'0','2','4','6','8','10','12','14','16','18','24','36','48'},'XTickLabelRotation',45)
set(gca,'FontName','Arial','FontSize',20)
