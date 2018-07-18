% declare 
imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure S5\FigS5H_5G\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = [2]; numrows = numel(rows);
cols = [8:9]; numcols = numel(cols);
sites = 1:9; numsites = numel(sites);
plates = [1:5]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';


for plate = 1:numplates


    for col = 1:numcols


        for row = 1:numrows
          
            yfptagwell = [];
            ppargifwell = [];
            cebpbifwell = [];
            dnawell = [];

            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                yfptag = fixdata(:,6);
                ppargif = fixdata(:,8);
                cebpbif = fixdata(:,7);               

                yfptagwell = [yfptagwell;yfptag];
                ppargifwell = [ppargifwell;ppargif];
                cebpbifwell = [cebpbifwell;cebpbif];
                dnawell = [dnawell;fixdata(:,4)];
                

            end
            badyfptag = yfptagwell>300;
            badppargif =  ppargifwell>300;
            badcebpbif = cebpbifwell>500;
            if cols(col)>4 && cols(col)<8
                
                baddna = dnawell<5*10^4 | dnawell>3*10^5;
            else
                baddna = dnawell<2*10^4 | dnawell>2*10^5;
            end
            badcells = badppargif | badcebpbif | baddna | badyfptag;
            yfptagwell(badcells) = [];
            ppargifwell(badcells) = [];
            cebpbifwell(badcells) = [];
            dnawell(badcells) = [];
            
            yfpmean(row,col,plate) = mean(yfptagwell);
            ppargmean(row,col,plate) = mean(ppargifwell);
            cebpbmean(row,col,plate) = mean(cebpbifwell);
%             
            figure(plate)
%             
            subtightplot(numrows,numcols,posmat(row,col),[0.05 0.05])
            dscatter(cebpbifwell,yfptagwell)
            set(gca,'YLim',[0 350],'XLim',[0 350])
            axis square
            ylabel 'CEBPB-citrine'
            xlabel 'CEBPB-IF'
            
        end

    end


end

citrine = squeeze(yfpmean);
IF = squeeze(cebpbmean);
figure,barwitherr([std(IF);std(citrine)]',[mean(IF);mean(citrine)]')
legend({'IF signal','Citrine Signal'})
set(gca,'FontName','Arial','FontSize',15,'XLim',[0 6],'XTick',1:5,...
    'XTickLabel',{'0hr','24hr','48hr','72hr','96hr'},'XTickLabelRotation',45)
