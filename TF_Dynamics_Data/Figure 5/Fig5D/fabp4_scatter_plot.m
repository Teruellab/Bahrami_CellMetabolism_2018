% declare 
imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure 5\Fig5D\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = [4]; numrows = numel(rows);
cols = [4,11]; numcols = numel(cols); %24 hrs: cols 2-6, 48hrs: cols 7-11
sites = 1:9; numsites = numel(sites);
plates = [1]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

pulsepparg = [];
pulsefabp4 = [];
contpparg = [];
contfabp4 = [];


for plate = 1:numplates
figure(1)
set(gcf,'color','white')

    for col = 1:numcols


        for row = 1:numrows
          
            ppargwell = [];
            fabp4well = [];
            dnawell = [];

            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                
                pparg = fixdata(:,6);
                fabp4 = fixdata(:,9);            

                ppargwell = [ppargwell;pparg];
                fabp4well = [fabp4well;fabp4];
                dnawell = [dnawell;fixdata(:,4)];
                

            end
            badpparg =  ppargwell>500;
            badfabp4 = isnan(fabp4well) | fabp4well<0 | fabp4well>3*prctile(fabp4well,99);
            baddna = dnawell<2*10^5 | dnawell>1.5*10^6;
            badcells = badpparg | badfabp4 | baddna;
            ppargwell(badcells) = [];
            fabp4well(badcells) = [];
            dnawell(badcells) = [];
            
            if cols(col)<7
                ppargThreshold = 40;
                pulsepparg = [pulsepparg; ppargwell];
                pulsefabp4 = [pulsefabp4;fabp4well];
            else
                ppargThreshold = 40;
                contpparg = [contpparg;ppargwell];
                contfabp4 = [contfabp4;fabp4well];
            end            
            
            figure(1), hold on
            
            subtightplot(numrows,numcols,posmat(row,col),[0.05 0.05])
            dscatter(ppargwell,log2(fabp4well))
            axis square
            hold on
            line([ppargThreshold ppargThreshold],get(gca,'YLim'),'Color','k','LineWidth',1)
            line(get(gca,'XLim'),[12 12],'Color','k','LineWidth',1)
            set(gca,'XLim',[0 300],'YLim',[8 15],'FontName','Arial')
            hold off

           percdiffplate(row,col) = sum(ppargwell>ppargThreshold)*100/numel(ppargwell);

        end

    end


end
