% declare 
imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Figure S1\FigS1A-S1C\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = [4]; numrows = numel(rows);
cols = [2:7]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

bodipyall = [];
fatgeneall = [];
ppargall = [];
for plate = 1:numplates
% figure(plate)
% set(gcf,'color','white')

    for col = 1:numcols


        for row = 1:numrows
            ppargwell = [];
            bodipywell = [];
            fatgenewell = [];
            dnawell = [];
            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                
                fatgene = fixdata(:,13);
                bodipy = fixdata(:,9); 
                pparg = fixdata(:,14);
                
                fatgenewell = [fatgenewell;fatgene];
                bodipywell = [bodipywell;bodipy];
                ppargwell= [ppargwell;pparg];
                
                dnawell = [dnawell;fixdata(:,4)];
                

            end
            badbodipy = isnan(bodipywell) | bodipywell<0 | bodipywell >2^15;
            badfatgene = isnan(fatgenewell) | fatgenewell<0 | fatgenewell>2^12;
            badpparg = ppargwell<0 | ppargwell > 4000;
            baddna = dnawell<4e5 | dnawell >2e6;
            badcells = badfatgene | badbodipy  | badpparg | baddna;
            
            fatgenewell(badcells) = [];
            bodipywell(badcells) = [];
            ppargwell(badcells) = [];
            dnawell(badcells) = [];
            
            figure(3)
            subtightplot(numrows,numcols,posmat(row,col),[.01 .04],[.05 .05],[.04 .01])
            dscatter(ppargwell,log2(bodipywell))
            set(gca,'FontName','Arial','FontSize',18)
            if cols(col)>4
                set(gca,'YLim',[5 13],'YTick',[5:2:13],'XLim',[0 2500])
            else
                set(gca,'YLim',[5 13],'YTick',[5:2:13],'XLim',[0 1000])
            end
            
            xlabel 'PPARG (a.u.)'
            ylabel 'log2(Bodipy)'
            axis square

            figure(4)
            subtightplot(numrows,numcols,posmat(row,col),[.01 .04],[.05 .05],[.04 .01])
            dscatter(ppargwell,log2(fatgenewell))
            set(gca,'FontName','Arial','FontSize',18)
            
            xlabel 'PPARG (a.u.)'
            axis square
            
            if cols(col)>4
                ylabel 'log2(Adiponectin)'
                set(gca,'YLim',[5 13],'YTick',[5:2:13],'XLim',[0 2500])
            else
                ylabel 'log2(Glut4)'
                set(gca,'YLim',[5 13],'YTick',[5:2:13],'XLim',[0 1000])
            end
            
            
            
            ppargall= [ppargall;ppargwell];
            bodipyall = [bodipyall;bodipywell];
            fatgeneall = [fatgeneall;fatgenewell];
        end

    end


end

