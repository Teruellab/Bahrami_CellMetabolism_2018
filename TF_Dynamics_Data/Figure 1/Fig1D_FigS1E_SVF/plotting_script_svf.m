% declare 
imagedir = 'z:\michael\';
experimentdir = [imagedir,'Fixed Cell Experiments\2017\20170307-SVF-diff\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 2:7;numrows = numel(rows);
cols = [2,7:10]; numcols = numel(cols); 
% Fig 1D: cols = [2,7:10] Fig S1E cols = [2:6];
sites = 1:12; numsites = numel(sites);
plates = [1:3]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

zeromat = zeros(numrows,numcols,numplates);
percdiffplate = zeromat;
numdiffplate = zeromat;
bodipycountplate = zeromat;
bodipysumplate = zeromat;
bodipymeanplate = zeromat;
boipyhighplate = zeromat;


for plate = 1:numplates
% figure(1)
% set(gcf,'color','white')

    for col = 1:numcols


        for row = 1:numrows
          
            ppargwell = [];
            bodipywell = [];
            dnawell = [];
            
            bodipycount = 0;
            bodipysum = 0;

            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                
                pparg = fixdata(:,10);
                bodipy = fixdata(:,9);

                ppargwell = [ppargwell;pparg];
                bodipywell = [bodipywell;bodipy];
                dnawell = [dnawell;fixdata(:,4)];
                
                bodipycount = bodipycount + bodipydata(1); % count of bright bodipy spots
                bodipysum = bodipysum + bodipydata(3); % sum signal of bodipy spots

                

            end
            badpparg = ppargwell<0 | ppargwell>1e5;
            badbodipy = isnan(bodipywell) | bodipywell<0 | bodipywell>3*prctile(bodipywell,99);
            baddna = dnawell<2.5e4 | dnawell>1.5e5;
            badcells =   baddna | badpparg | badbodipy;
            ppargwell(badcells) = [];
            bodipywell(badcells) = [];
            dnawell(badcells) = [];
            
%             figure(1), hold on
% % %             
%             subtightplot(numrows,numcols,posmat(row,col),[0.05 0.05])
%             
%             histogram(log2(bodipywell),[-5:.5:20])
%             set(gca,'XLim',[0 12],'YLim',[0 25])
%             dscatter(log2(ppargwell+1),log2(bodipywell+1))
%             line(get(gca,'XLim'),[13 13],'Color','k','LineWidth',1);
%             line([5 5],get(gca,'YLim'),'Color','k','LineWidth',1);
%             set(gca,'XLim',[0 12],'YLim',[0 22])

            percdiffplate(row,col,plate) = sum(ppargwell>2^5)*100/numel(ppargwell);
            bodipyhighplate(row,col,plate) = sum(bodipywell>2^13)*100/numel(ppargwell);
            numdiffplate(row,col,plate) = sum(ppargwell>2^5);
            bodipycountplate(row,col,plate) = bodipycount;
            bodipysumplate(row,col,plate) = bodipysum;
        end

    end

end



%% plot bar graphs for continuous versus pulse
r = [[231,76,60]./255];g = [77,175,74]./255; b = [[8,81,156]./255]; gy =[.4 .4 .4];

ylim1 = [0 20]; ystep1 = ylim1(1):(ylim1(2)-ylim1(1))/5:ylim1(2);
ylim2 = [0 3e9]; ystep2 = ylim2(1):(ylim2(2)-ylim2(1))/5:ylim2(2);
for p = 1:3
   sigmat = bodipyhighplate(:,:,p);
   bodipymat = bodipysumplate(:,:,p);
   
   control = mean(sigmat(:,1));
   cont = fliplr(sigmat(1:3,2:end));
   pulse = fliplr(sigmat(4:6,2:end));
   
   meancont = mean(cont,1); semcont = std(cont)/sqrt(3);
   meanpulse = mean(pulse,1); sempulse = std(pulse)/sqrt(3);
   
   figure
   subtightplot(2,2,1,[.07 .07],[.05 .01],[.01 .01])
   h1 = barwitherr(semcont,meancont);
   h1.FaceColor = r;
   set(gca,'YLim',ylim1,'YTick',ystep1,'XLim',[0.5 4.5],'XTick',[1:4],'XTickLabel',{'12hr','24hr','36hr','48hr'})
   line(get(gca,'XLim'),[control control],'LineStyle','--','Color','k','LineWidth',1)
   axis square
   ylabel '%Bodipy high cells'
   subtightplot(2,2,2,[.07 .07],[.05 .01],[.01 .01])
   h2 = barwitherr(sempulse,meanpulse);
   h2.FaceColor = b;
   set(gca,'YLim',ylim1,'YTick',ystep1,'XLim',[0.5 4.5],'XTick',[1:4],'XTickLabel',{'12hr','24hr','36hr','48hr'})
   line(get(gca,'XLim'),[control control],'LineStyle','--','Color','k','LineWidth',1)
   axis square
   ylabel '%Bodipy high cells'
   
   
   control = mean(bodipymat(:,1));
   cont = fliplr(bodipymat(1:3,2:end));
   pulse = fliplr(bodipymat(4:6,2:end));
   
   meancont = mean(cont,1); semcont = std(cont)/sqrt(3);
   meanpulse = mean(pulse,1); sempulse = std(pulse)/sqrt(3);
   
   
   subtightplot(2,2,3,[.07 .07],[.05 .01],[.01 .01])
   h1 = barwitherr(semcont,meancont);
   h1.FaceColor = r;
   set(gca,'YLim',ylim2,'YTick',ystep2,'XLim',[0.5 4.5],'XTick',[1:4],'XTickLabel',{'12hr','24hr','36hr','48hr'})
   line(get(gca,'XLim'),[control control],'LineStyle','--','Color','k','LineWidth',1)
   axis square
   ylabel 'Total Bodipy Intensity'
   
   subtightplot(2,2,4,[.07 .07],[.05 .01],[.01 .01])
   h2 = barwitherr(sempulse,meanpulse);
   h2.FaceColor = b;
   set(gca,'YLim',ylim2,'YTick',ystep2,'XLim',[0.5 4.5],'XTick',[1:4],'XTickLabel',{'12hr','24hr','36hr','48hr'})
   line(get(gca,'XLim'),[control control],'LineStyle','--','Color','k','LineWidth',1)
   axis square
   ylabel 'Total Bodipy Intensity'  
   
   
end

%% plot graphs for different pulsatile patterns
r = [[231,76,60]./255];g = [77,175,74]./255; b = [[8,81,156]./255]; gy =[.4 .4 .4];
ylim1 = [0 10]; ystep1 = ylim1(1):(ylim1(2)-ylim1(1))/5:ylim1(2);
ylim2 = [0 .5e9]; ystep2 = ylim2(1):(ylim2(2)-ylim2(1))/5:ylim2(2);

for p = 3
    ppargmat = bodipyhighplate(1:6,:,p);
    bodipymat = bodipysumplate(1:6,:,p);
    
    control = mean(ppargmat(4:6,1));
    cont = ppargmat(1:3,2:end);
    meancont = mean(cont,1); semcont = std(cont)/sqrt(3);
    figure
    subtightplot(1,2,1,[.01 .08],[.05 .05],[.05 .01])
    h1 = barwitherr(semcont,meancont);
    h1.FaceColor = b;
    set(gca,'YLim',ylim1,'YTick',ystep1,'XLim',[0.5 4.5],'XTick',[1:4],...
        'XTickLabel',{'1x, 1 pulse\newline48 hours','1x, 4 pulses\newline12 hours each','2x, 4 pulses\newline6 hours each','2x, 2 pulses\newline12 hours each'},...
        'FontName','Arial','XTickLabelRotation',45)
    line(get(gca,'XLim'),[control control],'LineStyle','--','Color','k','LineWidth',1)
    axis square
    ylabel '%Bodipy high cells'
    
    control = mean(bodipymat(4:6,1));
    cont = bodipymat(1:3,2:end);
    meancont = mean(cont,1); semcont = std(cont)/sqrt(3);
    
    subtightplot(1,2,2,[.01 .08],[.05 .05],[.05 .01])
    h1 = barwitherr(semcont,meancont);
    h1.FaceColor = b;
    set(gca,'YLim',ylim2,'YTick',ystep2,'XLim',[0.5 4.5],'XTick',[1:4],...
        'XTickLabel',{'1x, 1 pulse\newline48 hours','1x, 4 pulses\newline12 hours each','2x, 4 pulses\newline6 hours each','2x, 2 pulses\newline12 hours each'},...
        'FontName','Arial','XTickLabelRotation',45)
    line(get(gca,'XLim'),[control control],'LineStyle','--','Color','k','LineWidth',1)
    axis square
    ylabel 'Total Bodipy Intensity'
end

