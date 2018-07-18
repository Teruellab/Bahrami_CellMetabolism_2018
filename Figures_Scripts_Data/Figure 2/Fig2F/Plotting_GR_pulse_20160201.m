imagedir = 'z:\michael\';
experimentdir = [imagedir,'Fixed Cell Experiments\2016\20160121-GC-Pulse-3T3L1\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze

plates = 1:5; numplates = numel(plates);
colorchoice = {'r','r','b','b'};

rows = 3:6;
cols = 4:9;
sites = 1:4; numsites = numel(sites);

numrows = numel(rows); numcols = numel(cols);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';


yfpmat = zeros(numrows,numcols,numplates);
cy5mat = yfpmat;

allyfp = cell(numrows,numcols,numplates);
allcy5 = allyfp;

for plate = 1:numplates
    
%     if plates(plate) == 1
%         rows = 3:4; numrows = numel(rows);
%         cols = 4:6; numcols = numel(cols);
%     else
%         rows = 3:6; numrows = numel(rows);
%         cols = 4:9; numcols = numel(cols);
%     end
    
    
%     sites = 1:4; numsites = numel(sites);
    
    for col = 1:numcols

        for row = 1:numrows
            yfpwell = [];
            cy5well = [];
            
            for site = 1:numsites
                
%                 if rows(row) == 2 && cols(col) ==6 && sites(site)==1 && plates(plate) == 6
%                     continue
%                 end
                
                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',num2str(plate),'_',shot,'.mat']); 
                yfp = fixdata(:,6);
                cy5 = fixdata(:,7);
                
                
                yfpwell = [yfpwell;yfp];
                cy5well = [cy5well;cy5];
    

            end
            
%             figure(plate)
%             subtightplot(numrows,numcols,posmat(row,col),[0.03 0.03])
%             histogram(yfpwell,'EdgeColor','none');
%             [f,x] = ksdensity(yfpwell);
%             plot(x,f,'LineWidth',2);
%             set(gca,'XLim',[0 250])
            
%             figure(plate+numplates)
%             subtightplot(numrows,numcols,posmat(row,col),[0.03 0.03])
%             histogram(cy5well,'EdgeColor','none');
%             [f,x] = ksdensity(cy5well);
%             plot(x,f,'LineWidth',2);
%             if row == 1 || row == 3
%                 set(gca,'XLim',[-50 500])
%             else
%                 set(gca,'XLim',[-100 1300])
%             end
            
            badyfp = yfpwell>prctile(yfpwell,95) | yfpwell<prctile(yfpwell,5);
            badcy5 = cy5well>prctile(cy5well,95) | cy5well<prctile(cy5well,5);
            badcells = badyfp | badcy5;
            
            yfpwell(badcells) = [];
            cy5well(badcells) = [];
            
            yfpmat(row,col,plate) = mean(yfpwell);
            cy5mat(row,col,plate) = mean(cy5well);
            
            allyfp{row,col,plate} = yfpwell;
            allcy5{row,col,plate} = cy5well;

            

        end

    end



end

% yfpmat(:,:,1) = repmat(yfpmat(1:2,1:3,1),2,2);
% cy5mat(:,:,1) = repmat(cy5mat(1:2,1:3,1),2,2);

%%
GRpulse = zeros(numplates,3); GRcont = GRpulse;
cebpbpulse = zeros(numplates,3); cebpbcont = cebpbpulse;
cebpapulse = zeros(numplates,3); cebpacont = cebpapulse;
ppargpulse = zeros(numplates,3); ppargcont = ppargpulse;

for p = 1:numplates
    GRpulse(p,:) = cy5mat(3,4:6,p); GRcont(p,:) = cy5mat(3,1:3,p);
    cebpbpulse(p,:) = yfpmat(3,4:6,p); cebpbcont(p,:) = yfpmat(3,1:3,p);
    cebpapulse(p,:) = cy5mat(4,4:6,p); cebpacont(p,:) = cy5mat(4,1:3,p);
    ppargpulse(p,:) = yfpmat(4,4:6,p); ppargcont(p,:) = yfpmat(4,1:3,p);
end

% for p = 1:numplates
%     GRpulse(p,:) = cy5mat(3,4:6,p)./cy5mat(1,4:6,p); GRcont(p,:) = cy5mat(3,1:3,p)./cy5mat(1,1:3,p);
%     cebpbpulse(p,:) = yfpmat(3,4:6,p)./yfpmat(1,4:6,p); cebpbcont(p,:) = yfpmat(3,1:3,p)./yfpmat(1,1:3,p);
%     cebpapulse(p,:) = cy5mat(4,4:6,p)./cy5mat(2,4:6,p); cebpacont(p,:) = cy5mat(4,1:3,p)./cy5mat(2,1:3,p);
%     ppargpulse(p,:) = yfpmat(4,4:6,p)./yfpmat(2,4:6,p); ppargcont(p,:) = yfpmat(4,1:3,p)./yfpmat(2,1:3,p);
% end

% GRcont(1:2,:) = [];
% cebpbcont(1:2,:) = [];
% cebpacont(1:2,:) = [];
% ppargcont(1:2,:) = [];

figure
set(gcf,'color','white')
indmat = [1 1 1; 2 2 2; 3 3 3; 4 4 4; 5 5 5;6 6 6];
indmat2 = [3 3 3; 4 4 4; 5 5 5;6 6 6];
subtightplot(2,2,1,[0.1 0.05],[.05 .05], [.05 .02])
% scatter(indmat(:),GRpulse(:),100,'filled','MarkerEdgeColor','k')
barwitherr(std(GRpulse,0,2),mean(GRpulse,2))
% barwitherr(std(GRcont,0,2),mean(GRcont,2))
hold on
% scatter(indmat(:),GRcont(:),100,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980])
set(gca,'XLim', [0 6],'XTickLabel',{'0h','12h','24h','36h','48h',''},'FontSize',20)
ylabel('Mean Nuclear Intensity')
title('GR')

subtightplot(2,2,2,[0.1 0.05],[.05 .05], [.05 .02])
% scatter(indmat(:),cebpbpulse(:),100,'filled','MarkerEdgeColor','k')
barwitherr(std(cebpbpulse,0,2)./sqrt(3),mean(cebpbpulse,2))
errorbar(mean(cebpbpulse,2),std(cebpbpulse,0,2)./sqrt(3),'LineWidth',1,'Color','k')
% barwitherr(std(cebpbcont,0,2),mean(cebpbcont,2))
hold on
% scatter(indmat(:),cebpbcont(:),100,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980])
set(gca,'XLim', [0 6],'XTickLabel',{'','0h','12h','24h','36h','48h',''},'FontSize',14)
ylabel('Mean Nuclear Intensity')
xlabel('Time (hours)')
title('CEBPB')

subtightplot(2,2,3,[0.1 0.05],[.05 .05], [.05 .02])
% scatter(indmat(:),cebpapulse(:),100,'filled','MarkerEdgeColor','k')
barwitherr(std(cebpapulse,0,2),mean(cebpapulse,2))
% barwitherr(std(cebpacont,0,2),mean(cebpacont,2))
hold on
% scatter(indmat(:),cebpacont(:),100,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980])
set(gca,'XLim', [0 6],'XTickLabel',{'0h','12h','24h','36h','48h',''},'FontSize',20)
ylabel('Mean Nuclear Intensity')
title('CEBPA')

subtightplot(2,2,4,[0.1 0.05],[.05 .05], [.05 .02])
% scatter(indmat(:),ppargpulse(:),100,'filled','MarkerEdgeColor','k')
barwitherr(std(ppargpulse,0,2),mean(ppargpulse,2))
% barwitherr(std(ppargcont,0,2),mean(ppargcont,2))
hold on
% scatter(indmat(:),ppargcont(:),100,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980])
set(gca,'XLim', [0 6],'XTickLabel',{'0h','12h','24h','36h','48h',''},'FontSize',20)
ylabel('Mean Nuclear Intensity')
title('PPARG')
