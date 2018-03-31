%% data extraction

imagepath = 'Y:\michael\';
experimentpath='Fixed Cell Experiments\2016\20160720-CRISPR-Degradation\'; 
datadir = [imagepath,experimentpath,'Data\'];


rows = 2:4; numrows = numel(rows);
cols = [1:3]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1:5]; numplates = numel(plates);

numwells = numrows*numcols;
postmat = reshape(1:numwells,numcols,numrows)';
mchyplate = zeros(numrows,numcols,numplates);
yfpplate = mchyplate;
figure(1)
set(gcf,'color','white')
figure(2)
set(gcf,'color','white')


for plate = 1:numplates
% figure(2+plate)
    for row = 1:numrows
        for col = 1:numcols
            
            
            mchywell = [];
            yfpwell = [];
            dnawell = [];
            for site = 1:numsites

                shot = [num2str(plates(plate)),'_',num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',shot,'.mat'])
                mchy = fixdata(:,7); %antibody staining signal
                yfp = fixdata(:,6); % CRISPR signal, if relevant

                mchywell = [mchywell;mchy];
                yfpwell = [yfpwell;yfp];
                dnawell = [dnawell;fixdata(:,4)];

            end
            
            % Eliminate ridiculously large values
            
            badmchy = mchywell<0 | mchywell>600;
            badyfp = yfpwell<0 | yfpwell > 150;
            badcells = badmchy | badyfp;
            mchywell(badcells) = [];
            yfpwell(badcells) = [];
            
            mchyplate(row,col,plate) = nanmean(mchywell);
            yfpplate(row,col,plate) = nanmean(yfpwell);
            
%             figure(1)
%             hold on
%             subtightplot(numrows,numcols,postmat(row,col)) 
            

                
%                 histogram(yfpwell,0:2.5:150,'EdgeColor','none','FaceAlpha',0.4)
%                 [f,x] = ksdensity(yfpwell);
%                 plot(x,f,'LineWidth',1);
%                 axis square
%                 set(gca,'XLim',[0 xlimmat(col)])
%              hold off
%              figure(2)   
%              hold on
%              subtightplot(numrows,numcols,postmat(row,col)) 
%                 histogram(mchywell,0:12:300,'EdgeColor','none','FaceAlpha',0.4)
                
%                 [f,x] = ksdensity(mchywell);
%                 plot(x,f,'LineWidth',1);
%                 axis square
%                 set(gca,'XLim',[0 xlimmat(col)])
%             figure(2+plate)
%             subtightplot(numrows,numcols,postmat(row,col)) 
%             dscatter(mchywell,yfpwell)
        end
            
%             hold off

    end
end
%% Plot by column
time = [0 .5 1 1.5 3];
timemat = repmat(time,3,1);
timerange = time(1):.1:time(end);
interval = [1:5];
for ind = 1 %ind 1 = CHX, ind 2 = MG132, ind 3 = DMSO
    yfptemp = squeeze(yfpplate(:,ind,interval));
    yfptempnorm = yfptemp./mean(yfptemp(:,1));
    yfpmean = mean(yfptemp); 
    yfpstd = std(yfptemp); 
    yfpstd = yfpstd./yfpmean(1);
    yfpmean = yfpmean./yfpmean(1);
    
    mchytemp = squeeze(mchyplate(:,ind,interval));
    mchytempnorm = mchytemp./mean(mchytemp(:,1));
    mchymean = mean(mchytemp);
    mchystd = std(mchytemp);
    mchystd = mchystd./mchymean(1);
    mchymean = mchymean./mchymean(1);
    
%     figure(5)
%     subtightplot(4,3,1,[0.03 0.03])
%     hold on
%     errorbar(time,yfpmean,yfpstd,'LineWidth',2)
%     set(gca,'YLim',[0 2])
%     hold off
%     figure(6)
%     subtightplot(4,3,1,[0.03 0.03])
%     hold on
%     errorbar(time,mchymean,mchystd,'LineWidth',2)
%     set(gca,'YLim',[0 2])
%     hold off
    
    modelfit = fit(timemat(:),mchytempnorm(:),'exp1');    
    confidence = confint(modelfit);
    lowerint = confidence(1,:);
    upperint = confidence(2,:);
    lb  = lowerint(1).*exp(lowerint(2).*timerange);
    ub = upperint(1).*exp(upperint(2).*timerange);
    halflife = log(2)/-modelfit.b;
    llife = log(2)/-lowerint(2);
    ulife = log(2)/-upperint(2);
    X = [timerange, fliplr(timerange)];
    F = [ub, fliplr(lb)];
    figure
    modelhandle = plot(modelfit);
    hold on 
%     plot(timerange,lb,'LineStyle','--','Color','k','LineWidth',1)
%     plot(timerange,ub,'LineStyle','--','Color','k','LineWidth',1)
    patch(X,F,'r','FaceAlpha',0.2,'EdgeColor','none')
    plot(timemat',mchytempnorm','LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
    modelhandle.LineWidth = 1.5;
    set(gca,'YLim',[0 1.2])
    text(0.3, 0.2, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'])
    
end



