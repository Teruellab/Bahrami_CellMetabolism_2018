%% data extraction

imagepath = 'Y:\michael\';
experimentpath='TF dynamics data\Figure S5\FigS5I\'; 
datadir = [imagepath,experimentpath,'Data\'];


rows = 2:4; numrows = numel(rows);
cols = [10]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1:6]; numplates = numel(plates);

numwells = numrows*numcols;
postmat = reshape(1:numwells,numcols,numrows)';
mchyplate = zeros(numrows,numcols,numplates);

for plate = 1:numplates
% figure(2+plate)
    for row = 1:numrows
        for col = 1:numcols
            
            
            mchywell = [];
            dnawell = [];
            for site = 1:numsites

                shot = [num2str(plates(plate)),'_',num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',shot,'.mat'])
                mchy = fixdata(:,7);

                mchywell = [mchywell;mchy];
                dnawell = [dnawell;fixdata(:,4)];

            end
            
            badmchy = mchywell<0 | mchywell>600;
            badcells = badmchy;
            mchywell(badcells) = [];
            
            mchyplate(row,col,plate) = nanmean(mchywell);
            
        end
            
%             hold off

    end
end
%% Plot by column
time = [0 .5 1 1.5 3];
timemat = repmat(time,3,1);
timerange = time(1):.1:time(end);
interval = [1:5];
for ind = 1 
    
    mchytemp = squeeze(mchyplate(:,ind,interval));
    mchytempnorm = mchytemp./mean(mchytemp(:,1));
    mchymean = mean(mchytemp);
    mchystd = std(mchytemp);
    mchystd = mchystd./mchymean(1);
    mchymean = mchymean./mchymean(1);
    
    fitnorm = mchytempnorm;

    modelfit = fit(timemat(:),fitnorm(:),'exp1');    
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
    patch(X,F,'r','FaceAlpha',0.2,'EdgeColor','none')
    plot(timemat',fitnorm','LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
    modelhandle.LineWidth = 1.5;
    set(gca,'YLim',[0 1.2])
    text(0.3, 0.2, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'])
    
end


xlabel 'Time (hrs)'
ylabel 'Normalized CEBPB Intensity'
