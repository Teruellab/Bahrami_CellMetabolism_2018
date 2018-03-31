%% data extraction

imagepath = 'Y:\michael\';
experimentpath='TF dynamics data\protein_48hr_decay\'; 
datadir = [imagepath,experimentpath,'Data\'];

%% Forming plate 1 (time 0) matrix
rows = 2:3; numrows = numel(rows);
cols = [2,3,7:11]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1]; numplates = numel(plates);

numwells = numrows*numcols;
postmat = reshape(1:numwells,numcols,numrows)';
mchyplate = zeros(numrows,numcols,numplates);
mchycytoplate = mchyplate;

for plate = 1:numplates
    for row = 1:numrows
        for col = 1:numcols
            
            if plates(plate) == 1 && rows(row)==3 && cols(col) > 4
                continue
            end
            
            mchywell = [];
            mchycytowell = [];

            for site = 1:numsites

                shot = [num2str(plates(plate)),'_',num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',shot,'.mat'])
                mchy = fixdata(:,6);
                mchycyto = fixdata(:,8);

                mchywell = [mchywell;mchy];
                mchycytowell = [mchycytowell;mchycyto];

            end
            
            % Eliminate ridiculously large values
            
            badmchy = mchywell > prctile(mchywell,99)*3;
            badmchycyto = mchycytowell > prctile(mchycytowell,99)*3 | isnan(mchycytowell);
            badcells = badmchy | badmchycyto;
            mchywell(badcells) = [];
            mchycytowell(badcells) = [];
            
            ctonwell(row,col,plate) = mean(mchycytowell./mchywell);
            mchyplate(row,col,plate) = nanmean(mchywell);
            mchycytoplate(row,col,plate) = nanmean(mchycytowell);

        end
    end

end
pparg0 = mchyplate(1,1:2)';
cebpa0 = mchyplate(1,3:5)';
fabp4_ab0 = mchycytoplate(1,6:7)';
fabp4_rd0 = mchycytoplate(2,1:2)';

controlmat = [[pparg0; nan],cebpa0,[fabp4_ab0;nan],[fabp4_rd0;nan]];
%% Forming other time points (0.5hr to 3hr)
clearvars -except controlmat imagepath experimentpath datadir
rows = 2:4; numrows = numel(rows);
cols = [2,4,5,6]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [2:5]; numplates = numel(plates);

numwells = numrows*numcols;
postmat = reshape(1:numwells,numcols,numrows)';
mchyplate = zeros(numrows,numcols,numplates);
mchycytoplate = mchyplate;

for plate = 1:numplates
    for row = 1:numrows
        for col = 1:numcols
            
            if plates(plate) == 1 && rows(row)==3 && cols(col) > 4
                continue
            end
            
            mchywell = [];
            mchycytowell = [];

            for site = 1:numsites

                shot = [num2str(plates(plate)),'_',num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',shot,'.mat'])
                mchy = fixdata(:,6);
                mchycyto = fixdata(:,8);

                mchywell = [mchywell;mchy];
                mchycytowell = [mchycytowell;mchycyto];

            end
            
            % Eliminate ridiculously large values
            
            badmchy = mchywell > prctile(mchywell,99)*3;
            badmchycyto = mchycytowell > prctile(mchycytowell,99)*3 | isnan(mchycytowell);
            badcells = badmchy | badmchycyto;
            mchywell(badcells) = [];
            mchycytowell(badcells) = [];
            
            ctonwell(row,col,plate) = mean(mchycytowell./mchywell);
            mchyplate(row,col,plate) = nanmean(mchywell);
            mchycytoplate(row,col,plate) = nanmean(mchycytowell);

        end
    end

end

%% Plotting Normalized Values with etimated confidence bounds
condind = [1:5];
dmsoind = [6:10];
numcond = numel(condind);
% time = [0 0.5 1 1.5 3 5 7 9];
time = [0 0.5 1 1.5 3];
timerange = 0:.01:3;
timemat = repmat(time,3,1);
protein = {'PPARG','CEBPA','FABP4-abcam','FABP4-RandD'};
for c = [1,2,3,4]
    dmii = condind(c);
    dmsoi = dmsoind(c);
    interval = [1:4];
    if c<3
        chx = [controlmat(:,c),squeeze(mchyplate(:,dmii,interval))]./nanmean(controlmat(:,c));
    else
        chx = [controlmat(:,c),squeeze(mchycytoplate(:,dmii,interval))]./nanmean(controlmat(:,c));
    end
    timevec = timemat(:);
    chxvec = chx(:);
    badpoint = isnan(chxvec);
    timevec(badpoint) = [];
    chxvec(badpoint) = [];
    modelfit = fit(timevec,chxvec,'exp1');
    confidence = confint(modelfit);
    lowerint = confidence(1,:);
    upperint = confidence(2,:);
    lb  = lowerint(1).*exp(lowerint(2).*timerange);
    ub = upperint(1).*exp(upperint(2).*timerange);
%     pinterval = predint(modelfit,timerange,.95,'functional','on');
%     lb = pinterval(:,1)';
%     ub = pinterval(:,2)';
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
    plot(timevec,chxvec,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
    modelhandle.LineWidth = 1.5;
    text(0.3, 0.2, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'])
    legend off
    title(protein{c})
    ylabel 'Normalized Value'
    xlabel 'Time (hours)'
    set(gca,'YLim',[0 1.2],'YTick',[0:0.2:1],'XLim',[-0.1 3.1],'XTick',[0:.5:3],'FontName','Arial','FontSize',20)
    
end


