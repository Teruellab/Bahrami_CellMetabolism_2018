%% data extraction

imagepath = 'Y:\michael\';
experimentpath='TF dynamics data\Figure S6\FigS6I\'; 
datadir = [imagepath,experimentpath,'Data\'];


rows = 2:4; numrows = numel(rows);
cols = [10:12]; numcols = numel(cols);
sites = 1:4; numsites = numel(sites);
plates = [1:6]; numplates = numel(plates);

numwells = numrows*numcols;
postmat = reshape(1:numwells,numcols,numrows)';
cy5plate = zeros(numrows,numcols,numplates);


for plate = 1:numplates

    for row = 1:numrows
        for col = 1:numcols
            
            
            cy5well = [];
            dnawell = [];
            for site = 1:numsites

                shot = [num2str(plates(plate)),'_',num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                datafile = [datadir,'fixdata_Plate_',shot,'.mat'];
                if ~exist(datafile,'file')
                    continue
                end
                load(datafile)
                cy5 = fixdata(:,7);

                cy5well = [cy5well;cy5];
                dnawell = [dnawell;fixdata(:,4)];

            end
            
            % Eliminate ridiculously large values
            
            badcy5 = cy5well<0 | cy5well>prctile(cy5well,99)*3;
            baddna = dnawell<80000 | dnawell>800000;
            badcells = badcy5 | baddna;
          
            cy5plate(row,col,plate) = nanmean(cy5well);           

        end
            


    end
end

iimat = 1:3;
cy5normvec = cy5plate(1,:,1);
for ii = 1
     cy5plate(:,iimat(ii):iimat(ii)+2,1) = repmat(cy5normvec(iimat(ii):iimat(ii)+2)',1,3);
end

%% Plot by column
time = [0 0.5 1 1.5 3 ];
timerange = time(1):.1:time(end);
timemat = repmat(time,3,1);
indlabel = {'CHX-48h','MG132-48h','DMSO-48h'};
for ind = 1
   
    cy5temp = squeeze(cy5plate(:,ind,1:5));
    cy5tempnorm = cy5temp./nanmean(cy5temp(:,1));
    cy5mean = nanmean(cy5temp);
    cy5std = std(cy5temp);
    cy5std = cy5std./cy5mean(1);
    cy5mean = cy5mean./cy5mean(1);
    
    modelfit = fit(timemat(:),cy5tempnorm(:),'exp1');    
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
    plot(timemat',cy5tempnorm','LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
    modelhandle.LineWidth = 1.5;
    set(gca,'YLim',[0 1.2])
    text(0.3, 0.2, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'])
    title(indlabel{ind})
end
xlabel 'Time (hrs)'
ylabel 'Normalized PPARG Intensity'

