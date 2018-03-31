clear all 
close all

imagedir = 'Y:\michael\';
experimentdir = [imagedir,'TF dynamics data\Fig5B\'];
datadir = [experimentdir,'Data\'];
% declare row, col, and sites to analyze
rows = 3; numrows = numel(rows);
cols = [3:5]; numcols = numel(cols);
sites = 1:9; numsites = numel(sites);
plates = [1:6]; numplates = numel(plates);
posmat = 1:numrows*numcols;
posmat = reshape(posmat,numcols,numrows)';

% figure(1)
% fishcell = cell(numplates);
for plate = 1:numplates
%     figure
    for col = 1:numcols


        for row = 1:numrows

            fish = [];
            tempmat = zeros(3,3);

            for site = 1:numsites

                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                load([datadir,'fixdata_Plate_',num2str(plates(plate)),'_',shot,'.mat']); 
                fish = [fish; IFdata(:,10)];
                [i,j] = ind2sub([3,3],site);
                tempmat(i,j) = nanmean(IFdata(:,10));
                
            end
%             qcmat(row,col) = sum(fish>0)*100/numel(fish);
            figure(100)
            hold on
            subtightplot(numrows,numcols,posmat(row,col),[0.06,0.01])
            hold on
            histogram(fish,'EdgeColor','none','FaceAlpha',0.3,'normalization','probability')
%             set(gca,'XLim',[0 500])
            hold off

            fishmat(row,col,plate) = nanmean(fish);
            fishmatstd(row,col,plate) = nanmean(fish);  
%             fishcell(plate) = [fishcell{plate};fish];
%             figure(plate)
%             subtightplot(numrows,numcols,posmat(row,col),[0.06,0.01])
%             imagesc(tempmat)
%             axis square
%             colorbar
%             caxis([150 300])

        end

    end

end

%% plotting degradation with confidence bounds 2016-11-21
time = [0 1 3 6 9 12];
interval = 1:6;
timemat = repmat(time,3,1);
timevec = timemat(:);
timerange = 0:.01:time(end);
gene = {'FABP4'};
for c = 1
    fishvals = squeeze(fishmat(c,:,interval));
    fishvec = fishvals(:);
    modelfit = fit(timevec,fishvec,'exp1');
    confidence = confint(modelfit);
    lowerint = confidence(1,:);
    upperint = confidence(2,:);
    lb  = lowerint(1).*exp(lowerint(2).*timerange);
    ub = upperint(1).*exp(upperint(2).*timerange);
    halflife = round(10*log(2)/-modelfit.b)/10;
    llife = round(10*log(2)/-lowerint(2))/10;
    ulife = round(10*log(2)/-upperint(2))/10;
    X = [timerange, fliplr(timerange)];
    F = [ub, fliplr(lb)];
    figure(2)
    modelhandle = plot(modelfit);
    modelhandle.LineWidth = 1.5;
    hold on
    patch(X,F,'r','FaceAlpha',0.2,'EdgeColor','none')
    plot(timevec,fishvec,'LineStyle','none','Marker','o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','none')
    title(gene{c})
    xlabel 'Time (hours)'
    ylabel 'mRNA puncta'
    set(gca,'YLim',[0 150],'XLim',[-0.1 12.1],'XTick',[0:1:12],'FontName','Arial','FontSize',20)
    text(0,5, ['Half-life = ',num2str(halflife),' (',num2str(llife),', ',num2str(ulife),')'],'FontName','Arial','FontSize',16)
    
end

