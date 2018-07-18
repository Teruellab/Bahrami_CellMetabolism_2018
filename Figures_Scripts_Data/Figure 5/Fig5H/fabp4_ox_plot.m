
imagepath='Y:\michael\';

experimentpath='TF dynamics data\Fig5H\';
datadir=([imagepath,experimentpath,'Data\']);
siName = {'mchy','fabp4-mchy'};

numsi = numel(siName);

r = [231,76,60]./255;
g = [77,175,74]./255;
b = [8,81,156]./255;
gy = [90,90,90]./255;
gd = [253,192,16]./255;
colorchoices = [r;b];

fignum = 60;

for sin = 1:numsi

load([datadir,siName{sin},'_MackTrackCombined.mat'])


gatetraces = alltraces1;

numgated=length(alltracestats(:,1));
allIDs = 1:numgated;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr = 8;
framesperhr = 5;
drugspike = 0; %114
spikeTime = [338]./framesperhr; %,142,205,284,352,424,493
numspikes = numel(spikeTime);
startframe = 1;
endframe = 596;
frames=startframe:endframe;

for i = 1:numgated
    lastind = find(~isnan(alltraces1(i,:)),1,'last');
    lastvalue(i) = alltraces1(i,lastind);
end
% histogram(lastvalue,20)
startmedian = zeros(numgated,1);
endmedian = startmedian;
endvalue = startmedian;
sumthresh = startmedian;
maxdiff = startmedian;
midval = startmedian;
for i = 1:numgated
   tracestart = alltracestats(i,1);
   traceend = alltracestats(i,2);
   smoothedtrace = smoothignorenans(gatetraces(i,:),1);
   smoothedtrace2 = smoothignorenans(alltraces1(i,:),3);
   maxdiff(i) = min((diff(alltraces1(i,tracestart:traceend)))./alltraces1(i,tracestart:traceend-1));
%    startmedian(i) = nanmedian(smoothedtrace(lastmitosisframe(numgated):113)); 
   startmedian(i) = nanmedian(smoothedtrace(1:9)); 
   endmedian(i) = median(smoothedtrace(traceend-6:traceend));
   midval(i) = nanmean(smoothedtrace(285:288));
   endvalue(i) = smoothedtrace(traceend);
   sortedvalues = sort(smoothedtrace(tracestart:traceend));
   sumthresh(i) = sum(smoothedtrace>2^7.3);
%    sumthresh(i) = sum(gatetraces(i,traceend-110:traceend)>200);
end

foldchange = endmedian./startmedian;
% sumthresh = logical(sumthresh<50);
highcells = endmedian > 2^9;
% highcells = midval>2^7;
diffcells = endmedian>2^7.4;
badfoldchange = foldchange<-10000 | foldchange>(prctile(foldchange,99));
lowcells = endmedian<200&sumthresh;
goodstart = startmedian<150;
noisy = maxdiff<-0.55;

disp((sum(highcells)/length(highcells))*100)
disp((sum(diffcells)/length(diffcells))*100)

startdiff = nanmedian(alltraces1(highcells,52)) - nanmedian(alltraces1(~highcells,52));

mitosiscount = cellfun(@(x) numel(x),allmarkedmitosis)-1;
startpoint = 50;
endpoint = 550;
gatestartpoint = alltracestats(:,1) < startpoint;
gateendpoint = alltracestats(:,2)> endpoint;
goodlength = endpoint-startpoint;
toohigh = max(alltraces1,[],2)<1000;
nantotal = sum(isnan(alltraces1),2);
badnan = nantotal > (endframe-goodlength);
gatein = gatestartpoint & gateendpoint;


alltraces1 = alltraces1(gatein,:) ;
alltraces2 = alltraces2(gatein,:) ;
alltraces3 = alltraces3(gatein,:) ;
allwellID = allwellID(gatein,:);
allIDs = allIDs(gatein);
alltracestats = alltracestats(gatein,:);
allmarkedmitosis = allmarkedmitosis(gatein);
highcells = highcells(gatein);
nantotal = nantotal(gatein);
numgated = length(allIDs);

% plot median + percentile bounds
tracematrix = alltraces1(highcells,:);

colorchoice = colorchoices(sin,:);
xtime=(((1:endframe))/framesperhr)-25;
middleTrace = nanmedian(tracematrix,1);
bottomBound = prctile(tracematrix,25);
topBound = prctile(tracematrix,75);

firstgoodval = find(~isnan(middleTrace),1,'first');
lastgoodval =  find(~isnan(middleTrace),1,'last');

middleTrace = middleTrace(firstgoodval:lastgoodval);
bottomBound = bottomBound(firstgoodval:lastgoodval);
topBound = topBound(firstgoodval:lastgoodval);
xtime = xtime(firstgoodval:lastgoodval);

X = [xtime, fliplr(xtime)];
F = [topBound, fliplr(bottomBound)];

figure(fignum)
set(gcf,'Renderer','painter')
hold on

fill(X,F,colorchoice,'FaceAlpha',0.4,'EdgeColor','none');

plot(xtime,middleTrace,'LineWidth',3,'Color',colorchoice);
set(gca,'Layer','Top');

hold off

end
set(gca,'XLim',[0 47],'Xtick',0:20:40,'YLim',[0 2.1e3]);