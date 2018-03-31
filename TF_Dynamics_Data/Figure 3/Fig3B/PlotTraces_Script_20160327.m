
imagepath='Z:\Michael\';

experimentpath='Live Cell Experiments\20160327-alphapparg-pulse-ZB\';
datadir=([imagepath,experimentpath,'Data\']);

load([datadir,'48hr_CombinedWithLinking.mat'])
gatetraces = alltraces1;

%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 1;    %0:one trace per plotn 1:all traces in one plot
selectmode = 0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2 = 0;
plotsignal3 = 0;
immunoframe=0;
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
removequiescent = 0;
%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr = 8;
framesperhr = 6;
drugspike = 1; %114
spikeTime = 286./framesperhr;%[69,141,217,286,357,430,509]./framesperhr; 
numspikes = numel(spikeTime);
startframe = 1;
endframe = 577;
frames=startframe:endframe;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1= 0; ymax1 = 1000;
ymin2= -2000; ymax2= 2000; %ymax2=1000;

numgated=length(alltracestats(:,1));
allIDs = 1:numgated;
% remove = zeros(numgated,1);
% remove(exclude) = 1;
% keep = remove==0;

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
maxloss = startmedian;
smoothedmat = zeros(numgated,endframe);
for i = 1:numgated
   tracestart = alltracestats(i,1);
   traceend = alltracestats(i,2);
   smoothedtrace = smoothignorenans(gatetraces(i,:),1);
   smoothedtrace2 = smoothignorenans_butter(alltraces1(i,:),3,0.1);
   maxdiff(i) = max(diff(alltraces1(i,tracestart:traceend)));
   maxloss(i) = min(diff(alltraces1(i,tracestart:traceend))./alltraces1(i,tracestart:traceend-1));
   endmedian(i) = mean(smoothedtrace(traceend-6:traceend));
   endvalue(i) = smoothedtrace(traceend);
   sortedvalues = sort(smoothedtrace(tracestart:traceend));
   sumthresh(i) = sum(bwareaopen(smoothedtrace>181,3));
   smoothedmat(i,:) = smoothedtrace2;
end

highcells = endmedian>181;
lowcells = endmedian<200&sumsig;
goodstart = startmedian<150;
noisy = maxloss<-0.6 | maxdiff>120;

disp(100*sum(highcells)/numel(highcells))

% startdiff = nanmedian(alltraces1(highcells,52)) - nanmedian(alltraces1(~highcells,52));

% gatemitosisframe = logical(lastmitosisframe<300);
startpoint = 20;
endpoint = 570;
gatestartpoint = alltracestats(:,1) < startpoint;
gateendpoint = alltracestats(:,2)> endpoint;

gatein = gatestartpoint & gateendpoint;


alltraces1 = alltraces1(gatein,:) ;
alltraces2 = alltraces2(gatein,:) ;
alltraces3 = alltraces3(gatein,:) ;
allwellID = allwellID(gatein,:);
allIDs = allIDs(gatein);
alltracestats = alltracestats(gatein,:);
allmarkedmitosis = allmarkedmitosis(gatein);
highcells = highcells(gatein);
smoothedmat = smoothedmat(gatein,:);
numgated = length(allIDs);

diffmean = nanmean(alltraces1(highcells,:));
undiffmean = nanmean(alltraces1(~highcells,:));

if numgated > 50
    indexvector = 1:numgated;
    numgated = 50;
    randomized = sort(datasample(indexvector,numgated,'Replace',false));
    alltraces1 = alltraces1(randomized,:) ;
    alltraces2 = alltraces2(randomized,:) ;
    alltraces3 = alltraces3(randomized,:) ;
    alltracestats = alltracestats(randomized,:);
    allwellID = allwellID(randomized,:);
    allIDs = allIDs(randomized);
    allmarkedmitosis = allmarkedmitosis(randomized);
    highcells = highcells(randomized);
    smoothedmat = smoothedmat(randomized,:);
    lowcells = lowcells(randomized);
end

r = [[231,76,60]./255]; b = [[8,81,156]./255]; 

figure,plot(frames./framesperhr,smoothedmat(highcells,:)','Color',[r, 0.4],'LineWidth',1)
hold on
plot(frames./framesperhr,diffmean,'Color',r,'LineWidth',4)
plot(frames./framesperhr,smoothedmat(~highcells,:)','Color',[b .4],'LineWidth',1)
plot(frames./framesperhr,undiffmean,'Color',b,'LineWidth',4)

xlabel 'Time (hrs)'
ylabel 'PPARG'