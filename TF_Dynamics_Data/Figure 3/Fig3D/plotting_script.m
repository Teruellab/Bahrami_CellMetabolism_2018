
imagepath='Z:\Michael\';

experimentpath='Live Cell Experiments\20151203-pparg-stim-duration\';
datadir=([imagepath,experimentpath,'Data\']);

load([datadir,'48h_CombinedWithLinking.mat'])
badcellIDs = csvread([imagepath,experimentpath,'\badcellIDs.csv']);

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
endframe = 618;
frames=startframe:endframe;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1= 0; ymax1 = 1000;
ymin2= -2000; ymax2= 2000; %ymax2=1000;

numgated=length(alltracestats(:,1));
allIDs = 1:numgated;

for i = 1:numgated
    lastind = find(~isnan(alltraces1(i,:)),1,'last');
    lastvalue(i) = alltraces1(i,lastind);
end
startmedian = zeros(numgated,1);
endmedian = startmedian;
endvalue = startmedian;
sumthresh = startmedian;
maxdiff = startmedian;
maxloss = startmedian;
midval = startmedian;
smoothedmat = zeros(numgated,endframe);
for i = 1:numgated
   tracestart = alltracestats(i,1);
   traceend = alltracestats(i,2);
   smoothedtrace = smoothignorenans(gatetraces(i,:),1);
   smoothedtrace2 = smoothignorenans(alltraces1(i,:),3);
   maxdiff(i) = max(diff(alltraces1(i,tracestart:traceend)));
   maxloss(i) = min(diff(alltraces1(i,tracestart:traceend))./alltraces1(i,tracestart:traceend-1));

   endmedian(i) = mean(smoothedtrace(traceend-3:traceend));
   midval(i) = mean(smoothedtrace(287-3:287));
   endvalue(i) = smoothedtrace(traceend);
   sortedvalues = sort(smoothedtrace(tracestart:traceend));
   sumthresh(i) = sum(bwareaopen(smoothedtrace>128,12));
   smoothedmat(i,:) = smoothedtrace;
end

foldchange = endmedian./startmedian;
sumthresh = logical(sumthresh>0);
highcells = endmedian > 2^7.4;
diffcells = endmedian> 2^7.4;
badfoldchange = foldchange<0 | foldchange>(prctile(foldchange,99));
lowcells = endmedian<200&sumthresh;
goodstart = startmedian<150;
noisy = maxloss<-0.6 | maxdiff>120;


disp(100*sum(highcells)/numel(highcells))
startpoint = 72;
endpoint = 602;
gatestartpoint = alltracestats(:,1) < startpoint;
gateendpoint = alltracestats(:,2)> endpoint;
goodlength = endpoint-startpoint;
toohigh = max(alltraces1,[],2)<1000;
nantotal = sum(isnan(alltraces1),2);
badnan = nantotal > (endframe-goodlength);

[~,~,badtraceind] = intersect(badcellIDs,allIDs);
gatein = gatestartpoint & gateendpoint;
gatein(badtraceind) = 0;

hightemptraces = alltraces1(highcells,:);
lowtemptraces = alltraces1(~highcells,:);
highmean = nanmedian(hightemptraces);
lowmean = nanmedian(lowtemptraces); 

alltraces1 = alltraces1(gatein,:) ;
alltraces2 = alltraces2(gatein,:) ;
alltraces3 = alltraces3(gatein,:) ;
allwellID = allwellID(gatein,:);
allIDs = allIDs(gatein);
alltracestats = alltracestats(gatein,:);
allmarkedmitosis = allmarkedmitosis(gatein);
highcells = highcells(gatein);
diffcells = diffcells(gatein);
smoothedmat = smoothedmat(gatein,:);
numgated = length(allIDs);

% plot 48 versus 96 hour
day2values = alltraces1(:,300);
day4values = alltraces1(:,600);
bad2 = day2values<8;
bad4 = day4values<8;
badcells = bad2|bad4|isnan(day2values)|isnan(day4values);
day2values(badcells) = [];
day4values(badcells) = [];
figure
dscatter_diffcells(log2(day2values),log2(day4values),7.4,3:.3:11)

