
imagepath='Y:\michael\Live Cell Experiments\';

experimentpath='20150813-CEBPb-live\';
datadir=([imagepath,experimentpath,'Data\']);
conditionname = 'C13 - 4x12hr pulse';
tic     
load([datadir,conditionname,'_NoLinking.mat'])
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alltraces1: CEBPB integrated intensity
%alltraces2: H2B mean intensity 
%alltraces3: H2B integrated intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startpoint = 20;
endpoint = 529;
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
lastmitosisframe = lastmitosisframe(gatein);
numgated = length(allIDs);

sampleNumber = 10;
if numgated > sampleNumber
    indexvector = 1:numgated;
    numgated = sampleNumber;
    randomized = sort(datasample(indexvector,numgated,'Replace',false));
    alltraces1 = alltraces1(randomized,:) ;
    alltraces2 = alltraces2(randomized,:) ;
    alltraces3 = alltraces3(randomized,:) ;
    alltracestats = alltracestats(randomized,:);
    allwellID = allwellID(randomized,:);
    allmotherstats = allmotherstats(randomized);
    allIDs = allIDs(randomized);
    allmarkedmitosis = allmarkedmitosis(randomized);
end


%%
numframes = numel(alltraces1(1,:));
timevector = [1:numframes]./6;
figure
plot(timevector,alltraces1','Color','b','LineWidth',1);
hold on
plot(timevector,nanmean(alltraces1),'Color','k','LineWidth',5)
ylabel 'CEBPB'
xlabel 'Time (hrs)'
