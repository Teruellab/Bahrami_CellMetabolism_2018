%  function ViewTraces
row=4;col=5;site=4;

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
% projectpath='D:\Documents\Projects\';
imagepath='z:\devon\';
%imagepath='E:\';
experimentpath='20170116-CEBPB-PPARG-tnfa-live\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
wellname = nameandsite(shot);
%shot=wellnum2str(row,col,site);
datadir=([imagepath,experimentpath,'Data\']);
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
%     rawdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=8;
framesperhr=5; % 5.5 before
drugspike=0/framesperhr;
drugspike2 = 0/framesperhr;
startframe = 1;
endframe = 481;
frames=startframe:endframe;
channelnames={'CFP_', 'YFP_', 'mChy_'};
edgemask=1; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 0;    %0:one trace per plotn 1:all traces in one plot
selectmode = 0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2 = 1;
plotsignal3 = 0;
IFoption=0; %0:no IFdata 1:IFdata
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1 =1500;
ymin2=0; ymax2=2; %ymax2=1000;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daug1hters 2:no daughters
quiescentanalysis=0;
removequiescent = 0;
%%%%%% Load Data %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datafile =[datadir,'tracedata_',shot,'_nolink',a','genealogy','jitters');'.mat'];
datafile =[datadir,'tracedata_',shot,'_WithLinkSigBound','.mat'];
% load([datadir,'tracedata_',shot,'.mat'],'tracedat
load(datafile,'tracedata','genealogy','jitters');
% [tracedata,tracestats,motherstats,IFdata,IDs]=gathertracedata_3(datadir,datafile,shot,motheroption,daughteroption,IFoption);
[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother]=gathertracedata_mz_1(datadir,datafile,shot,motheroption,daughteroption,IFoption);
% hist(tracedata(:,end,7),50)
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
% reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace= 350;

[tempcellcount,~]=size(tracestats);
alteredlength  = zeros(tempcellcount,1);
alteredstart = zeros(tempcellcount,1);
alteredend = alteredstart;
for i = 1:tempcellcount
    alteredlength(i) = sum(~isnan(tracedata(i,:,1)));
    alteredstart(i) = find(~isnan(tracedata(i,:,1)),1,'first');
    alteredend(i) = find(~isnan(tracedata(i,:,1)),1,'last');
end
tracestats(:,1) = alteredstart;
tracestats(:,2) = alteredend;
tracestats(:,3) = alteredlength;
badlengths=tracestats(:,3)<minlengthtrace;
%%% gate PPARg data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelPPARg=6;nucareachannel = 3;
noisethresh = 10^5;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=10^10;  %threshold above which max of each trace must be %50
minpos = 0;      %0:anywhere 1:firstframe 2:lastframe
minthresh = 0; %threshold below which min of each trace must be %50
% [traces1,badtraces1]=gate_pparg_2_singlesite(tracedata,tracestats,noisethresh,channelPPARg,nucareachannel,minthresh,minpos,maxthresh,maxpos);
[traces1,badtraces1]=gate_pparg_1_allsites(tracedata,tracestats,noisethresh,channelPPARg,nucareachannel,minthresh,minpos,maxthresh,maxpos);
% [traces2,badtraces2]=gate_pparg_8_mother(tracedata,tracestats,noisethresh,channelPPARg,minthresh,minpos,maxthresh,maxpos);
%%% gate Degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelGem=7; 
% maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
% maxthresh=100;  %threshold above which max of each trace must be %50
% minpos=0;      %0:anywhere 1:mothertrace 2:daughtertrace
% minthresh=25; %threshold below which min of each trace must be %50
% [traces1,badtraces1]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
% [traces2,badtraces2]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nucchannel=6; cytochannel=7;
% nucthreshoption=0;
%0:normalize each trace by percentile of entire trace
%1:normalize by first frame of trace (regardless of whether mitosis or not)
% nucthresh=100;    %threshold for nuc intensity according to nucthreshoption

% motherthresh=0;   %threshold for max DHB ratio of mother. default gating=1.0.  0:ignore
% noisethresh=0.4;  %threshold for max positive jump in cyto/nuc ratio
% [traces3,badtraces3]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis,removequiescent);
%%% gate miscellaneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces2=tracedata(:,:,3); % Area is 3, Mass(total h2b intensity) is 4
traces3 = tracedata(:,:,4);
% traces3 = tracedata(:,:,8);

% orphan = isnan(tracestats(:,4));
badstart = tracestats(:,1)> 72;
badend = tracestats(:,2) <432;
overlap = zeros(tempcellcount,1);
repeatedcelltrace = unique(lastcellmother,'sorted');
repeatedcelltrace(repeatedcelltrace==-1)=[];
overlap(repeatedcelltrace) = 1; overlap(lastcellmother==-1)=0;
toohigh = max(traces1,[],2) > ymax1;
%%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces = badstart | badend | overlap ;% | badtraces2;% | orphan;  
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
traces3=traces3(~badtraces,:);
tracedata=tracedata(~badtraces,:,:);
tracestats=tracestats(~badtraces,:);
markedmitosis = markedmitosis(~badtraces);
IDs=IDs(~badtraces,:);

% traces1 = traces1(:,144:end); frames = frames(144:end); 
% choicecells =[21,68,76,80,97,98,103,108,123,127,131,136,140,146,161,162,168,175,179,181,215,219,222,229,233]; clone1
% choicecells = [33,38,40,49,50,51,52,54,60,62,80,95,105,112,131,136,154,171,200,202,206,216,253,274]; %clone2
% choicecells = [6,7,16,20,27,28,29,36,39,40,49,61,66,98,108,109,110,115,116,112,128,140,141,151,156,195]; % clone3
% numcells = numel(choicecells);
% newtraces = zeros(numcells,130);
% for cell = 1:numcells
%     idx = find(choicecells(cell)==IDs);
%     newtraces(cell,:) = traces1(idx,:);
% end
% clear traces1
% traces1 = newtraces;

%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
% traces1=normalizetraces_3(traces1,tracestats,0);
% traces2=normalizetraces_3(traces2,tracestats,0);
% traces3=normalizetraces_3(traces3,tracestats,0);

numgated=size(traces1,1);
if numgated>100
    numgated=100;
end

%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=1; %default=2 presentation=4
steps=5; ystep1=round(((ymax1-ymin1)/steps)*10)/10; ystep2=round(((ymax2-ymin2)/steps)*10)/10;
trace2color=[0 0 1];
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
%     drugtime2=drugspike2;
%     drugplot2=0;
else
    drugtime=0;
    drugplot=drugspike;
%     drugtime2=0;
%     drugplot2=drugspike2;
end
xtime = frames./framesperhr;
fitTime=frames(143:end)/framesperhr;
% xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(IDs,selectin));
end
if ensemble
    figure; hold on;
end
for counter=1:length(selection)
    i=selection(counter);
    if selectmode
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        plotfignum=gcf;
        %set(gca,'Position',[0.03 0.1 0.95 0.8]);
    elseif ~ensemble
        figure(ceil(counter/24));
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(traces1(i,:),5); %default 3
    if ~plotsignal2
%         truncatedSig = ysig1(143:end);
%         yStart = find(~isnan(truncatedSig),1,'first');
%         yEnd = find(~isnan(truncatedSig),1,'last');
%         fitline = fit(fitTime(yStart:yEnd)',truncatedSig(yStart:yEnd)','exp2');
%         plot(fitline,fitTime,truncatedSig)
%         line(xtime,ysig1,'DisplayName',num2str(i),'LineStyle','none','Marker','.','MarkerSize',2);
        line(xtime,ysig1,'DisplayName',num2str(i),'Color',[0 .7 0],'LineWidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
    elseif plotsignal2
        ysig2=smoothignorenans(traces2(i,:),3);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1,'YColor','k','FontName','Arial','FontSize',8);
        set(hline1,'DisplayName',num2str(i),'color',[0.9 .75 0],'linewidth',tracewidth);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike>0
%         line([drugplot drugplot],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
        %line([drugtime drugtime],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    
    if drugspike2>0
%         line([drugplot2 drugplot2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
        line([drugtime2 drugtime2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
    end
    %%%%% mark mitosis (only the current cell's birth) %%%%%%%%%%%%%%%%%%%%
    if ~isnan(tracestats(i,4))
%         plot(xtime(markedmitosis{i}),ysig1(markedmitosis{i}),'ro','markerfacecolor', 'g','markersize',dotsize);
    end
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2,'FontName','Arial','FontSize',8);
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if plotsignal3
        trace3color = [0 0.8 0];
        ysig3=smoothignorenans(traces3(i,:),6);
        line(xtime,ysig3,'color',trace3color,'DisplayName',num2str(i),'linewidth',tracewidth);
    end
    title(strcat('IDs=',num2str(IDs(i)),', row=', num2str(row),', col=',num2str(col),', site=',num2str(site)));
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        selectmodedisplay(edgemask);
    end
end
fprintf([num2str(numgated),'\n']);
%{
title('Sample Trace: TFEB-EYFP');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 4 3]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}