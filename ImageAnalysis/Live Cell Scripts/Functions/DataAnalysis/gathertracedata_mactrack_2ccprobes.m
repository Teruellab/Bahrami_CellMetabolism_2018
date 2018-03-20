function [tracedata,tracestats,motherstats,IFdata,samplecellsID,markedmitosis,lastcellmother]=gathertracedata_mactrack_2ccprobes(datadir,datafile,shot,motheroption,daughteroption,IFoption)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(datafile); %CellMeasurements.MedianIntensity_nuc2,CellMeasurements.IntegratedIntensity_nuc2,...
tracedata = cat(3, CellMeasurements.MedianNuc1,CellMeasurements.IntegratedNuc1,...
    CellMeasurements.MedianNuc3,CellMeasurements.IntegratedCyto3,...
    CellMeasurements.IntegratedNuc4);
genealogy = CellMeasurements.CellData(:,5);
genealogy(genealogy==0) = nan;
% load(
if IFoption
    load([datadir,'IF_',shot,'.mat'],'IFdata'); % change back to IF instead of IF2
else
    IFdata=ones(size(genealogy))*NaN;
end
%%% gate by genealogy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gatemother=zeros(size(genealogy));
mothers=genealogy(~isnan(genealogy));
cellid=(1:length(genealogy))';
nonmothers=cellid(~ismember(cellid,mothers));
if motheroption==0 %no gating
    gatemother(:)=1;
elseif motheroption==1 %only traces that end in mitosis
    gatemother(mothers)=1;
elseif motheroption==2 %only traces that don't end in mitosis
    gatemother(nonmothers)=1;
end
if daughteroption==0
    gatedaughter=ones(size(genealogy)); %no gating
elseif daughteroption==1
    gatedaughter=~isnan(genealogy); %only traces that start with mitosis
elseif daughteroption==2
    gatedaughter=isnan(genealogy); %only traces that don't start with mitosis
end
gategenealogy=gatemother & gatedaughter;
%%% gate by presence of IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    gateIF=~isnan(IFdata(:,1));
else
    gateIF=ones(size(genealogy));
end
%%% combine gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=gategenealogy & gateIF;
samplecellsID=find(samplecells);
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
tracedata=linkancestry(tracedata,tracestats,samplecellsID);
[markedmitosis,lastcellmother] = trackmitosis(tracestats,samplecellsID);
%%% record mother stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if daughteroption==1
    motherstats=getmotherstatsonly(tracedata,tracestats,samplecellsID);
else
    motherstats=ones(size(genealogy))*NaN;
end
%%% remove gated out cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=tracedata(samplecells,:,:);
tracestats=tracestats(samplecells,:);
motherstats=motherstats(samplecells,:);
IFdata=IFdata(samplecells,:);
% markedmitosis = markedmitosis(samplecells);
% lastcellmother = lastcellmother(samplecells);
end