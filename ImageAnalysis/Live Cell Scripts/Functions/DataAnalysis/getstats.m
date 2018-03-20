function tracestats = getstats(tracedata,genealogy)
numcells=size(tracedata,1);
tracestats=ones(numcells,4)*NaN;
for c=1:numcells
    firsttemp = find(~isnan(tracedata(c,:,1)),1,'first');
    lasttemp = find(~isnan(tracedata(c,:,1)),1,'last');
    if isempty(firsttemp) || isempty(lasttemp)
        firsttemp = 1;
        lasttemp = 2;
    end
    
    tracestats(c,1)=firsttemp;
    tracestats(c,2)=lasttemp;
end
tracestats(:,3)=tracestats(:,2)-tracestats(:,1)+1;
tracestats(:,4)=genealogy;
