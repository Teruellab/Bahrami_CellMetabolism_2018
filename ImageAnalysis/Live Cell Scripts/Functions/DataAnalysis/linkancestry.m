function linkedtracedata=linkancestry(tracedata,tracestats,samplecells)
linkedtracedata=ones(size(tracedata))*NaN;
for i=1:numel(samplecells)
    orgcellid=samplecells(i);
    cellid=orgcellid;
    condition=true;
    while condition
        goodframes=find(~isnan(tracedata(cellid,:,1)));
        linkedtracedata(orgcellid,goodframes,:)=tracedata(cellid,goodframes,:);
        oldcellid = cellid;
        cellid=tracestats(cellid,4);
        if cellid == oldcellid
            break
        end
        condition=~isnan(cellid);
    end
end