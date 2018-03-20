%Change these parameters
rows= 2:7;
cols= [2:4,6:9];
sites= 1:2;

manualwells = [
                2,4,1;
                2,4,2;
                7,4,3;
                7,6,1;
                7,6,3;
                7,7,3;
                7,8,3];

manualcontrol=0;
numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
shots=numrows*numcols*numsites;
if manualcontrol==1
    shots=size(manualwells,1);
end
time1=tic;
parpool(10)
parfor shot=1:shots
    if manualcontrol==1 
        row=manualwells(shot,1);
        col=manualwells(shot,2);
        site=manualwells(shot,3);   
    else
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites(siteidx);
        colidx=mod(ceil(shot/numsites),numcols);
        if colidx==0
            colidx=numcols;
        end
        col=cols(colidx);
        rowidx=ceil(shot/(numcols*numsites));
        row=rows(rowidx);
    end
    fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'\n']);
    try 
%         FormatFiles(row,col,site);
        %%% Timelapse %%%%%%%%%%%%%%%%%%%%%%%%%
%         Timelapse(row,col,site);
        Timelapse_bin2(row,col,site);
%         Timelapse_bin2_Original_CodeBase(row,col,site);
%         Timelapse_ppargg_marys_exp(row,col,site);
        %%% Immunostain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Immunostain(row,col,site);
%         Immunostain_AddToTimelapse(row,col,site) 
    catch
        disp(['Error: ',num2str(row),'_',num2str(col),'_',num2str(site)]);
    end
end
toc(time1)