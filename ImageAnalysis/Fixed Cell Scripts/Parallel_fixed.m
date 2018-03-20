rows = [2:7];
cols = [2:11];
sites = 1:4;
plates = 1;

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
numplates = length(plates);
shots = numrows*numcols*numsites*numplates;

parfor shot = 1:shots

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
    
    rowidx=mod(ceil(shot/(numcols*numsites)),numrows);
    if rowidx==0
        rowidx=numrows;
    end
    row=rows(rowidx);
    
    plateidx = ceil(shot/(numcols*numsites*numrows));
    plate = plates(plateidx);
    
    
 try

     Fixed_Cell_analysis_multiple_plate(row,col,site,plate);
%     Fixed_Cell_2Ring(row,col,site,plate);
%     FatSection(row,col,site,plate)
%     FISH_multistain(row,col,site,plate)


    fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'_',num2str(plate),'\n']);  
 catch
     disp(['Error: ',num2str(row),'_',num2str(col),'_',num2str(site),'_Plate_',num2str(plate)]);
 end

end
