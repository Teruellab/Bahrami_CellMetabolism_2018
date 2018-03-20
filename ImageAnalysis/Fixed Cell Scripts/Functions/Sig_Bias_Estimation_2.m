function Sig_Bias_Estimation_2(ExperimentPath,nucr,nucName,sigName,rows,cols,bin,curplate,cursite)

blobthreshold =-0.03;
debrisarea = 150;


shadingpath = 'D:\ShadingImages\DCYTC 10x\';

biasdir=[ExperimentPath,'Raw\','Bias\'];
% make bias directory if one doesn't already exist
if ~exist(biasdir,'dir')
    mkdir(biasdir)
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowmat=rows;
colmat=cols;
numrows=numel(rowmat); 
numcols=numel(colmat); 
numshots = numrows*numcols;
maskdilation = nucr*2;
if bin == 1
    load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos=BG; 
elseif bin == 2
    load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG;  
end

[height,width]=size(bgcmos);
biasStack = cell(numshots,1);

load([biasdir,nucName,'_',num2str(cursite),'_',num2str(curplate),'.mat'],'bias');
nucbias=bias;
clear bias

parfor shotInd = 1:numshots

    colidx=mod(shotInd,numcols);
    if colidx==0
        colidx=numcols;
    end

    col=colmat(colidx);

    rowidx=ceil(shotInd/(numcols));
    row=rowmat(rowidx);

    shot=[num2str(row),'_',num2str(col),'_',num2str(cursite)];
%             wellname = nameandsite(shot);

    % read in image
    rawdir=[ExperimentPath,'Raw\',shot,'_'];
    nucraw = double(imread([rawdir,nucName,'_',num2str(curplate),'.tif']));
    sigraw = double(imread([rawdir,sigName,'_',num2str(curplate),'.tif']));
%     nucraw = double(imread([rawdir,nucName,'.tif']));
%     sigraw = double(imread([rawdir,sigName,'.tif']));
    
    nucraw=(nucraw-bgcmos)./nucbias;
    sigraw = sigraw-bgcmos;

    % find nuclei
%     nucmask = threshmask(nucraw);  
    nucmask = blobdetector_3(log(nucraw),nucr,blobthreshold,debrisarea);
    sigblur = imfilter(sigraw,fspecial('gaussian',round(nucr/2)),'symmetric');    
    backgroundonly = sigblur;
    % Remove foreground pixels by setting them to NaN (Not a Number)
    nanmask = imdilate(nucmask,strel('disk',maskdilation,0));
    backgroundonly(nanmask) = NaN;
    % find percentile scores in each image block 
    blocknum=12;
    prctilethresh=10;
    %background=blockpercentile(backgroundonly,blocknum,prctilethresh);
    bgblock=blockpercentile_blockimage(backgroundonly,blocknum,prctilethresh);

    midrc=ceil(blocknum/2);
    refval=bgblock(midrc);
    if ~isnan(refval)
        %bgblocknorm=bgblock/max(bgblock(:));
        bgblocknorm=bgblock/refval;
        biasStack{shotInd} = bgblocknorm;
    end
end
biasStack = cat(3,biasStack{:});
blockbias=nanmedian(biasStack,3);
bias=imresize(blockbias,[height width],'bicubic');
%%%%save data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([ExperimentPath,'Raw\','Bias\',sigName,'_',num2str(cursite),'_',num2str(curplate),'.mat'],'bias');



