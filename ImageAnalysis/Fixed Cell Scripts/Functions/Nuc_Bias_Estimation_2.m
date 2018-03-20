function Nuc_Bias_Estimation_2(ExperimentPath,nucr,nucName,rows,cols,bin,curplate,cursite)

blobthreshold = -0.03;


shadingpath = 'D:\ShadingImages\DCYTC 10x\';

biasdir=[ExperimentPath,'Raw\','Bias\'];
% make bias directory if one doesn't already exist
if ~exist(biasdir,'dir')
    mkdir(biasdir)
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelName = nucName;
rowmat=rows;
colmat=cols;
numrows=numel(rowmat); 
numcols=numel(colmat); 
numshots = numrows*numcols;

if bin == 1
    load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos=BG; 
    debrisarea = 100;
elseif bin == 2
    load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG; 
    debrisarea = 150;
end

[height,width]=size(bgcmos);
biasStack = cell(numshots,1);

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
    nucraw=double(imread([rawdir,channelName,'_',num2str(curplate),'.tif']));
%     nucraw=double(imread([rawdir,channelName,'.tif']));
    nucraw=nucraw-bgcmos;

    % find nuclei
%     nucmask=threshmask(nucraw);
%     nucmask=imdilate(nucmask,strel('disk',round(nucr/2),0));

    nucmask = blobdetector_3(log(nucraw),nucr,blobthreshold,debrisarea);
    nucmask=imdilate(nucmask,strel('disk',round(nucr/2),0));

    % find percentile scores in each image block 
    blocknum=12;
    prctilethresh=10;
    blurnan=imfilter(nucraw,fspecial('disk',3),'symmetric'); %10x:3, 10xbin2:2, 20x:6
    blurnan(nucmask)=NaN;
%                     bg=blockpercentile(blurnan,blocknum,prctilethresh);
    bgblock=blockpercentile_blockimage(blurnan,blocknum,prctilethresh);

        midrc=ceil(blocknum/2);
        refval=bgblock(midrc);
%     refval = median(bgblock);
    if ~isnan(refval)
        %bgblocknorm=bgblock/max(bgblock(:));
        bgblocknorm=bgblock/refval;
%         biasStack=cat(3,biasStack,bgblocknorm);
        biasStack{shotInd} = bgblocknorm;
    end
end
biasStack = cat(3,biasStack{:});
blockbias=nanmedian(biasStack,3);
bias=imresize(blockbias,[height width],'bicubic');
%%%%save data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([ExperimentPath,'Raw\','Bias\',channelName,'_',num2str(cursite),'_',num2str(curplate),'.mat'],'bias');


% delete(poolobj)
