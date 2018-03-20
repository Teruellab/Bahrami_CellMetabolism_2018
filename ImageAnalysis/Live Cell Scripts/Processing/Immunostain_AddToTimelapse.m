function Immunostain_2_AddToTimelapse(row,col,site)
row=2;col=2;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projectpath='H:\Documents\Projects\';
imagepath='D:\Oscillations Paper Data\';
experimentpath='20160706-siRNA-cebpa-fabp4-mimi2\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';

% shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
wellname = nameandsite(shot);
datadir=[imagepath,experimentpath,'Data\'];
% datadir=([projectpath,experimentpath,'Data_20140606\']);
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'_nucedge_'];
end
IFdir = [imagepath,experimentpath,'Raw\IF\'];
maskwrite=0;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='CFP_';
% nucedgename='nucedge_';
% name2='YFP_';
% name3='mCherry_';
IFname1='YFP_';
IFname2='Cy5_';
ringcalc = 1;
moviebin=2;
if moviebin==1
    nucr = 12;
    debrisarea=250; %MCF-10A: 100
    %nucr=20;
    %debrisarea=600; %MCF-10A 20xBin1
elseif moviebin==2
    nucr = 8;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea = 700;
blobthreshold=-0.03; %MCF10A 10xBin1
% blobthreshold=-0.03; %MCF10A 20xBin1
timetotal=tic;
%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'_WithLinkSigBound.mat'],'tracedata','genealogy','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rawprev=double(imread([rawdir,name1,num2str(totalframes),'.tif']));
rawprev = logical(imread([maskdir,num2str(totalframes),'.tif']));
nuc_mask_prev = imfill(rawprev,'holes');
% [nuc_mask_prev,~]=blobdetector_foreground(log(rawprev),nucr,blobthreshold,debrisarea);
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG;
% load([imagepath,experimentpath,'Raw\Bias\CFP_IF_',num2str(site),'.mat']); sc1=bias;
% load([imagepath,experimentpath,'Raw\Bias\YFP_IF_',num2str(site),'.mat']); sc2=bias;
% load([imagepath,experimentpath,'Raw\Bias\mChy_IF_',num2str(site),'.mat']); sc3=bias;
% load([imagepath,experimentpath,'Raw\Bias\Cy5_IF_',num2str(site),'.mat']); sc4=bias;
% load([imagepath,experimentpath,'Raw\IBG\Cy5_stain_',num2str(site),'.mat']); sc5=inferredbg;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1 = double(imread([IFdir,shot,'_',name1,'IF.tif'])); 
% raw1=(raw1-bgcmos)./sc1;
raw1=(raw1-bgcmos); raw1 = imtophat(raw1,strel('disk',50));
% raw2 = double(imread([IFdir,shot,'_',name2,'IF.tif'])); 
% raw2=(raw2-bgcmos)./sc2;
% raw3=double(imread([IFdir,shot,'_',name3,'stain.tif'])); 
% raw3=(raw3-bgcmos)./sc3;
IFraw1=double(imread([IFdir,shot,'_',IFname1,'IF.tif'])); 
% IFraw1=(IFraw1-bgcmos)./sc3;
IFraw1=(IFraw1-bgcmos); IFraw1 = imtophat(IFraw1,strel('disk',50));
IFraw2=single(imread([IFdir,shot,'_',IFname2,'IF.tif'])); 
% IFraw2=(IFraw2-bgcmos)./sc4;
IFraw2=(IFraw2-bgcmos); IFraw2 = imtophat(IFraw2,strel('disk',50));
% NEfile=[maskdir,nucedgename,'stain.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
%[nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=blobdetector_3_bin2(log(raw1),nucr,blobthreshold,debrisarea);
foreground=nuc_mask;
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
% nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression=2;
% nanmask1=imdilate(foreground,strel('disk',nucr));
% nanmask2=imdilate(foreground,strel('disk',nucr*1.5));
% nanmaskcyto=imdilate(foreground,strel('disk',nucr));
blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
% real1=bgsubmasked_global_2(blur1,nanmask1,1,compression,50);
real1 = blur1;
% blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
% real2=bgsubmasked_global_2(blur2,nanmask2,1,compression,50);
% blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
% real3=bgsubmasked_global_2(blur3,nanmask,1,compression);
IFblur1=imfilter(IFraw1,fspecial('disk',3),'symmetric');
% IFreal1=bgsubmasked_global_2(IFblur1,nanmask1,1,compression,50);
IFreal1 = IFblur1;
IFblur2=imfilter(IFraw2,fspecial('disk',3),'symmetric');
% IFreal2=bgsubmasked_global_2(IFblur2,nanmask2,1,compression,50);
IFreal2 = IFblur2;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label,numcells]=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mass=nuc_density.*nuc_area;
%%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx,reljity]=detectGlobalShift(nuc_mask_prev,nuc_mask);
reljitx=-reljitx; reljity=-reljity;
% [reljitx,reljity]=registerimages(nuc_mask_prev,nuc_mask);
reljitter=[reljitx,reljity];
prevjitter=jitters(totalframes,:);
IFjitter=prevjitter+reljitter;
nuc_center(:,1)=nuc_center(:,1)+IFjitter(1);
nuc_center(:,2)=nuc_center(:,2)+IFjitter(2);
%%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
%%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
debugpackage={rawprev,nuc_mask_prev,nuc_mask,prevjitter,reljitter};
[tracked,nuc_label]=adaptivetrack_IF(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,debugpackage);
numcells=size(tracked,1);
nuc_center=tracked(:,[1 2]);
nuc_area=tracked(:,3);
nuc_mass=tracked(:,4);
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if maskwrite
    extractmask = bwmorph(nuc_label,'remove');
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellid=find(~isnan(nuc_mass));
nuc_info=regionprops(nuc_label,'PixelIdxList');
innerrad=0.5; outerrad=3; %10xB1|20xB2: 1/5
ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,IFreal1);
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN; sig1=nanvec; %sig2=nanvec; sig3=nanvec;
IFsig1ring_50th=nanvec;% IFsig2ring_50th=nanvec;
IFsig1=nanvec; %IFsig2=nanvec;
for cc=cellid'
    
    if cc>numel(ring_info)
        break;
    end
    
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
%     sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
%     sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    IFsig1(cc)=median(IFreal1(nuc_info(cc).PixelIdxList));
%     IFsig2(cc)=median(IFreal2(nuc_info(cc).PixelIdxList));
    if ringcalc==1
    ringall1 = IFreal1(ring_info(cc).PixelIdxList); 
    ringall1(ringall1>prctile(ringall1,95) | ringall1<50) = []; %hist(ringall,25);
%     ringall2 = IFreal2(ring_info(cc).PixelIdxList); 
%     ringall2(ringall2>prctile(ringall2,95) | ringall2<50) = []; %hist(ringall,25);  
    
    IFsig1ring_50th(cc) = prctile(ringall1,50);
%     IFsig2ring_50th(cc) = prctile(ringall2,50); %previously mean

%     truering=ringall(ringall>5);
% 
%     if numel(truering)<50
%          truering=ringall;
%     end
%     if numel(ringall)<2
%          IFsig2ring_mode(cc)=NaN;
%          continue;
%     end             
% 
%     bmin=min(truering); bmax=max(truering); bstep=(bmax-bmin)/25; bins=bmin:bstep:bmax;
%     [kval,xval]=ksdensity(truering,bins);
%     maxidx=find(kval==max(kval),1); %first mode
%     IFsig2ring_mode(cc)=xval(maxidx);
    end
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([imagepath,experimentpath,'Raw\','bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
% IF_bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
% IFsig2=IFsig2-IF_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,IFsig1,IFsig1ring_50th];
save([datadir,'IF2_',shot,'.mat'],'IFdata','IFjitter');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(IFraw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}