function Immunostain_2(row,col,site)
plate='9hr';
row=3;col=5;site=9;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='Z:\zahra\';
%imagepath='D:\Images\';
imagepath='Z:\zahra\';
%cmosoffsetfile='D:\Images\CMOS_Offset\cmosoffset_bin1.mat';
cmosoffsetfile='D:\ShadingImages\DCYTC 10x\CameraNoise_1080_bin2.mat';
experimentpath='ZB_20160321_mRNA_decay_FISH\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'Data\',plate,'_'];
%biasdir=[imagepath,experimentpath,'Bias\'];
separatedirectories=0;
if separatedirectories
    rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir = [imagepath,experimentpath,'Mask\',shot];
else
    rawdir = [imagepath,experimentpath,'Images\','Raw ',plate,'\',shot,'_'];
    maskdir = [imagepath,experimentpath,'Mask'];
end
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nameA1='Hoechst'; %nuc
nameA2='FISH';

nameB1='Hoechst';
nameB2='IF';

nucedgename='nucedge';
ringcalc=1;
%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
thickenradius=2*nucr;
debrisarea=200;
boulderarea=1500;
timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(cmosoffsetfile,'BG'); %load CMOS offset
cmosoffset = BG;
[height,width]=size(cmosoffset);
%load([biasdir,name1,'_Stain_',num2str(site),'.mat']); bias1=bias;
%load([biasdir,name2,'_Stain_',num2str(site),'.mat']); bias2=bias;
%load([biasdir,name3,'_Stain_',num2str(site),'.mat']); bias3=bias;
%load([biasdir,name4,'_Stain_',num2str(site),'.mat']); bias4=bias;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawA1=double(imread([rawdir,nameA1,'_stain_1.tif'])); %raw1=(raw1-bgcmos)./bias1;
rawA2=double(imread([rawdir,nameA2,'_stain_1.tif'])); %raw2=(raw2-bgcmos)./bias2;

rawB1=double(imread([rawdir,nameB1,'_stain_2.tif'])); %raw1=(raw1-bgcmos)./bias1;
rawB2=double(imread([rawdir,nameB2,'_stain_2.tif'])); %raw2=(raw2-bgcmos)./bias2;

if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
% %%% remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% foregroundthresh=1000;
% areathresh=20000;
% bgprctile=20;
% [rawA1,stainflag]=removesmears_2(rawA1,foregroundthresh,areathresh,bgprctile);
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segmethod='single';
switch segmethod
    case 'concavity'
        blurradius=5;
        blur=imfilter(rawA1,fspecial('gaussian',blurradius),'symmetric');
        nuc_maskA=ThreshImage(blur);
        nuc_maskA=markershed(nuc_maskA,round(nucr*2/3));
    case 'log'
        blobthreshold=-0.03;
        nuc_maskA=blobdetector_3(log(rawA1),nucr,blobthreshold,debrisarea);
    case 'single'
        blurradius=3;
        nuc_maskA=threshmask(rawA1,blurradius);
        nuc_maskA=markershed(nuc_maskA,round(nucr*2/3));
    case 'double'
        blurradius=3;
        nuc_maskA=threshmask(rawA1,blurradius);
        nuc_maskA=markershed(nuc_maskA,round(nucr*2/3));
        nuc_maskA=secondthresh(rawA1,blurradius,nuc_maskA,boulderarea*2);
end
nuc_maskA=imfill(nuc_maskA,'holes');
nuc_maskA=bwareaopen(nuc_maskA,debrisarea);
nuc_maskA=segmentdeflections_bwboundaries(nuc_maskA,nucr,debrisarea);
nuc_maskA=logical(nuc_maskA);
%%% segment second nuclear image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_maskB=threshmask(rawB1,3);
nuc_maskB=markershed(nuc_maskB,nucr*2/3);
%%% align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitxAB,reljityAB]=registerimages(nuc_maskA,nuc_maskB);
jitmatx=[reljitxAB];
jitmaty=[reljityAB];
cropcoors=getcropcoors([height width],jitmatx,jitmaty);
nuc_maskA=nuc_maskA(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
rawA1=rawA1(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
rawA2=rawA2(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
rawB2=rawB2(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%%% Note border/large/warped objects for later removal %%%%%%%%%%%%%%%%%%%%
antiborder_mask=imclearborder(nuc_maskA);
border_mask=nuc_maskA-antiborder_mask;
antilargewarped_mask=excludelargeandwarped(nuc_maskA,boulderarea);
largewarped_mask=nuc_maskA-antilargewarped_mask;
badmask=border_mask | largewarped_mask;
goodmask=logical(nuc_maskA-badmask);
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanmaskcyto=imdilate(nuc_maskA,strel('disk',nucr*2));
blurA1=imfilter(rawA1,fspecial('disk',3),'symmetric');
blurB2=imfilter(rawB2,fspecial('disk',3),'symmetric');
realA1=bgsubmasked_global(blurA1,nanmaskcyto,1,1);
realA2=imtophat(rawA2,strel('disk',2,0));
realB2=bgsubmasked_global(blurB2,nanmaskcyto,1,1);
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_568to647.mat'],'bleedthroughrate');
% real4bleedthrough=real3*bleedthroughrate(2)+bleedthroughrate(1);
% real4=real4-real4bleedthrough;
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% label both good and bad nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label_good,numgood]=bwlabel(goodmask);
nuc_label_bad=bwlabel(badmask);
nuc_label_bad=nuc_label_bad+numgood;
nuc_label_bad(nuc_label_bad==numgood)=0;
nuc_label_all=nuc_label_good+nuc_label_bad;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=regionprops(nuc_label_all,'Centroid','Area','PixelIdxList');
outer_label=labelthicken(nuc_label_all,thickenradius);
puncta_info_A2_low=getpuncta(outer_label,realA2,350); %org:500
puncta_info_A2_med=getpuncta(outer_label,realA2,400); %org:1000
puncta_info_A2_high=getpuncta(outer_label,realA2,1000); %org:2000
nanvec=ones(numgood,1)*NaN;
xcoor=nanvec; ycoor=nanvec; nuc_area=nanvec;
sigA1=nanvec;
sigA2low=nanvec; sigA2med=nanvec; sigA2high=nanvec;
sigB2=nanvec;
for cc=1:numgood
    xcoor(cc)=nuc_info(cc).Centroid(1);
    ycoor(cc)=nuc_info(cc).Centroid(2);
    nuc_area(cc)=nuc_info(cc).Area;
    sigA1(cc)=mean(realA1(nuc_info(cc).PixelIdxList));
    if (cc<=numel(puncta_info_A2_low)), sigA2low(cc)=puncta_info_A2_low(cc).Area; end
    if (cc<=numel(puncta_info_A2_med)), sigA2med(cc)=puncta_info_A2_med(cc).Area; end
    if (cc<=numel(puncta_info_A2_high)), sigA2high(cc)=puncta_info_A2_high(cc).Area; end
    sigB2(cc)=median(realB2(nuc_info(cc).PixelIdxList));
end
if ringcalc==1
    innerrad=1; outerrad=5; %10xB1|20xB2: 1/5
    ring_label=getcytoring_thicken(nuc_label_all,innerrad,outerrad,realB2);
    ring_info=regionprops(ring_label,'PixelIdxList'); 
    sigB2ring=nanvec;
    ringthresh=50;
    for cell=1:numgood
        if cell>numel(ring_info)
            break;
        end

        ringB2all=real2(ring_info(cell).PixelIdxList);
        ringB2all(ringB2all>prctile(ringB2all,99))=[];
        ring2foreground = ringB2all(ringB2all>0);

        if numel(ringB2all)>ringthresh
             sigB2ring(cell)=nanmedian(ring2foreground);
        end
         
    end
end
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[xcoor,ycoor,nuc_area,sigA1,sigA2low,sigA2med,sigA2high,sigB2,sigB2ring];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
end

function puncta_info=getpuncta(outer_label,raw,thresh)
puncta_mask=raw>thresh;
puncta_label=outer_label.*puncta_mask;
puncta_info=regionprops(puncta_label,'PixelIdxList','Area');
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(goodmask,'remove');
tempframe=imadjust(mat2gray(rawA1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(smear_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

anti_mask=bwareaopen(nuc_mask,debrisarea);
temp_mask=nuc_mask-anti_mask;
extractmask=bwmorph(temp_mask,'remove');

anti_mask=bwareaopen(smear_mask,10000);
extractmask=bwmorph(anti_mask,'remove');
%}