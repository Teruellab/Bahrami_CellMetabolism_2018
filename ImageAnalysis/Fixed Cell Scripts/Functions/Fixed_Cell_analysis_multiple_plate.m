% Image extraction
function Fixed_Cell_analysis_multiple_plate(row,col,site,plate)

% row = 7;col = 11;site = 1;plate = 1;
tic
imagepath = 'Z:\kyle\CRISPRi\';
shadingpath='C:\Users\Teruel Lab\Desktop\Fixed Cell Scripts\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='160601_CRISPRi_Knockdown\';
datadir = [imagepath,experimentpath,'Data\'];

if ~exist(datadir,'dir')
    mkdir(datadir)
end

rawdir = [imagepath, experimentpath,'Raw\'];
% rawdir = [rawdir,'Day',num2str(day),'\'];
biasdir=[rawdir,'Bias\'];
%load the bias correction
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos = BG;
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos = BG;
% define sites to visit

% names for signals, name1 is nuclear and names 2-4 are signals
name1 = 'Hoechst';
name2 = '488';
name3 = 'mCherry';
name4 = '647';

nucr = 16;
blobthreshold = -0.1;
debrisarea = 150;
boulderarea = 2000;
ringcalc = 0;

% Loading Bias Estimations
load([biasdir,name1,'_',num2str(site),'_',num2str(plate),'.mat']);bias1=bias;
load([biasdir,name2,'_',num2str(site),'_',num2str(plate),'.mat']);bias2=bias;
load([biasdir,name3,'_',num2str(site),'_',num2str(plate),'.mat']);bias3=bias;
load([biasdir,name4,'_',num2str(site),'_',num2str(plate),'.mat']);bias4=bias;

shot = [num2str(row),'_',num2str(col),'_',num2str(site)];

raw1 = double(imread([rawdir,shot,'_',name1,'_',num2str(plate),'.tif']));
raw2 = double(imread([rawdir,shot,'_',name2,'_',num2str(plate),'.tif']));
raw3 = double(imread([rawdir,shot,'_',name3,'_',num2str(plate),'.tif']));
raw4 = double(imread([rawdir,shot,'_',name4,'_',num2str(plate),'.tif']));

% subtract camera noise and correct for illumination bias
raw1=(raw1-bgcmos)./bias1;
raw2=(raw2-bgcmos)./bias2;
raw3=(raw3-bgcmos)./bias3;
raw4=(raw4-bgcmos)./bias4;

% nuclear mask extraction
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground = nuc_mask;


nuc_mask=bwareaopen(nuc_mask,debrisarea);

nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);

[nuc_mask,~] = excludelargeandwarped(nuc_mask,boulderarea);

nuc_mask = imclearborder(nuc_mask);

%Run the code below to assess the mask
%assessmask(raw1,nuc_mask)

%establish nanmask with foreground and bodipy foreground
nanmask = imdilate(foreground,strel('disk',round(nucr)));

compression = 2;
% background subtraction to get 'real' intensity values
blurradius = 3;

blur1=imfilter(raw1,fspecial('disk',blurradius),'symmetric');
blur2=imfilter(raw2,fspecial('disk',blurradius),'symmetric');
blur3=imfilter(raw3,fspecial('disk',blurradius),'symmetric');
blur4=imfilter(raw4,fspecial('disk',blurradius),'symmetric');

real1 = bgsubmasked_global_2(blur1,nanmask,1,compression,50);
real2 = bgsubmasked_global_2(blur2,nanmask,1,compression,50);
real3 = bgsubmasked_global_2(blur3,nanmask,1,compression,50);
real4 = bgsubmasked_global_2(blur4,nanmask,1,compression,40);

[nuc_label,numcells]=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
nuc_mass=nuc_density.*nuc_area;
nuc_info=regionprops(nuc_label,'PixelIdxList');
sig1 = zeros(numcells,1);
sig2 = zeros(numcells,1); %sig2mean = zeros(numcells,1);
sig3 = zeros(numcells,1); %sig3mean = zeros(numcells,1);
sig4 = zeros(numcells,1);
% sig4median = zeros(numcells,1); sig4mean = zeros(numcells,1);
% sig3 = zeros(numcells,1);
nanvec=ones(numcells,1)*NaN;
for cell=1:numcells
    sig1(cell) = mean(real1(nuc_info(cell).PixelIdxList));
    
    sig2(cell) = median(real2(nuc_info(cell).PixelIdxList));
    % sig2mean(cell) = mean(real2(nuc_info(cell).PixelIdxList));
    
    sig3(cell) = median(real3(nuc_info(cell).PixelIdxList));
    % sig3mean(cell) = mean(real3(nuc_info(cell).PixelIdxList));
    sig4(cell)= median(real4(nuc_info(cell).PixelIdxList));
    %sig4median(cell)=median(real4(nuc_info(cell).PixelIdxList));
    % sig4mean(cell) = mean(real4(nuc_info(cell).PixelIdxList));
end

if ringcalc==1
    innerrad=1; outerrad=5; %10xB1|20xB2: 1/5
    ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
    ring_info=regionprops(ring_label,'PixelIdxList');
    ring_vals = regionprops(ring_label,real2,'PixelValues');
    sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_fgmode=nanvec;
    for cell=1:numcells
        if cell>numel(ring_info)
            break;
        end
        ring2all=real2(ring_info(cell).PixelIdxList);
        ring2all(ring2all>prctile(ring2all,98))=[];
        sig2ring_75th(cell)=prctile(ring2all,75);
        ring2foreground=ring2all(ring2all>20);
        histogram(ring2all,25)
        if numel(ring2foreground)<100
            ring2foreground=ring2all;
        end
        if numel(ring2all)>100
            sig2ring_fgmedian(cell)=nanmedian(ring2foreground);
            numbins=25;
            sig2ring_fgmode(cell)=getmode(ring2foreground,numbins);
        end
    end
end


% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,repmat(numcells,numel(sig3),1)];
fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig4];
%fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode,sig3,sig3ring_75th,sig3ring_fgmedian,sig3ring_fgmode];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2sum,sig2median,sig3,sig4];
save([datadir,'fixdata_','Plate_',num2str(plate),'_',shot,'.mat'],'fixdata');
disp([shot,'_Plate_',num2str(plate)])

% save([datadir,'fixdata_Day4_',shot,'.mat'],'fixdata');
% disp(shot)
toc
% end

