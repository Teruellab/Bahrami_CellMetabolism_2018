% Image extraction
function Fixed_Cell_Cytosol(row,col,site,day)
row = 2;col = 2;site = 2;day = 1;
tic
imagepath = 'Z:\kyle\';

shadingpath='D:\ShadingImages\DCYTC 10x\';

experimentpath='Fixed_Cell_Data\170331_siMEK_siERK_siARAF_Hela_HUVEC\170403-Hela-HUVEC-MCF10A-Many-Correlations\'; 
datadir = [imagepath,experimentpath,'Data\'];

if ~exist(datadir,'dir')
    mkdir(datadir)
end

rawdir = [imagepath, experimentpath,'Raw\'];
biasdir=[rawdir,'Bias\'];
%load the bias correction
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos = BG;
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos = BG;  
% bgcmos = zeros(1080,1080);
% define sites to visit

% channels = {'hoechst','bodipy','pparg','cebpb'};
name1 = 'Hoechst';
name2 = 'YFP';
name3 = 'mChy';
name4 = 'Cy5';

nucr = 4;
blobthreshold = -0.08;
debrisarea = 50;
boulderarea = 500;
% By Plate and Site
load([biasdir,name1,'_',num2str(site),'_',num2str(day),'.mat']);bias1=bias;
load([biasdir,name2,'_',num2str(site),'_',num2str(day),'.mat']);bias2=bias;
load([biasdir,name3,'_',num2str(site),'_',num2str(day),'.mat']);bias3=bias;
load([biasdir,name4,'_',num2str(site),'_',num2str(day),'.mat']);bias4=bias;

% By site only
% load([biasdir,name1,'_',num2str(site),'.mat']);bias1=bias;
% load([biasdir,name2,'_',num2str(site),'.mat']);bias2=bias;
% load([biasdir,name3,'_',num2str(site),'.mat']);bias3=bias;
% load([biasdir,name4,'_',num2str(site),'.mat']);bias4=bias;

shot = [num2str(row),'_',num2str(col),'_',num2str(site)];

raw1 = double(imread([rawdir,shot,'_',name1,'_',num2str(day),'.tif'])); 
raw2 = double(imread([rawdir,shot,'_',name2,'_',num2str(day),'.tif'])); 
raw3 = double(imread([rawdir,shot,'_',name3,'_',num2str(day),'.tif'])); 
raw4 = double(imread([rawdir,shot,'_',name4,'_',num2str(day),'.tif'])); 


raw1=(raw1-bgcmos)./bias1; 
raw2=(raw2-bgcmos)./bias2;
raw3=(raw3-bgcmos)./bias3; 
raw4=(raw4-bgcmos)./bias4;
% raw1= raw1-bgcmos;
% raw1 = imtophat(raw1,strel('disk',50));
% raw2 = raw2-bgcmos;
% raw2 = imtophat(raw2,strel('disk',50));
% raw3 = raw3-bgcmos;
% raw3 = imtophat(raw3,strel('disk',50));

cytochannel = raw4;
cytosegmethod  = 1;
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground = nuc_mask;

% nuc_mask_marker = blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
% nuc_mask_marker = segmentdeflections_bwboundaries(nuc_mask_marker,nucr,debrisarea);
% nuc_mask_marker = imerode(nuc_mask_marker,strel('disk',1));
% nuc_mask = threshmask(raw1);
% foreground = nuc_mask;
% raw1comp = imcomplement(raw1);
% raw1min = imimposemin(raw1comp,~foreground|nuc_mask_marker);
% raw1shed = watershed(raw1min);
% nuc_mask = raw1shed>1;
%        
nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);

[nuc_mask,bad_mask] =excludelargeandwarped(nuc_mask,boulderarea);

anti_border_mask = imclearborder(nuc_mask);
border_mask = nuc_mask-anti_border_mask;
bad_mask = bad_mask+border_mask;

[nuc_label,numcells]=bwlabel(nuc_mask);
bad_nuc_label = bwlabel(bad_mask);
bad_nuc_label = bad_nuc_label+numcells;
bad_nuc_label(bad_nuc_label==numcells) = 0;
nuc_label_all = nuc_label + bad_nuc_label;

% nuc_label_extend = bwmorph(nuc_label,'thicken',3);
% nuc_info_extended = regionprops(nuc_label_extend,'PixelIdxList');

markermask = nuc_label_all>0;

switch cytosegmethod
    case 0
        % do nothing
    case 1
        cytogray = mat2gray(cytochannel);
%         cytoadapt = adapthisteq(cytogray);
        [cytohist,bins] = imhist(cytogray,256);
        threshval = triangle_th(cytohist,numel(bins));
        cytobw = cytogray>threshval;
    case 2
        cytogray = mat2gray(log(cytochannel));
%         cytogray = adapthisteq(cytogray);
        threshval = minerrthresh(cytogray);        
%         threshval = graythresh(cytogray);
        cytobw = cytogray>threshval;
    case 3
        cytogray = mat2gray(cytochannel);
        [counts,edges] = histcounts(cytogray);
        [~,posi] = max(counts);
        threshval = edges(posi+1);
        cytobw = cytogray>threshval;
end
% clean up 
% remove background pixels that come from thresholding 
% (assumes background pixels that are above threshold are much smaller than objects)
cytobw2 = imopen(cytobw,strel('disk',2));
cytobw3 = bwareaopen(cytobw2,debrisarea*2);
cytodilate = imdilate(cytobw3,strel('disk',1));
% cytoadaptsmooth = imfilter(cytologgray,fspecial('average',[3 3]));
% cytoadaptsmooth = imfilter(cytologgray,fspecial('gaussian',4,1));
cytogray = imsharpen(cytogray);
cytocomp = imcomplement(cytogray);
cytomod = imimposemin(cytocomp,~cytodilate|markermask);
cytoshed = watershed(cytomod);
cell_mask = (cytoshed-1)>0;
cell_mask = cell_mask & cytobw3;
cell_mask = imfill(cell_mask,'holes');
cell_mask = bwareaopen(cell_mask,debrisarea);
cell_mask = imclearborder(cell_mask);
cell_label = bwlabel(cell_mask);

good_mask = nuc_mask & cell_mask;
cell_label_temp = cell_label;
cell_label_temp(~good_mask) = 0;
good_cell_inds = unique(cell_label_temp(cell_label_temp>0));

allarray = 1:max(cell_label(:));
badlabelmarkers = setdiff(allarray,good_cell_inds);

cell_label(ismember(cell_label,badlabelmarkers)) = 0;

[cell_label,numcells] = bwlabel(cell_label>0);

cyto_label = cell_label;
cyto_label(markermask==1) = 0;

nuc_label = cell_label;
nuc_label(cyto_label>0) = 0;

cyto_info = regionprops(cyto_label,'Area','PixelIdxList');
cyto_area = [cyto_info.Area]';

%establish nanmask with foreground and bodipy foreground
nanmask = cytodilate;
% nanmaskbodipy = bodipymask;
compression = 2;
% background subtraction to get 'real' intensity values
blurradius = 3;

% blur1=imfilter(raw1,fspecial('disk',blurradius),'symmetric');
blur2=imfilter(raw2,fspecial('disk',blurradius),'symmetric');
blur3=imfilter(raw3,fspecial('disk',blurradius),'symmetric');
blur4=imfilter(raw4,fspecial('disk',blurradius),'symmetric');

real1= bgsubmasked_global_2(raw1,nanmask,1,compression,50);
real2= bgsubmasked_global_2(blur2,nanmask,1,compression,50);
real3= bgsubmasked_global_2(blur3,nanmask,1,compression,50);
real4= bgsubmasked_global_2(blur4,nanmask,1,compression,50);
% real1 = raw1;
% real2 = blur2;
% real3 = blur3;

nuc_info=struct2cell(regionprops(nuc_label,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
nuc_mass=nuc_density.*nuc_area;
nuc_info=regionprops(nuc_label,'PixelIdxList');


nanvec=ones(numcells,1)*NaN;
sig1_nucmean = nanvec; %sig1_nucmedian = nanvec;
sig2_nucmean = nanvec; sig2_nucmedian = nanvec;
sig3_nucmean = nanvec; sig3_nucmedian = nanvec;
sig4_nucmean = nanvec; sig4_nucmedian = nanvec;
% sig1_cytomean = nanvec; sig1_cytomedian = nanvec;
sig2_cytomean = nanvec; sig2_cytomedian = nanvec;
sig3_cytomean = nanvec; sig3_cytomedian = nanvec;
sig4_cytomean = nanvec; sig4_cytomedian = nanvec;
% sig1_cytototal = nanvec;
sig2_cytototal = nanvec;
sig3_cytototal = nanvec;
sig4_cytototal = nanvec;
% sig1_celltotal = nanvec;
sig2_celltotal = nanvec;
sig3_celltotal = nanvec;
sig4_celltotal = nanvec;
for cell=1:numcells
    nucpixels = nuc_info(cell).PixelIdxList;
    cytopixels = cyto_info(cell).PixelIdxList;
    cellpixels = union(nucpixels,cytopixels);
    
    sig1_nucmean(cell) = mean(real1(nucpixels)); %sig1_nucmedian(cell) = median(real1(nucpixels));
    sig2_nucmean(cell) = mean(real2(nucpixels)); sig2_nucmedian(cell) = median(real2(nucpixels));
    sig3_nucmean(cell) = mean(real3(nucpixels)); sig3_nucmedian(cell) = median(real3(nucpixels));
    sig4_nucmean(cell) = mean(real4(nucpixels)); sig4_nucmedian(cell) = median(real4(nucpixels));
    
%     sig1_cytomean(cell) = mean(real1(cytopixels)); sig1_cytomedian(cell) = median(real1(cytopixels));
    sig2_cytomean(cell) = mean(real2(cytopixels)); sig2_cytomedian(cell) = median(real2(cytopixels));
    sig3_cytomean(cell) = mean(real3(cytopixels)); sig3_cytomedian(cell) = median(real3(cytopixels));
    sig4_cytomean(cell) = mean(real4(cytopixels)); sig4_cytomedian(cell) = median(real4(cytopixels));
    
%     sig1_cytototal(cell) = sum(real1(cytopixels));
    sig2_cytototal(cell) = sum(real2(cytopixels));
    sig3_cytototal(cell) = sum(real3(cytopixels));
    sig4_cytototal(cell) = sum(real4(cytopixels));
    
%     sig1_celltotal(cell) = sum(real1(cellpixels));
    sig2_celltotal(cell) = sum(real2(cellpixels));
    sig3_celltotal(cell) = sum(real3(cellpixels));
    sig4_celltotal(cell) = sum(real4(cellpixels));

end


fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,cyto_area,sig1_nucmean,...
    sig2_nucmean,sig2_nucmedian,sig2_cytomean,sig2_cytomedian,sig2_cytototal,sig2_celltotal,...
    sig3_nucmean,sig3_nucmedian,sig3_cytomean,sig3_cytomedian,sig3_cytototal,sig3_celltotal,...
    sig4_nucmean,sig4_nucmedian,sig4_cytomean,sig4_cytomedian,sig4_cytototal,sig4_celltotal,];
save([datadir,'fixdata_','Plate_',num2str(day),'_',shot,'.mat'],'fixdata');



disp([shot,'_Plate_',num2str(day)])

% save([datadir,'fixdata_Day4_',shot,'.mat'],'fixdata');
% disp(shot)
toc
% end

