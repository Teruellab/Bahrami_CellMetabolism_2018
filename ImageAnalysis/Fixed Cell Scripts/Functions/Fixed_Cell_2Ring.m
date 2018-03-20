% Image extraction
function Fixed_Cell_2Ring(row,col,site,day)

row = 4;col = 10;site = 1;day = 2;
tic
imagepath = 'D:\Michael\';
%imagepath='E:\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
% experimentpath='11202014-Michael-CellCycle-48hour-glass\';
experimentpath='20150116-Michael-CDK4-titration\'; 
datadir = [imagepath,experimentpath,'Data\'];

if ~exist(datadir,'dir')
    mkdir(datadir)
end

rawdir = [imagepath, experimentpath,'Raw\Day4'];
% rawdir = [rawdir,'Day',num2str(day),'\'];
biasdir=[rawdir,'Bias\'];
%load the bias correction
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos = BG;
% load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos = BG;  
load([shadingpath,'cameranoise_4sites.mat'],'BG'); bgcmos = BG;  
% define sites to visit

% channels = {'hoechst','bodipy','pparg','cebpb'};
name1 = 'Hoechst';
name2 = 'bodipy';
name3 = 'pparg';
% name4 = 'Cy5';

nucr = 12;
blobthreshold = -0.06;
debrisarea = 75;
boulderarea = 500;
ringcalc = 1;
% By Plate and Site
% load([biasdir,name1,'_',num2str(site),'_',num2str(day),'.mat']);bias1=bias;
% load([biasdir,name2,'_',num2str(site),'_',num2str(day),'.mat']);bias2=bias;
% load([biasdir,name3,'_',num2str(site),'_',num2str(day),'.mat']);bias3=bias;
% load([biasdir,name4,'_',num2str(site),'_',num2str(day),'.mat']);bias4=bias;

% By site only
load([biasdir,name1,'_',num2str(site),'.mat']);bias1=bias;
load([biasdir,name2,'_',num2str(site),'.mat']);bias2=bias;
load([biasdir,name3,'_',num2str(site),'.mat']);bias3=bias;
% load([biasdir,name4,'_',num2str(site),'_',num2str(day),'.mat']);bias4=bias;

shot = [num2str(row),'_',num2str(col),'_',num2str(site)];

raw1 = double(imread([rawdir,shot,'_',name1,'_',num2str(day),'.tif'])); 
raw2 = double(imread([rawdir,shot,'_',name2,'_',num2str(day),'.tif'])); 
raw3 = double(imread([rawdir,shot,'_',name3,'_',num2str(day),'.tif'])); 
% raw4 = double(imread([rawdir,shot,'_',name4,'_',num2str(day),'.tif'])); 


% raw1=(raw1-bgcmos)./bias1; 
raw1=(raw1)./bias1; 
raw2=(raw2-bgcmos)./bias2;
raw3=(raw3-bgcmos)./bias3; 
% raw4=(raw4-bgcmos)./bias4;

% raw2 = raw2-bgcmos;
% raw2 = imtophat(raw2,strel('disk',100));
% 
% raw3 = raw3-bgcmos;
% raw3 = imtophat(raw3,strel('disk',50));

nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground = nuc_mask;

% nuc_mask = threshmask(raw1,3);
% foreground = nuc_mask;
% nuc_mask = markershed(nuc_mask,round(nucr*(2/3)));


       
nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);

[nuc_mask,bad_mask] =excludelargeandwarped(nuc_mask,boulderarea);

anti_border_mask = imclearborder(nuc_mask);
border_mask = nuc_mask-anti_border_mask;
bad_mask = bad_mask+border_mask;



%establish nanmask with foreground and bodipy foreground
nanmask = imdilate(foreground,strel('disk',round(nucr*2.5)));
% nanmaskbodipy = bodipymask;
compression = 2;
% background subtraction to get 'real' intensity values
blurradius = 3;

% blur1=imfilter(raw1,fspecial('disk',blurradius),'symmetric');
blur2=imfilter(raw2,fspecial('disk',blurradius),'symmetric');
blur3=imfilter(raw3,fspecial('disk',blurradius),'symmetric');
% blur4=imfilter(raw4,fspecial('disk',blurradius),'symmetric');

real1= bgsubmasked_global_2(raw1,nanmask,1,compression,25);
real2= bgsubmasked_global_2(blur2,nanmask,1,compression,25);
real3= bgsubmasked_global_2(blur3,nanmask,1,compression,25);
% real4=bgsubmasked_global_2(blur4,nanmask,1,compression,50);

% real2 = blur2;
% real3 = blur3;

[nuc_label,numcells]=bwlabel(nuc_mask);
bad_nuc_label = bwlabel(bad_mask);
bad_nuc_label = bad_nuc_label+numcells;
bad_nuc_label(bad_nuc_label==numcells) = 0;
nuc_label_all = nuc_label + bad_nuc_label;

nuc_info=struct2cell(regionprops(nuc_label,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
nuc_mass=nuc_density.*nuc_area;
nuc_info=regionprops(nuc_label,'PixelIdxList');
% nuc_label_extend = bwmorph(nuc_label,'thicken',3);
% nuc_info_extended = regionprops(nuc_label_extend,'PixelIdxList');
sig1 = zeros(numcells,1);
sig2 = zeros(numcells,1); %sig2tot = zeros(numcells,1);
sig3 = zeros(numcells,1); %sig3tot = zeros(numcells,1);
% sig4 = zeros(numcells,1);
% sig4median = zeros(numcells,1); sig4mean = zeros(numcells,1);
% sig3 = zeros(numcells,1);
nanvec=ones(numcells,1)*NaN;
for cell=1:numcells
sig1(cell)=mean(real1(nuc_info(cell).PixelIdxList));

sig2(cell)=median(real2(nuc_info(cell).PixelIdxList));
% sig2tot(cell) = mean(real2(nuc_info_extended(cell).PixelIdxList))*numel(nuc_info(cell).PixelIdxList);
% sig2mean(cell) = mean(real2(nuc_info(cell).PixelIdxList));

sig3(cell)=median(real3(nuc_info(cell).PixelIdxList));
% sig3tot(cell) = mean(real3(nuc_info_extended(cell).PixelIdxList))*numel(nuc_info(cell).PixelIdxList);
% sig3mean(cell) = mean(real3(nuc_info(cell).PixelIdxList));
% sig4(cell)= median(real4(nuc_info(cell).PixelIdxList));
% sig4median(cell)=median(real4(nuc_info(cell).PixelIdxList));
% sig4mean(cell) = mean(real4(nuc_info(cell).PixelIdxList));
end
if ringcalc==1
        innerrad=0.5; outerrad=3; %10xB1|20xB2: 1/5
        ring_label=getcytoring_thicken(nuc_label_all,innerrad,outerrad,real2);
        ring_label_all = imfill(ring_label,'holes');
        ring_info=regionprops(ring_label,'PixelIdxList'); 
        ring_info_all = regionprops(ring_label_all,'PixelIdxList');
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_tot=nanvec;
%         sig3ring_75th=nanvec; sig3ring_fgmedian=nanvec; sig3ring_tot=nanvec;        
        for cell=1:numcells
            if cell>numel(ring_info)
                break;
            end
            ring2all=real2(ring_info(cell).PixelIdxList);
            ring2cell = real2(ring_info_all(cell).PixelIdxList);
            ring2cell(ring2cell<0) = [];
            sig2ring_tot(cell) = sum(ring2cell);
            ring2all(ring2all>prctile(ring2all,99))=[];
            sig2ring_75th(cell)=prctile(ring2all,75);
            ring2foreground = ring2all(ring2all>70);
            
%             ring3all=real3(ring_info(cell).PixelIdxList);
%             ring3cell = real3(ring_info_all(cell).PixelIdxList);
%             ring3cell(ring3cell<0) = [];
%             sig3ring_tot(cell) = sum(ring3cell);            
%             ring3all(ring3all>prctile(ring3all,99))=[];
%             sig3ring_75th(cell)=prctile(ring3all,75);
%             ring3foreground = ring3all(ring3all>5);
            
%             histogram(ring2all,25)
%             if numel(ring2foreground)<100
%                  ring2foreground=ring2all;
%             end
            if numel(ring2all)>50
                 sig2ring_fgmedian(cell)=nanmedian(ring2foreground);
%                  numbins=15;
%                  sig2ring_fgmode(cell)=getmode(ring2foreground,numbins);
            end
%             if numel(ring3all)>50
%                  sig3ring_fgmedian(cell)=nanmedian(ring3foreground);
%                  numbins=15;
%                  sig2ring_fgmode(cell)=getmode(ring2foreground,numbins);
%             end            
        end
end


% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,repmat(numcells,numel(sig3),1)];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_tot,sig2ring_75th,sig2ring_fgmedian,sig3,sig3ring_tot,sig3ring_75th,sig3ring_fgmedian];
fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig3];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2sum,sig2median,sig3,sig4];
save([datadir,'fixdata_','Plate_',num2str(day),'_',shot,'.mat'],'fixdata');
disp([shot,'_Plate_',num2str(day)])

% save([datadir,'fixdata_Day4_',shot,'.mat'],'fixdata');
% disp(shot)
toc
% end

