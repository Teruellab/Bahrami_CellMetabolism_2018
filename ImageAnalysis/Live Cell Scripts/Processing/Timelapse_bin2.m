function Timelapse_bin2(row,col,site)

% Test specific site by uncommenting this
% row = 2; col = 6; site = 2;

% File Path Setup (Windows)
% imagepath = 'Z:\michael\';
% shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
% experimentpath = 'Live Cell Experiments\20161218-PPARg-FABP4-siRNA-CDKi\';
% shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
% wellname = nameandsite(shot);
% datadir=[imagepath,experimentpath,'Data\'];
% biasdir=[imagepath,experimentpath,'Raw\Bias2\'];
% separatedirectories=0;
% if separatedirectories==1
%     rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
%     maskdir=[imagepath,experimentpath,'Mask2\',shot];
% else
%     rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
%     maskdir=[imagepath,experimentpath,'Mask2\'];
% end
% File Path Setup (Linux)
imagepath = '/mnt/fluffy/michael/';
shadingpath='/mnt/fluffy/michael/Scripts/ShadingImages/DCYTC 10x/';
experimentpath = 'Live Cell Experiments/20170421-CC-Diff-mimi2/';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
wellname = nameandsite(shot,0);
datadir=[imagepath,experimentpath,'Data/'];
biasdir=[imagepath,experimentpath,'Raw/Bias/'];
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw/',shot,'/'];
    maskdir=[imagepath,experimentpath,'Mask/',shot];
else
    rawdir=[imagepath,experimentpath,'Raw/',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask/'];
end

maskwrite = 1;

if ~exist(datadir,'dir')
    mkdir(datadir)
end
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end

if separatedirectories==0
    maskdir=[maskdir,shot,'_'];
end
%%% Tracking setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF = 1;
EF = 508; %SF=starting frame EF=ending frame CHANGE EACH TIME
nucr = 8;
blobr = 12;
deflectionstep = 12;
maxjump=nucr*3; %default nucr*3
debrisarea= 50;
boulderarea= 1500; %MCF-10A: H2B:20 NLS:10
blobthreshold = -0.025;  %default 10xbin1=-0.02 20xbin2=-0.03; NLS=-0.01;2
% blurradius=3;
ringcalc=1; %set to zero if ring is unneeded
name1='Cy5_'; % nuclei channel 
name2='CFP_';
name3='YFP_';
name4='mChy_';
namenucedge='nucedge_';
% adjustframe=[114,386]; %default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
framelength = length(SF:EF);
totalframes=numel(frames);
badframes=ones(framelength,1)*NaN;
if SF>1
    badframes(1:SF-1)=0;
end

jitters=zeros(EF,2);
blocksize=10000;
maxcellnum=blocksize;
if ringcalc==1
    parameternum=10; % 13 for two rings
    tracedata=ones(maxcellnum,EF,parameternum)*NaN;
else
    parameternum = 7; % 6 for 1 signal 7 for 2 signals
    tracedata=ones(maxcellnum,EF,parameternum)*NaN;
end
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
cellcounter = ones(framelength,1);
timetotal=tic;
%[firstgoodindex,blurthreshhigh,numthresh,badframes,height,width]=timelapsesetup(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_3(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
% [firstgoodindex,blurthreshhigh,blurthres;hlow,numthresh,badframes,height,width]=timelapsesetup_4thresh(rawdir,name1,frames,nucr,debrisarea,badframes,maskwrite);
%blurthreshlow=0;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numsets = 1;
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos=BG;
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG;
bias1all = zeros(height,width,numsets);
bias2all = zeros(height,width,numsets);
bias3all = zeros(height,width,numsets);
bias4all = zeros(height,width,numsets);
for seti = 1:numsets
load([biasdir,name1,num2str(site),'_',num2str(1),'.mat']); bias1all(:,:,seti)=bias;
load([biasdir,name2,num2str(site),'_',num2str(1),'.mat']); bias2all(:,:,seti)=bias;
load([biasdir,name3,num2str(site),'_',num2str(1),'.mat']); bias3all(:,:,seti)=bias;
load([biasdir,name4,num2str(site),'_',num2str(1),'.mat']); bias4all(:,:,seti)=bias;
end
% 
% for seti = 1:numsets
%     load([biasdir,name1,num2str(site),'.mat']); bias1all(:,:,seti) = bias;
%     load([biasdir,name2,num2str(site),'.mat']); bias2all(:,:,seti) = bias;
%     load([biasdir,name3,num2str(site),'.mat']); bias3all(:,:,seti) = bias;
%     load([biasdir,name4,num2str(site),'.mat']); bias4all(:,:,seti) = bias;
% end
rhfrac = 1; rwfrac = 1;
regheight=1:rhfrac*height; regwidth=1:rwfrac*width;
for i=firstgoodindex:totalframes
    
    if i <= 690
        bias1 = bias1all(:,:,1);
        bias2 = bias2all(:,:,1);
        bias3 = bias3all(:,:,1);
        bias4 = bias4all(:,:,1);
%     elseif i >= 53 && i <= 120
%         bias1 = bias1all(:,:,2);
%         bias2 = bias2all(:,:,2);
%     elseif i >= 121 && i <= 195
%         bias1 = bias1all(:,:,3);
%         bias2 = bias2all(:,:,3);
%     elseif i >= 196 && i <= 242
%         bias1 = bias1all(:,:,4);
%         bias2 = bias2all(:,:,4);
%     elseif i >= 243 && i <= 253
%         bias1 = bias1all(:,:,5);
%         bias2 = bias2all(:,:,5);
%     elseif i >= 254 && i <= 268
%         bias1 = bias1all(:,:,6);
%         bias2 = bias2all(:,:,6);
%     elseif i >= 269 && i <= 338
%         bias1 = bias1all(:,:,7);
%         bias2 = bias2all(:,:,7);
%     elseif i >= 412 && i <= 482
%         bias1 = bias1all(:,:,8);
%         bias2 = bias2all(:,:,8);
%     elseif i >= 483 && i <= 554
%         bias1 = bias1all(:,:,9);
%         bias2 = bias2all(:,:,9);
%     elseif i >= 555;
%         bias1 = bias1all(:,:,10);
%         bias2 = bias2all(:,:,10);
    end

    f=frames(i); fprintf('frame %0.0f\n',f);
    
    timeframe=tic;
        
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw1 = double(imread([rawdir,name1,num2str(f),'.tif'])); 
    raw1 = (raw1-bgcmos)./bias1;
%     fprintf('1');
    raw2 = double(imread([rawdir,name2,num2str(f),'.tif'])); 
    raw2 = (raw2-bgcmos)./bias2;
%     fprintf('2');
    raw3 = double(imread([rawdir,name3,num2str(f),'.tif'])); 
    raw3 = (raw3-bgcmos)./bias3;
%     fprintf('3');
    raw4 = single(imread([rawdir,name4,num2str(f),'.tif']));
    raw4 = (raw4-bgcmos)./bias4;
%     fprintf('end');486
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if i==firstgoodindex
        %%% Triangle-Threshold filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rawgray = mat2gray(raw1);
%         histvar = imhist(rawgray);
%         thresh = triangle_th(histvar,numel(histvar));
%         
%         [counts,edges] = histcounts(rawgray);
%         [~,maxpos] = max(counts);
%         nucmode = edges(maxpos+1);
%         
%         realthresh = nucmode + (thresh-nucmode)/2;
%         
%         nuc_mask = im2bw(rawgray,realthresh);
% %         nuc_mask = foreground;
% %         nuc_mask=imfill(foreground,'holes');
%         nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %bin1:4 bin2:2
%         nuc_mask=~bwmorph(~nuc_mask,'diag');
%         nuc_mask=~bwmorph(~nuc_mask,'bridge');

%         nuc_mask=markershed(nuc_mask,round(nucr*2/3));
        %%% LoG filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rawgray = mat2gray(raw1);
%         adaptgray = adapthisteq(rawgray);        
        nuc_mask=blobdetector_3_bin2(log(raw1),blobr,blobthreshold,debrisarea);
%         nuc_mask2=blobdetector_3_bin2(log(raw1),nucr,blobthreshold,debrisarea);
%         foreground=logical(nuc_mask);

%         %%% Double-Threshold filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=threshmask(raw1,3);
%         nuc_mask=markershed(nuc_mask,round(nucr*2/3));
%         foreground=nuc_mask;
%         nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea); % org= 2*boulderarea
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask = imfill(nuc_mask,'holes');
%         nuc_mask = bwareaopen(nuc_mask,debrisarea);
        %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask = segmentdeflections_bwboundaries(nuc_mask,deflectionstep,debrisarea);
%         nuc_mask = excludelargeandwarped(nuc_mask,boulderarea);
        %%%% marker controlled water shed %%%%%%
%         markermask = blobdetector_3_bin2(log(raw1),nucr,blobthreshold,debrisarea);
%         markermask = segmentdeflections_bwboundaries(markermask,deflectionstep,debrisarea);
%         maskcents = regionprops(logical(markermask),'Centroid');
%         maskcents = cat(1,maskcents.Centroid);
%         nuc_mask = aprior_markermask(nuc_mask,maskcents,[0 0]);
%         nuc_mask = aprior_markershed_grad(nuc_mask,maskcents,[0 0],imgradient(raw1));
%         nuc_mask=bwareaopen(nuc_mask,debrisarea);
    else
        %%% Triangle-Threshold filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         rawgray = mat2gray(raw1);
%         histvar = imhist(rawgray);
%         thresh = triangle_th(histvar,numel(histvar));
%         
%         [counts,edges] = histcounts(rawgray);
%         [~,maxpos] = max(counts);
%         nucmode = edges(maxpos+1);
%         
%         realthresh = nucmode + (thresh-nucmode)/2;
%         
%         nuc_mask = rawgray>realthresh;
% %         nuc_mask = foreground;
% %         nuc_mask=imfill(foreground,'holes');
%         nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %bin1:4 bin2:2
%         nuc_mask=~bwmorph(~nuc_mask,'diag');
%         nuc_mask=~bwmorph(~nuc_mask,'bridge');
        %%% LoG filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=blobdetector_3_bin2(log(raw1),blobr,blobthreshold,debrisarea);
%         foreground=logical(nuc_mask);
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
        jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        
        %%% LoG filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
%         %%% Double-Threshold filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=threshmask(raw1,3);
%         nuc_mask=markershed(nuc_mask,round(nucr*2/3));
%         foreground=nuc_mask;
%         nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
%         %%% a priori marker-based watershed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         testim = nuc_mask;
%         [nuc_mask,~]=aprior_markermask(nuc_mask,nuc_center,jitters(f,:));
%         markermask = blobdetector_3_bin2(log(raw1),nucr,blobthreshold,debrisarea);
%         markermask = segmentdeflections_bwboundaries(markermask,deflectionstep,debrisarea);
%         maskcents = regionprops(logical(markermask),'Centroid');
%         maskcents = cat(1,maskcents.Centroid);
%         nuc_mask = aprior_markermask(nuc_mask,maskcents,[0 0]);
%         nuc_mask=bwareaopen(nuc_mask,debrisarea);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=bwareaopen(nuc_mask,debrisarea);
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucnanmask = nuc_mask;
    nuc_mask=imclearborder(nuc_mask);
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression1 = 2; % compression 4 for 2160x2160 and 2 for 1080x1080
    siggray = mat2gray(raw2);
    sighist = imhist(siggray);
    sigtrithresh = triangle_oneside(sighist,256);
    sigmaskint = siggray > sigtrithresh;
%     sigmaskfinal = imdilate(sigmask,strel('disk',3));

    [counts,edges] = histcounts(siggray);
    [~,posi] = max(counts);
    sigmode = edges(posi-1);
    sigmask = siggray>sigmode;
    sigmaskstruct= imopen(sigmask,strel('disk',2));
    sigmaskfinal = sigmaskint | sigmaskstruct;
%     nucnanmask = imdilate(foreground,strel('disk',round(1)));
    cellnanmask = nucnanmask | sigmaskfinal;
    cellnanmask = imdilate(cellnanmask,strel('disk',2));
    
%     blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
    blur2=imfilter(raw2,fspecial('disk',3),'symmetric'); % default is blur by 3 for 10x bin1, but 2 for 10x bin2
    blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
    blur4=imfilter(raw4,fspecial('disk',3),'symmetric');
    [real1] = bgsubmasked_global_2(raw1,cellnanmask,1,compression1,50);    
    [real2] = bgsubmasked_global_2(blur2,cellnanmask,1,compression1,50);
%     bgstore2(f) = estbg2;
    [real3] = bgsubmasked_global_2(blur3,cellnanmask,1,compression1,50);
    [real4] = bgsubmasked_global_2(blur4,cellnanmask,1,compression1,50);
    
    
    % save smoothed and bg subtracted images 
%     if (row==4|| row==5)&&(col==2||col==3) 
%     imwrite(uint16(real1),[realdir,name1,num2str(f),'.tif']);
%     imwrite(uint16(real2),[realdir,name2,num2str(f),'.tif']);
%     imwrite(uint16(real3),[realdir,name3,num2str(f),'.tif']);
%     imwrite(uint16(real4),[realdir,name4,num2str(f),'.tif']);
%     end

    % set pixels with negative value to 0
    real1(real1<0) = 0;
    real2(real2<0) = 0;
    real3(real3<0) = 0;
    real4(real4<0) = 0;
    
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)   
        fprintf('badframe: frame %0.0f\n',f);
        badframes(f)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[maskdir,namenucedge,num2str(f),'.tif']);
            %imwrite(uint16(real2),[maskdir,'\',name2,num2str(f),'.tif']);
            %imwrite(uint16(real3),[maskdir,'\',name3,num2str(f),'.tif']);
        end
        continue;
    end
    blurthreshhigh=1.2*mednuc;  %H2B:1.1 NLS:1.08
    blurthreshlow=0.8*mednuc;   %H2B:0.8 NLS:0.95
    numthresh=0.5*numcells;     %H2B:0.5 NLS:0.85
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=nuc_density.*nuc_area;
    curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>firstgoodindex
        nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        H2BvsNLS=1; %1:H2B 2:NLS
        debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
        %[tracedata,curdata,tracking,nuc_label]=adaptivetrack_5(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),H2BvsNLS,debugpackage);
        %[tracedata,curdata,tracking,nuc_label]=adaptivetrack_6(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),H2BvsNLS,debugpackage);
        %[tracedata,curdata,tracking,nuc_label]=adaptivetrack_7_winrad(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),maxjump,H2BvsNLS,debugpackage);
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_8_split(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),maxjump,debrisarea,H2BvsNLS,debugpackage);
%         [tracedata,curdata,tracking,nuc_label]=adaptivetrack_11_LAP(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),maxjump,debrisarea*2,H2BvsNLS,debugpackage);
        badframes(f)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
%     %nuc_inner=imerode(nuc_label,strel('disk',1,0));
%     %extractmask=nuc_label-nuc_inner;
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
%         %imwrite(uint16(real2),[maskdir,name2,num2str(f),'.tif']);
%         %imwrite(uint16(real3),[maskdir,name3,num2str(f),'.tif']);
    end
   
    
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    cellcounter(f) = numlivecells;
    disp(numlivecells)
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec; sig4=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=mean(real1(nuc_info(cc).PixelIdxList));
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
        sig3(n)=mean(real3(nuc_info(cc).PixelIdxList));
        sig4(n)=mean(real4(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        innerrad=1; outerrad=4; %10xB1|20xB2: 1/5
        ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
        ring_info=regionprops(ring_label,'PixelIdxList');
%         bgthreshold = abs(prctile(real2(~cellnanmask),5)); % use only with median bg subtraction, not with mode bg subtraction
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec;% sig2ring_fgmode=nanvec;
%         sig3ring_75th=nanvec; sig3ring_fgmedian=nanvec; sig3ring_fgmode=nanvec;
        for n=1:numlivecells
            cc=cellid(n);
            if cc>numel(ring_info)
                break;
            end
            ring2all=real2(ring_info(cc).PixelIdxList);
             ring2all(ring2all>prctile(ring2all,98))=[];
%             figure,histogram(ring2all)
            sig2ring_75th(n)=prctile(ring2all,75);
            ring2foreground=ring2all(ring2all>0);
            if numel(ring2foreground)<50 %% was 100
                 ring2foreground=ring2all;
            end
%             ring3all=real3(ring_info(cc).PixelIdxList);
%             ring3all(ring3all>prctile(ring3all,98))=[];
%             sig3ring_75th(n)=prctile(ring3all,75);
%             ring3foreground=ring3all(ring3all>0);
%             if numel(ring3foreground)<100
%                  ring3foreground=ring3all;
%             end
            if numel(ring2all)>50
                 sig2ring_fgmedian(n)=nanmedian(ring2foreground);
%                  numbins=12;
%                  sig2ring_fgmode(n)=getmode(ring2foreground,numbins);
%                  sig3ring_fgmedian(n)=nanmedian(ring3foreground);
%                  sig3ring_fgmode(n)=getmode(ring3foreground,numbins);
            end
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc == 1
%         tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode,sig3,sig3ring_75th,sig3ring_fgmedian,sig3ring_fgmode];
%         tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig3];
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig3,sig4];
    else 
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3];  
%         tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2]; 
%         tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2]; 
    end
    if maxcellnum-max(cellid)<blocksize
        tempdata=ones(blocksize,EF,parameternum)*NaN;
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
tracedata4linking = tracedata; jitters4linking = jitters;
% [tracedata,genealogy,jitters] = postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
% save([datadir,'tracedata_',shot,'_NoLink','.mat'],'tracedata','genealogy','jitters','cellcounter');
%[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr,H2BvsNLS);
% [tracedata,genealogy,jitters]=postprocessing(tracedata4linking,cellid,jitters4linking,badframes,tracking,maxcellnum,nucr);
% save([datadir,'tracedata_',shot,'_WithLink','.mat'],'tracedata','genealogy','jitters','cellcounter');
[tracedata,genealogy,jitters] = postprocessing_SigMatch(tracedata4linking,cellid,jitters4linking,badframes,tracking,maxcellnum,nucr);
save([datadir,'tracedata_',shot,'_WithLinkSigBound','.mat'],'tracedata','genealogy','jitters','cellcounter');
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save([datadir,'tracedata_',shot,'_withlink','.mat'],'tracedata','genealogy','jitters','cellcounter');
toc(timetotal);
clear all; 
clear mex;

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}