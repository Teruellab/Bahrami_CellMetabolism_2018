%%% Windows file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imagepath='Z:\michael\';
% shadingpath='D:\Michael\ShadingImages\DCYTC 10x\'; % Never to be changed
% %shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
% experimentpath='Live Cell Experiments\20161218-PPARg-FABP4-siRNA-CDKi\';
% % staindir = [experimentpath,'Stains\'];
% biasdir=[imagepath,experimentpath,'Raw\Bias2\'];
% %biasdir=[imagepath,experimentpath,'Raw\','Bias\'];

%%% linux file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath='/home/michael/Experiments/';
shadingpath='/mnt/fluffy/michael/Scripts/ShadingImages/DCYTC 10x/'; % Never to be changed
experimentpath='20170524-mimi2-cebpb/';
biasdir=[imagepath,experimentpath,'Raw/Bias/'];

if ~exist(biasdir,'dir')
    mkdir(biasdir)
end

separatedirectories=0;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
% names={
%     % 'hoechst';
%     'CFP';
% };
name = 'Cy5';
rowmat = [2:7];
colmat = [2:10];
sitemat = 1:3;
framemat = [
   1 25 50 75 100 125 150 175 200 225 250 275 300 325];
% framesperset = [
%    1 25 50 75 100 125 150 175 200 250 300 325 350 375 400 450 500 550 600];
% [frameSets,~] = size(framesperset);
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat);% numframes=numel(framemat);
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos=BG; 
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG; 
% bgcmos = imresize(bgcmos,0.5); % for bin 2 images 1080x1080
[height,width]=size(bgcmos);
debrisarea = 50;
blobthreshold = -.02;
% for set = 1:frameSets
%     framemat = framesperset(set,:);
%     framemat = [1 10 20 30 40 50 60];
% for i=1:length(names)
set=1;
    for s=1:numsites  %for or parfor
        site = sitemat(s);
        nucName = name;
%         nucWellBiasCalc(imagepath,experimentpath,biasdir,nucName,site,rowmat,colmat,framemat,nucr,blobthreshold,debrisarea,bgcmos,height,width,separatedirectories,set);
        
        nucWellBiasParallel(imagepath,experimentpath,biasdir,nucName,site,rowmat,colmat,framemat,nucr,blobthreshold,debrisarea,bgcmos,height,width,separatedirectories,set);
        
%         biasstack=[];
%         for c=1:numcols
%             col=colmat(c);
%             for r=1:numrows
%                 row=rowmat(r);
%                 for f=1:numframes
%                     frame=framemat(f);
%                     shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%                     wellname = nameandsite(shot);
%                     if separatedirectories
%                         rawdir=[imagepath,experimentpath,'Raw\Day4\'];
%                         
%                         %raw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
%                         %raw=single(imread([rawdir,shot,'\',names{i},'_',num2str(frame),'.tif']));
%                         nucraw=double(imread([rawdir,shot,'_',names{i},'.tif']));
%                     else
%                         % use this for regular time lapse stuff
%                         rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
%                         nucraw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
%                         % use this for fixed cell measurements of bias
% %                         rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
% %                         rawdir=[imagepath,staindir,shot,'_'];
% %                         nucraw=double(imread([rawdir,names{i},'.tif']));
%                         
%                     end
%                     nucraw=nucraw-bgcmos;
%                     nucmask=blobdetector_3_bin2(log(nucraw),nucr,blobthreshold,debrisarea);
% %                     nucmask=threshmask(nucraw,3);
%                     nucmask=imdilate(nucmask,strel('disk',round(nucr/2),0));
% 
%                     blocknum=31;
%                     prctilethresh=10;
%                     blurnan=imfilter(nucraw,fspecial('disk',3),'symmetric'); %10x:3, 10xbin2:2, 20x:6
%                     blurnan(nucmask)=NaN;
%                     %bg=blockpercentile(blurnan,blocknum,prctilethresh);
%                     bgblock=blockpercentile_blockimage(blurnan,blocknum,prctilethresh);
% 
%                     midrc=ceil(blocknum/2);
%                     refval=bgblock(midrc);
%                     if ~isnan(refval)
%                         %bgblocknorm=bgblock/max(bgblock(:));
%                         bgblocknorm=bgblock/refval;
%                         biasstack=cat(3,biasstack,bgblocknorm);
%                     end
%                 end
%             end
%         end
%         blockbias=nanmedian(biasstack,3);
%         bias=imresize(blockbias,[height width],'bicubic');
%         %%%%for immunostain
% %         save([biasdir,names{i},'_Stain_',num2str(site),'.mat'],'bias');
% %         save([biasdir,names{i},'_',num2str(site),num2str(1),'.mat'],'bias');
%         %%%%for regular timelapse bias
%         save([biasdir,names{i},'_',num2str(site),'.mat'],'bias');
    end
% end
% end
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nucmask,'remove');
[height,width]=size(nucmask);
RGB=zeros(height,width,3);
RGB(:,:,1)=imadjust(mat2gray(nucraw));
RGB(:,:,2)=extractmask;
figure,imshow(RGB);
%}