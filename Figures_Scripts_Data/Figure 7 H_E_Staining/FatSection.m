
function FatSection(sample,section,site,set)

% sample = 2;section = 3; site = 2;

sampleList = {'1B','1E','2B','2E','3B','3E','4B','4E','5B','5E','6B','6E','7B','7E','8B','8E',...
    '9B','9E','10B','10E','11B','11E','12B','12E','13B','13E','14B','14E','15B','15E','16B','16E',...
    '17B','17E','18B','18E','19B','19E','20B','20E','21B','21E','22B','22E','23B','23E','24B','24E'};

sampleName = sampleList{sample};

shot = [sampleName,'_',num2str(section),'_',num2str(site)];
% shot = '1B_1_2';

experimentdir = 'D:\Experiments\';
imagedir = [experimentdir,'FatHistology\'];
rawdir = [imagedir,'Raw\'];
maskdir = [imagedir,'SegmentedImages\'];
datadir = [imagedir,'Data\'];

if ~exist(maskdir,'dir')
    mkdir(maskdir)
end

if ~exist(datadir,'dir')
    mkdir(datadir,'dir')
end

% Read and smooth image
fatim = imread([rawdir,shot,'.tiff']);
fatim = double(fatim);
fat = imfilter(fatim,fspecial('average',3),'symmetric');

% Mask intense objects
grayfat = mat2gray(fat);
th = minerrthresh(grayfat);
% th = graythresh(grayfat);
fatmask = grayfat<th;
% figure,imshow(fatmask)

% Mask borders with more sensitive thresholding
graylog = mat2gray(log(fat));
loghist = imhist(graylog);
th2 = triangle_oneside(loghist,256);
logmask = graylog<th2;
logmask = imopen(logmask,strel('disk',3));
% 
% logmasktemp = graylog>th2;
% logmasktemp = bwareaopen(logmasktemp,8);
% BWconnect = bwmorph(logmasktemp,'skel','Inf');
% NoiseSize = 64;
% for i = 3:2:floor(sqrt(NoiseSize))
%     BWconnect = imdilate(BWconnect,ones(i));
%     BWconnect = bwmorph(BWconnect,'skel','Inf');
%     % Reduce back result of dilation
%     for j = 1:floor(i/2)
%         BWendpoints = bwmorph(bwmorph(BWconnect,'endpoints'),'shrink',Inf);
%         BWconnect(BWendpoints) = 0;
%     end
% end

% Distance transform mask for borders
distmap = bwdist(logmask);
L = watershed(distmap);
watermap = L>0;

% figure,imshow(watermap);

finalmask = fatmask & watermap;
finalmask = imfill(finalmask,'holes');
% finalmask = bwareaopen(finalmask,100);

labelmask = bwlabel(finalmask);
badlabel = labelmask;
cellprops = regionprops(labelmask,'Eccentricity','Solidity','Area','EquivDiameter','MajorAxisLength','MinorAxisLength');
cellEcc = [cellprops.Eccentricity];
cellSol = [cellprops.Solidity];
cellArea = [cellprops.Area];
cellDiameter = [cellprops.EquivDiameter];
cellMajor = [cellprops.MajorAxisLength];
cellMinor = [cellprops.MinorAxisLength];
gatecriteria = cellEcc >0.9 | cellSol<0.85 | cellArea>30000 | cellArea<200;
% gatecriteria = cellArea<200 & cellArea>150;

badcells = find(gatecriteria);
goodcells = find(~gatecriteria);

badlabel(ismember(labelmask,goodcells)) = 0;
badmask = badlabel > 0;

labelmask(ismember(labelmask,badcells)) = 0;
finalmask = labelmask>0;

anti_border_mask = imclearborder(finalmask);
border_mask = finalmask-anti_border_mask;
badmask = badmask | border_mask;

finalmask = anti_border_mask;

labelIM = goodbadLabel(fatim,finalmask,badmask);
imwrite(labelIM,[maskdir,shot,'.jpg'])

finaldata = cellprops(~gatecriteria);

save([datadir,shot,'_FinalData.mat'],'finaldata')
% end

% finallabel = bwlabel(finalmask);
% figure,assessmask(graylog,finalmask)

% figure,imagesc(fatim,[0 1e4])
% figure,imagesc(finalmask)



% save

% fatshed = markershed(~fatmask,15);

% figure,imshow(fatshed)
