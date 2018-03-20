function vIdx=getdeflections(orderedset,nucr)

%%% set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsetshort=nucr/4;
if offsetshort==1
    offsetshort=2;
end
offsetlong=2*offsetshort;
gradientoffsetshort=offsetshort;
gradientoffsetlong=offsetlong;
% shortgradthresh = pi/6; % for bin 1 OP9
% longgradthresh = pi/6; % for bin 1 OP9
shortgradthresh = pi/12; % for bin 2 OP9
longgradthresh = pi/12; % for bin 2 OP9

%%% calculate angular deflection at each point of the boundary %%%%%%%%%%%%
perilength=size(orderedset,1);
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
% find the vector that points from (x) to point (x+offsetshort)
% orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
orderedsetoffsetshort = circshift(orderedset,[-offsetshort, 0]);
shortdiff=orderedsetoffsetshort-orderedset;

% find the angle of the vector based on the Four-Quadrant Inverse from [-pi,pi];
shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians

% find angle difference between x and x+offset
% shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
shortgradoffset= circshift(shortgrad,[-gradientoffsetshort,0]);
shortgraddiff=shortgradoffset-shortgrad;

% make everyting from [0,2pi] instead of [-pi,pi];
shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);  %account for 4 quadrants

%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
% find the vector that points from x to x+offsetlong
% orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
orderedsetoffsetlong = circshift(orderedset, [-offsetlong, 0]);
longdiff = orderedsetoffsetlong-orderedset;

%calculate the angle of the vector [-pi,pi];
longgrad=atan2(longdiff(:,2),longdiff(:,1));

% find angle difference between x and x+offsetlong
% longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
longgradoffset =  circshift(longgrad, [-gradientoffsetlong, 0]);
longgraddiff=longgradoffset-longgrad;

% make everything in the interval [0,2pi]
longgraddiff=longgraddiff+2*pi*(longgraddiff<0);

%%% find deflections above threshold %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% short steps %%%%%%%%%%%%%%%%%%
% shortgraddiff(shortgraddiff>=pi)=0;  %exclude
vIdxmaskshort = shortgraddiff>shortgradthresh & shortgraddiff<pi;

%%%%%%% long steps %%%%%%%%%%%%%%%%%%%
vIdxmasklong = longgraddiff>longgradthresh & longgraddiff<pi;

vIdxmasklong = [zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
vIdxmasklong = imdilate(vIdxmasklong,strel('square',1+nucr/2));
% vIdxmasklong=imdilate(vIdxmasklong,strel('line',1+nucr/2,90));


%%% find local maxima of short steps %%%%%%%%%%%%%%%%%%%%%%%
vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));  %join proximal deflection islands
vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
maxmask=zeros(size(vIdxmaskshort));
for rpc=1:length(vIdxobs)
    pix=vIdxobs(rpc).PixelIdxList;
    [~,index]=max(shortgraddiff(pix));
    maxmask(pix(index)+offsetshort)=1;
end
maxmask=maxmask(1:perilength);  %remove any overhang

%%% find coincidence of long mask & local maxima of short mask %%%%%%%%%%%%
vIdxmask=vIdxmasklong & maxmask;
vIdx=find(vIdxmask);

end