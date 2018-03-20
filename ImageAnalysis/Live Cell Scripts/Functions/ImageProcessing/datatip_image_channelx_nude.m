function output_txt = datatip_image_channelx(~,event_obj)

global rawdir maskdir gatedtraces jitters plotfignum immunoframe nucr channel
winrad=3*nucr;    %window radius: default 5
imagescale=4;     %resize image so it's bigger
framesperhr=5;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagefignum=plotfignum+1;
pos=get(event_obj,'Position');
frame=pos(1)*framesperhr;
if frame==immunoframe
    framestring='poststain';
else
    framestring=num2str(frame);
end
signal=num2str(pos(2),4);
output_txt={['X: ',framestring],['Y: ',signal]};
%%% get x and y coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tar=get(event_obj,'Target');
traceid=str2double(get(tar,'DisplayName'));  %returns cellid
allx=gatedtraces(:,frame,1)-jitters(frame,1); ally=gatedtraces(:,frame,2)-jitters(frame,2);
cx=int16(allx(traceid));
cy=int16(ally(traceid));
%%%%% load image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(imagefignum);  %makes image figure the current figure
clf;
mask=single(imread([maskdir,'nucedge_',framestring,'.tif']));
switch channel
    case 1
        filecode='CFP_';
        maxval=max(gatedtraces(traceid,:,4)./gatedtraces(traceid,:,3),[],2);
        bgperctile=50;
    case 2
        filecode='YFP_';
        maxval=max(gatedtraces(traceid,:,5),[],2);
        bgperctile=50;
    case 3
        filecode='TexasRed_';
        maxval=max(gatedtraces(traceid,:,6),[],2);
        bgperctile=50;
end
raw=single(imread([rawdir,filecode,framestring,'.tif']));
%%%%% crop image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(mask);
%minpix=median(prctile(raw,0));  %sets absolute intensities. default 1
%maxpix=median(prctile(raw,100)); %default 99
minx=max([cx-winrad 1]);
maxx=min([cx+winrad width]);
miny=max([cy-winrad 1]);
maxy=min([cy+winrad height]);
cropmask=mask(miny:maxy,minx:maxx);
cropraw=raw(miny:maxy,minx:maxx);
%%% calc background %%%
bgmask=imfill(cropmask,'holes');
bgraw=cropraw.*~bgmask;
bgraw(bgraw==0)=NaN;
cropbg=prctile(bgraw(:),bgperctile);
%cropraw=cropraw-prctile(cropraw(:),10);
cropraw=cropraw-cropbg;
cropraw(cropraw>maxval)=maxval;
%cropmask=mat2gray(imresize(cropmask,imagescale));
cropmask=mat2gray(imresize(cropraw,imagescale));
cropmask(:,:,2)=mat2gray(imresize(cropraw,imagescale));
%cropmask(:,:,2)=0;
cropmask(:,:,3)=0;
cropmask(:,:,1)=0;
image(cropmask);
%%%%% add cellid over cells in cropped image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsinimage=find(ally>miny & ally<maxy & allx>minx & allx<maxx);
for i=cellsinimage
    %text(double(int16(allx(i))-minx)*imagescale,double(int16(ally(i))-miny)*imagescale,num2str(gatedtraces(i,frame,end-1)),'horizontalalignment','center','color','r','fontsize',10,'fontweight','bold');
end
%%%%% return control to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(plotfignum);            %returns control to plot