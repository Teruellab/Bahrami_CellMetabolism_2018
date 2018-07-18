%%
clear all
close all

load 'weights' WeightsTable

% figure,plot(WeightsTable(1:6,2),WeightsTable(1:6,1),'ro')
% hold on,plot(WeightsTable(7:12,2),WeightsTable(7:12,1),'bo')
% figure,plot(WeightsTable(1:6,3),WeightsTable(1:6,1),'ro')
% hold on,plot(WeightsTable(7:12,3),WeightsTable(7:12,1),'bo')
figure,plot(WeightsTable(1:6,3),WeightsTable(1:6,2),'ro')
hold on,plot(WeightsTable(7:12,3),WeightsTable(7:12,2),'bo')
title('Correlation between visceral and subq fat, red Cort, blue sham')

% E=visceral, gonadal, epidydma
% B=subq, inguinal, leg fat
WeightCB(1:6)=[621.2 397.1 670.5 686.9 375.7 340.4];
WeightSB(1:6)=[226.4 208.9 202 212.5 328.9 330.7];
WeightCE(1:6)=[1357.2 840.4 1210.8 1366.6 1000.3 646];
WeightSE(1:6)=[459.5 385 474.6 401.8 594.5 597.5];

WSB(1:15)=WeightSB(1);  %Lookup table to show which mouse the images(slicesxsites) belong to
WSB(16:30)=WeightSB(2);
WSB(31:45)=WeightSB(3);
WSB(46:59)=WeightSB(4);
WSB(60:74)=WeightSB(5);
WSB(75:89)=WeightSB(6);

WCB(1:15)=WeightCB(1);  %Lookup table to show which mouse the images(slicesxsites) belong to
WCB(16:30)=WeightCB(2);
WCB(31:45)=WeightCB(3);
WCB(46:60)=WeightCB(4);
WCB(61:75)=WeightCB(5);
WCB(76:90)=WeightCB(6);

WSE(1:15)=WeightSE(1);  %Lookup table to show which mouse the images(slicesxsites) belong to
WSE(16:30)=WeightSE(2);
WSE(31:45)=WeightSE(3);
WSE(46:60)=WeightSE(4);
WSE(61:75)=WeightSE(5);
WSE(76:90)=WeightSE(6);

WCE(1:15)=WeightCE(1);  %Lookup table to show which mouse the images(slicesxsites) belong to
WCE(16:30)=WeightCE(2);
WCE(31:45)=WeightCE(3);
WCE(46:60)=WeightCE(4);
WCE(61:75)=WeightCE(5);
WCE(76:90)=WeightCE(6);

imagedir = '/Volumes/labdata/mary/';
experimentdir = [imagedir,'171220_FatHistology/'];
datadir = [experimentdir,'Data/20171221/'];
imagedir1 = [datadir,'ShamB/'];
Files=dir(imagedir1);
for i=3:length(Files)
    load([imagedir1 Files(i).name]);
    s1=size(finaldata);
    clear d
    for num1=1:s1
        d(num1)=finaldata(num1).Area;
    end
    VolSB(i-2)=1.28^3 * 4/3*sqrt(pi)*median(d.^1.5); %median cell volume in each image
end

imagedir = '/Volumes/labdata/mary/';
experimentdir = [imagedir,'171220_FatHistology/'];
datadir = [experimentdir,'Data/20171221/'];
imagedir1 = [datadir,'CortB/'];
Files=dir(imagedir1);
for i=3:length(Files)
    load([imagedir1 Files(i).name]);
    s1=size(finaldata);
    clear d
    for num1=1:s1
        d(num1)=finaldata(num1).Area;
    end
    VolCB(i-2)=1.28^3 * 4/3*sqrt(pi)*median(d.^1.5);
end

figure,hold on
for i=1:88
    plot(VolSB(i),WSB(i),'b.')
end
for i=1:90
    plot(VolCB(i),WCB(i),'r.')
end
for i=1:length(VolSB)
    NumSB(i)=10^9*WSB(i)/VolSB(i);
end
for i=1:length(VolCB)
    NumCB(i)=10^9*WCB(i)/VolCB(i);
end

['SB, volume: ' num2str(round(mean(VolSB))) ' +- ' num2str(round(std(VolSB)/sqrt(88)))]
['CB, volume: ' num2str(round(mean(VolCB))) ' +- ' num2str(round(std(VolCB)/sqrt(90)))]

['SB, #cells: ' num2str(round(mean(NumSB))) ' +- ' num2str(round(std(NumSB)/sqrt(88)))]
['CB, #cells: ' num2str(round(mean(NumCB))) ' +- ' num2str(round(std(NumCB)/sqrt(90)))]

%%
imagedir = '/Volumes/labdata/mary/';
experimentdir = [imagedir,'171220_FatHistology/'];
datadir = [experimentdir,'Data/20171221/'];
imagedir1 = [datadir,'ShamE/'];
Files=dir(imagedir1);
for i=3:length(Files)
    load([imagedir1 Files(i).name]);
    s1=size(finaldata);
    clear d
    for num1=1:s1
        d(num1)=finaldata(num1).Area;
    end
    VolSE(i-2)=1.28^3 * 4/3*sqrt(pi)*median(d.^1.5);
end

imagedir = '/Volumes/labdata/mary/';
experimentdir = [imagedir,'171220_FatHistology/'];
datadir = [experimentdir,'Data/20171221/'];
imagedir1 = [datadir,'CortE/'];
Files=dir(imagedir1);
for i=3:length(Files)
    load([imagedir1 Files(i).name]);
    s1=size(finaldata);
    clear d
    for num1=1:s1
        d(num1)=finaldata(num1).Area;
    end
    VolCE(i-2)=1.28^3 * 4/3*sqrt(pi)*median(d.^1.5);
end

for i=1:90
    NumSE(i)=10^9*WSE(i)/VolSE(i);
    NumCE(i)=10^9*WCE(i)/VolCE(i);
end

['SE, volume: ' num2str(round(mean(VolSE))) ' +- ' num2str(round(std(VolSE)/sqrt(89)))]
['CE, volume: ' num2str(round(mean(VolCE))) ' +- ' num2str(round(std(VolCE)/sqrt(90)))]

['SE, #cells: ' num2str(round(mean(NumSE))) ' +- ' num2str(round(std(NumSE)/sqrt(90)))]
['CE, #cells: ' num2str(round(mean(NumCE))) ' +- ' num2str(round(std(NumCE)/sqrt(90)))]

figure,hold on
for i=1:90
    plot(VolSE(i),WSE(i),'b.')
end
for i=1:90
    plot(VolCE(i),WCE(i),'r.')
end
