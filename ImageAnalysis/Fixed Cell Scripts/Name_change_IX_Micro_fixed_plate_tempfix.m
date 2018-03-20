%% Name_change_IX_Micro_fixed_plate.m


%% Init
clear
close all
clc



directoryInput = {'C:\Users\Teruel Lab\Desktop\160601-KK\CRISPRi Knockdown\2016-06-02\12181\'};
directoryOutput = {'Z:\kyle\CRISPRi\160601_CRISPRi_Knockdown\Raw\'};
channelNames = {'Hoechst','488','mCherry','647'}; 

numdirs = numel(directoryInput);


time1=tic;  %Starts Timer
%% Load Directory with images

for dirnum = 1:numdirs

    dir= directoryInput{dirnum};

    alldir=getSubdirectories(dir); %Creates list of all the subfolders

 
    disp(dirnum) %Displays what timepoint you are on
    time2 = tic; %Starts timer to record time of each timepoint
    filenames = getFilenames(dir); %Returns list of all file names inside the folder
    filenames = filenames(boolRegExp(filenames,'.tif') & ~boolRegExp(filenames,'thumb')); %removes accessory files and thumbnails from the file list
    
    
    parfor j=1:length(filenames)
        % Designate Channels
        %Change this to match the correct channels used in the experiment  %e.g. {'Channel1','Channel2','Channel3'}  
        channels = channelNames;
        
        % target directory
        dirlist2 = directoryOutput;
        
        % extract wellname, site, and channel from original file name
        wellname=char(getTokens(filenames{j},'^[^_]+_([A-Z][0-9][0-9])'));  %Finds well name in format wellColumn (eg B04)
        row_column=wellNameToRowColumn(wellname); %Converts well name to format row_column (eg 2_4)
        site=char(getTokens(filenames{j},'^[^_]+_[A-Z][0-9][0-9]_s([0-9])')); % for multiple sites %Finds site image was taken (eg 1)
%         site = 1; % for 1 site
       
        channel = char(getTokens(filenames{j},'^[^_]+_[A-Z][0-9][0-9]_s[0-9]_w([1-9])'));  % For multiple Sites%Finds channel number (eg 1)
        
        if isempty(channel)
            site=char(getTokens(filenames{j},'^[^_]+_[A-Z][0-9][0-9]_s([0-9][0-9])'));
            channel = char(getTokens(filenames{j},'^[^_]+_[A-Z][0-9][0-9]_s[0-9][0-9]_w([1-9])')); 
        end
        
        
        
        
%         channel=char(getTokens(filenames{j},'^[^_]+_[A-Z][0-9][0-9]_w([1-9])')); %for one site

        % YOU HAVE TO CHANGE CHANNEL AND SITE IF YOU HAVE ONE SITE VERSUS
        % MULTIPLE SITE


         newFileName=[row_column,'_',num2str(site),'_',channels{str2double(channel)},'_',num2str(dirnum),'.tif'];   
         oldFileName=[filenames{j}];        %Designates old file name to be used in old file path
         
         %make sure file names actually exist before saving
         test = isempty(newFileName);
         if test == 0
%          copyfile([dir,oldFileName],[dir2,wellname,'\site_',num2str(site),'\',newFileName],'f'); %Moves file from old folder to new folder and changes the name. Use this line if you have multiple sites. 
         copyfile([dir,oldFileName],[dirlist2{1},newFileName],'f');
         else
             continue
         end
         %copyfile([dir,'\',oldFileName],[dir2,wellname,'\',newFileName],'f'); % for one site
         
    end
    toc(time2) %Displays time it takes to rename all the files in each Time point

  %Displays total time to rename all files

end
toc(time1)