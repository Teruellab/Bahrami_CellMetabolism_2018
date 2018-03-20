%% Illumination Bias estimation

% Directory Name
ExperimentPath ='D:\Experiments\20180127-e2f7-e2f8-ab-test\';

nucName = 'Hoechst';
sigName = {'pparg','e2f','prb'}; %,'mChy','Cy5'
numSigs = numel(sigName);
 
rows = [2:3];
cols = [2:10];
sites = 1:4; numsites = numel(sites);
plates = 1; numplates = numel(plates);
nucr = 8; % nuclear radius; nuclear radius = 12 for bin 1
bin = 2;

%% Don't touch
for plate = 1:numplates
    for site = 1:numsites
        Nuc_Bias_Estimation_2(ExperimentPath,nucr,nucName,rows,cols,bin,plates(plate),site)
    end
end

for sign = 1:numSigs
    for plate = 1:numplates
        for site = 1:numsites
            Sig_Bias_Estimation_2(ExperimentPath,nucr,nucName,sigName{sign},rows,cols,bin,plates(plate),site)
        end
    end
end
