saveopt = 0 ;
genenames = {'pparg','cebpa','fabp4'};
% genenames = {'pparg'};
time= [0,0.5,1,1.5,3,6];
% timemat = repmat(time,3,1);
timemat = time;
numgenes = numel(genenames);
datadir= 'Y:\michael\TF dynamics data\mRNA_decay_qPCR\Excel Files\';
timerange = 0:.1:6;       
for gene = 1:numgenes
    timevec = timemat(:);
    datafile = [datadir,genenames{gene},'-decay-analysis.xlsx'];
    datamat = xlsread(datafile,'Pfaffl','A13:M19');
    a = zeros(6,4);
    % Norm gene:rpl18//Condition: With DMI
    a(:,1)= datamat(3,1:6);
    %Norm gene:rpl0//Condition: With DMI
    a(:,2) = datamat(3,8:13);
    %Norm gene:rpl18//Condition: Without DMI
    a(:,3) = datamat(7,1:6);
    %Norm gene:rpl0//Condition: Without DMI
    a(:,4) = datamat(7,8:13);
    figure
    for ind = 1:2
        decayvals = a(:,ind);
        modelfit = fit(timevec,decayvals(:),'exp1');
        subtightplot(1,2,ind,[.1 .1])
        axis square
        plot_decay_curve(genenames{gene},[timevec,decayvals(:)],timerange,modelfit);
    end
    
end

