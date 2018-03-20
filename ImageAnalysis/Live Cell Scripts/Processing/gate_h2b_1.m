function [signals,badtraces]=gate_h2b_1(tracedata,tracestats,noisethresh,sigchannel,minthreshold,minpos,maxthreshold,maxpos)
%hist(max(signals,[],2),100);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
rawsignal=tracedata(:,:,sigchannel);
% nucareas = tracedata(:,:,areachannel); 
% signals = rawsignal.*nucareas;
signals = rawsignal;
% for i=1:numtraces
%     realframes=find(~isnan(signals(i,:)));
%     signals(i,realframes)=smooth(signals(i,realframes));
% end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toolow=ones(numtraces,1)*NaN;
switch minpos
    case 0
        toolow=max(signals,[],2)<minthreshold;
    case 1
        for i=1:numtraces
            toolow(i)=signals(i,tracestats(i,1))<minthreshold;
        end
    case 2
        for i=1:numtraces
            toolow(i)=signals(i,tracestats(i,2))<minthreshold;
        end
end

toohigh=ones(numtraces,1)*NaN;
switch maxpos
    case 0
        toohigh=min(signals,[],2)>maxthreshold;
    case 1
        for i=1:numtraces
            toohigh(i)=signals(i,tracestats(i,1))>maxthreshold;
        end
    case 2
        for i=1:numtraces
            toohigh(i)=signals(i,tracestats(i,2))>maxthreshold;
        end
end



%%%% Gate Noisy Traces %%%%%
noise=ones(numtraces,1)*NaN;
for i=1:numtraces
%     maxdiff=max(diff(signals(i,tracestats(i,1):tracestats(i,2)),1));
    normdiff = diff(signals(i,tracestats(i,1):tracestats(i,2)),1)./signals(i,tracestats(i,1):(tracestats(i,2)-1));
    noise(i) = sum(normdiff<-0.6 | normdiff > 0.2);
    
%     figure(1)
%     plot(normdiff)
%     hold on
%     plot(signals(i,:))
%     hold off
%     set(gca,'YLim',[-1 1]);
%     line(get(gca,'XLim'),[.2 .2],'Color','k','LineStyle','--')
%     line(get(gca,'XLim'),[-.6 -.6],'Color','k','LineStyle','--')
%     title(num2str(noise(i)))
%     if ~isempty(maxdiff)
%         noise(i)=maxdiff;
%     end
end
noisy= noise>noisethresh;

badtraces = noisy | toolow | toohigh;
% signals = tracedata(:,:,sigchannel).*tracedata(:,:,areachannel); %return raw signal (not smoothened)
signals = tracedata(:,:,sigchannel);