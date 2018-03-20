function [signals,badtraces]=gate_h2b_2(tracedata,tracestats,noisethresh,sigchannel,gatechannel,mitosistimes)
%hist(max(signals,[],2),100);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
rawsignal=tracedata(:,:,gatechannel);
% nucareas = tracedata(:,:,areachannel); 
% signals = rawsignal.*nucareas;
smoothsignals = ones(size(rawsignal))*nan;
for i=1:numtraces
    if tracestats(i,3)< 10
       smoothsignals(i,:)=medsmoothignorenans(rawsignal(i,:),9);
    else
        smoothsignals(i,:) = smoothignorenans_butter(rawsignal(i,:),3,0.1);
    end
end

nummitosis = cellfun(@(x) numel(x),mitosistimes);
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Gate Super Noisy Traces %%%%%
nanvec = ones(numtraces,1)*NaN;
noise      = nanvec;
jumps      = nanvec;
mitosisnoh2b = nanvec;
h2bnomitosis = nanvec;
for i=1:numtraces
    badmitosiscall = 0;
    % noisy
%     maxdiff=max(diff(signals(i,tracestats(i,1):tracestats(i,2)),1));
    normdiff = [0 diff(rawsignal(i,:),1)./rawsignal(i,1:end-1)];
    noise(i) = sum(normdiff<-0.1 | normdiff > 0.1);
    
    % jumpy
    smoothdiff = [0 diff(smoothsignals(i,:),1)./smoothsignals(i,1:end-1)];
    jumps(i) = sum(smoothdiff>0.15 | smoothdiff<-0.5);
    
%     % mitosis call without h2b drop
%     if nummitosis(i) > 1
%         h2bsplitcheck = zeros(nummitosis(i)-1,1);
%         for m = 2:nummitosis(i)
%             mtime = mitosistimes{i}(m);
%             if mtime<6
%                 h2bdiff = normdiff(1:mtime);
%             else
%                 h2bdiff = normdiff(mtime-5:mtime);
%             end
%             h2bsplitcheck(m-1) = any(h2bdiff<-0.4 & h2bdiff>-0.6);
%         end
%         mitosiscall = sum(h2bsplitcheck) == (nummitosis(i)-1);
%         if mitosiscall == 0 
%             mitosisnoh2b(i) = 1;
%         else 
%             mitosisnoh2b(i) = 0;
%         end
%     else
%         mitosisnoh2b(i) = 0;
%     end
    
    
%     % large h2b decreases without mitosis call
%     h2bhalves = find(smoothdiff<-0.2);
%     numh2bhalves  = numel(h2bhalves);
%     if numh2bhalves > 0 % if halves exist, check to see whether there is mitosis call nearby
%         framediff = zeros(numh2bhalves,1);
%         for h = 1:numh2bhalves
%             %does the h2b drop have a corresponding mitosis call within 5
%             %frames after?
%             framediff(h) = any((mitosistimes{i} - h2bhalves(h))<6 & (mitosistimes{i} - h2bhalves(h))>0);         
%         end
%         if any(framediff==0)
%             h2bnomitosis(i) = 1;
%         else 
%             h2bnomitosis(i) = 0;
%         end
%     else
%         h2bnomitosis(i) = 0; %if no large h2b decrease exist, then there is remains false
%     end
%     
    
end
noisy= noise>noisethresh;
jumpy = jumps>0;
% mismatchmitosis = mitosisnoh2b>0;
% mismatchh2b = a;

badtraces = noisy | jumpy ;
% signals = tracedata(:,:,sigchannel).*tracedata(:,:,areachannel); %return raw signal (not smoothened)
signals = tracedata(:,:,sigchannel);