function mask=segmentdeflections_bwboundaries(mask,nucr,debrisarea)
[B,L]=bwboundaries(mask,'noholes');
obnum=max(L(:));
bordermask=zeros(size(mask));
for ci=1:obnum
    orderedset=B{ci};
%     orderedset = B{943};
    % If using bwboundaries, reverse the order.
%     orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)];
    orderedset=[orderedset(end:-1:2,2) orderedset(end:-1:2,1)];
    
    bordermask=splitdeflections_4_bwboundaries(orderedset,bordermask,nucr);
%     bordermask=splitdeflections_5_bwboundaries(orderedset,bordermask,nucr);
end
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,debrisarea);
end

% nuc_info=struct2cell(regionprops(L,'Area','Centroid')');
% nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
% testim = insertText(double(mask),nuc_center,1:obnum,'AnchorPoint','Center','BoxOpacity',0,'TextColor','r','FontSize',10,'TextColor','green');
% % testim(:,:,2) = 0;
% % testim(:,:,3) = 0;
% imshow(testim)