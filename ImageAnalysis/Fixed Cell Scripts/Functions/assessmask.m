function assessmask(raw,mask)
[height,width]=size(raw);
ImRGB=zeros(height,width,3);
ImRGB(:,:,1)=imadjust(mat2gray(raw));
ImRGB(:,:,2)=bwperim(mask)+imadjust(mat2gray(raw));
ImRGB(:,:,3)=imadjust(mat2gray(raw));
imshow(ImRGB,'Border','tight');
end