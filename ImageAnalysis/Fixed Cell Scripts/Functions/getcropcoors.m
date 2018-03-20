function cropcoordinates = getcropcoors(imsize,jitmatx,jitmaty)
% imsize in [height,width]

% round to nearest whole pixel
jitmatx = round(jitmatx);
jitmaty = round(jitmaty);

if jitmatx <0    
    x1A = 1; 
    x2A = imsize(2) + jitmatx;
    x1B = 1 - jitmatx;
    x2B = imsize(2);
elseif jitmatx > 0 
    x1A = 1 + jitmatx; 
    x2A = imsize(2);
    x1B = 1;
    x2B = imsize(2) - jitmatx;
elseif jitmatx == 0
    x1A = 1; 
    x2A = imsize(2); 
    x1B = 1;
    x2B = imsize(2);
end

if jitmaty <0    
    y1A = 1; 
    y2A = imsize(1) + jitmaty;
    y1B = 1 - jitmaty;
    y2B = imsize(1);
elseif jitmaty > 0 
    y1A = 1 + jitmaty; 
    y2A = imsize(1);
    y1B = 1;
    y2B = imsize(1) - jitmaty;
elseif jitmaty == 0
    y1A = 1; 
    y2A = imsize(1); 
    y1B = 1;
    y2B = imsize(1);
end

cropcoordinates = [y1A,y2A,x1A,x2A;y1B,y2B,x1B,x2B];
end