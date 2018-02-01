% Copyright (c) 2017 The Johns Hopkins University/Applie?d Physics Laboratory
% This software was developed at The Johns Hopkins University/Applied Physics 
% Laboratory (â€œJHU/APLâ€?) that is the author thereof under the â€œwork made for hireâ€? 
% provisions of the copyright law. Permission is hereby granted, free of charge, 
% to any person obtaining a copy of this software and associated documentation 
% (the â€œSoftwareâ€?), to use the Software without restriction, including without 
% limitation the rights to copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit others to do so, subject to 
% the following conditions:
% 1. This LICENSE AND DISCLAIMER, including the copyright notice, 
% shall be included in all copies of the Software, 
% including copies of substantial portions of the Software;
% 
% 2. JHU/APL assumes no obligation to provide support of any kind with regard to the Software. 
% This includes no obligation to provide assistance in using the Software nor to 
% provide updated versions of the Software; and
% 
% 3. THE SOFTWARE AND ITS DOCUMENTATION ARE PROVIDED AS IS AND WITHOUT ANY 
% EXPRESS OR IMPLIED WARRANTIES WHATSOEVER. ALL WARRANTIES INCLUDING, BUT NOT 
% LIMITED TO, PERFORMANCE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
% AND NONINFRINGEMENT ARE HEREBY DISCLAIMED. USERS ASSUME THE ENTIRE RISK AND 
% LIABILITY OF USING THE SOFTWARE. USERS ARE ADVISED TO TEST THE SOFTWARE 
% THOROUGHLY BEFORE RELYING ON IT. IN NO EVENT SHALL THE JOHNS HOPKINS UNIVERSITY 
% BE LIABLE FOR ANY DAMAGES WHATSOEVER, INCLUDING, WITHOUT LIMITATION, ANY 
% LOST PROFITS, LOST SAVINGS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES, 
% ARISING OUT OF THE USE OR INABILITY TO USE THE SOFTWARE.â€?% 
% It generates a segmentation based only on connectivity of strong edges
% edge estimation => thresholding => connected components => bbox creating from cc
% 
% Version History:
%	09/30/2016  Marc Bosch   Initial version
function [blobIndIm2,boxes1,BW,lvep]=contour_screener_stal(img)
[M,N]=size(img);

[M,N,ch]=size(img);
[Lvep,lvep,tr,de] = structensorsvd(img(:,:,1),0.1);
[Lvep2,lvep2,tr2,de2] = structensorsvd(img(:,:,2),0.1);
[Lvep3,lvep3,tr3,de3] = structensorsvd(img(:,:,3),0.1);
lvep_array(:,:,1)=lvep;
lvep_array(:,:,2)=lvep2;
lvep_array(:,:,3)=lvep3;
lvep=max(lvep_array,[],3);
curdir=pwd;
addpath(genpath('./3rdparty/StrictiredEdgeDetection'))
setParametersSED;
lvep = 2*myrenormalize(edgesDetect(uint8(img),model))+myrenormalize(lvep);
cd(curdir)
gim = zeros(M,N);
gim(3:M-2,3:N-2)=sqrt(lvep(3:M-2,3:N-2));
BW = gim>(mean(gim(:))+2.5*std(gim(:)));
% h=ones(301,301);
% BW=gim>(imfilter(gim,h)./sum(h(:)))+0.75*stdfilt(gim,h);
BW=imresize(BW,[size(img(:,:,1),1),size(img(:,:,1),2)]);
% BW=medfilt2(BW,[5 5]);
se1=zeros(5,5);se1(3,:)=1;se1(:,3)=1; %it should break thicker edges
se2=zeros(3,3);se2(2,:)=1;se2(:,2)=1; %it should break edges but keep thicker ones
% BW = imerode(BW,se);
BW = bwmorph(BW,'thin',Inf);
CC=bwconncomp(BW,8);
Nsgms=1;
blobIndIm2=zeros(size(img,1),size(img,2));
for jji=1:CC.NumObjects
    pixpos=CC.PixelIdxList{jji};

        blobIndIm2(pixpos)=Nsgms+1;
        Nsgms=Nsgms+1;
    
end



    props = regionprops(blobIndIm2);   
    props = struct2cell(props);
    props = cell2mat(props(3,:));
    bboxtemp = reshape(props, [4,length(props)/4])';
    boxes1=[round(bboxtemp(:,2)) round(bboxtemp(:,1)) ...
        round(bboxtemp(:,2)+bboxtemp(:,4)) round(bboxtemp(:,3)+bboxtemp(:,1))];
