% Copyright (c) 2017 The Johns Hopkins University/Applie?d Physics Laboratory
% This software was developed at The Johns Hopkins University/Applied Physics 
% Laboratory (“JHU/APL”) that is the author thereof under the “work made for hire” 
% provisions of the copyright law. Permission is hereby granted, free of charge, 
% to any person obtaining a copy of this software and associated documentation 
% (the “Software”), to use the Software without restriction, including without 
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
% ARISING OUT OF THE USE OR INABILITY TO USE THE SOFTWARE.”
%It generates the initial guide to the affinity map for the FH segmentation
%(for all 4 neighbors right,bottom,bottomright,topright. There are two
%methods: (1) structure tensor eigen value (recommended), (2) difference of
%gaussians
% Modified Version History:
%	10/06/2016  Marc Bosch   Initial version
function [score1,score4,score2,score3,BW,blobIndIm2,boxes] = FH_guidescores_stal(imgstructure_1,method,screener_method1_on)
%method - 1 is using the small eigenvalue of the structure tensor 
%method - 2 uses the laplacian of gaussian 
[M,N,ch]=size(imgstructure_1);
if method==1
      

      [blobIndIm2,boxes,BW,lvep]=contour_screener_stal(imgstructure_1);
      blobIndIm2=BW;
%      showsegm(blobIndIm2),title('initial guide')
      lvep=(lvep-min(lvep(:)))/(max(lvep(:))-min(lvep(:)));
else
      dog=imfilter(imgstructure_1,filterDog2d(9,9,1,0),'same');
      partsmask=(dog>mean(dog(:))+std(dog(:)));
      BW=imresize(partsmask,[size(imgstructure_1(:,:,1),1),size(imgstructure_1(:,:,1),2)]);
      BW=medfilt2(BW,[5 5]);
      
      %guided filter
%       [BW]=bilatmedian(BW,imgstructure_1,7);
      CC = bwconncomp(BW, 8);
     blobIndIm2=zeros(M,N);
      Nsgms=1;
      for jji=2:CC.NumObjects
          pixpos=CC.PixelIdxList{jji};
          blobIndIm2(pixpos)=Nsgms+1;
          Nsgms=Nsgms+1;
      end
      labels=unique(blobIndIm2);
      cn=1;
      NN=length(unique(blobIndIm2));
      for n = 1:NN        
          pp=find(blobIndIm2==labels(n));
          [rw,cw]=ind2sub(size(imgstructure_1),pp);
          if length(rw)>0
              minr=min(rw);
              maxr=max(rw);
              minc=min(cw);
              maxc=max(cw);        
              boxes(cn,:)=[minr minc maxr maxc];        
              cn=cn+1;
          end
      end

      [Lvep,lvep,tr,de] = structensorsvd(imgstructure_1,0.1);
      lvep=(lvep-min(lvep(:)))/(max(lvep(:))-min(lvep(:)));
%       showsegm(blobIndIm2),title('initial guide2')
end

 % for each of the for edges get the scores for the FH framework
 if screener_method1_on==1
  score1 = segminitscore_stal(blobIndIm2,lvep,'right');
  score2 = segminitscore_stal(blobIndIm2,lvep,'bottomright');
  score3 = segminitscore_stal(blobIndIm2,lvep,'topright');
  score4 = segminitscore_stal(blobIndIm2,lvep,'bottom');
 else
  score1=[];score2=[];score3=[];score4=[];
 end