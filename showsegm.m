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
function [h,blobIndImdisplay]=showsegm(blobIndIm)

blobIndImdisplay1=reshape(blobIndIm,[size(blobIndIm,1)*size(blobIndIm,2) 1]);
blobIndImdisplay2=reshape(blobIndIm,[size(blobIndIm,1)*size(blobIndIm,2) 1]);
blobIndImdisplay3=reshape(blobIndIm,[size(blobIndIm,1)*size(blobIndIm,2) 1]);
labels=unique(blobIndIm);
for ij=1:length(labels)
    color(ij,1)=min(255,abs(125*randn(1)));
    color(ij,2)=min(255,abs(125*randn(1)));
    color(ij,3)=min(255,abs(125*randn(1)));
    po=find(blobIndIm==labels(ij));
    blobIndImdisplay1(po) = color(ij,1);
    blobIndImdisplay2(po) = color(ij,2);
    blobIndImdisplay3(po) = color(ij,3);
end
blobIndImdisplay(:,:,1)=reshape(blobIndImdisplay1,[size(blobIndIm,1) size(blobIndIm,2)]);
blobIndImdisplay(:,:,2)=reshape(blobIndImdisplay2,[size(blobIndIm,1) size(blobIndIm,2)]);
blobIndImdisplay(:,:,3)=reshape(blobIndImdisplay3,[size(blobIndIm,1) size(blobIndIm,2)]);
h=imshu8(blobIndImdisplay);
