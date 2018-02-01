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
function [score] = segminitscore_stal(blobIndIm2,lvep,direction)
[M,N]=size(lvep);


if strcmp(direction,'bottom')==1

    blobIndIm2_ref=blobIndIm2;
    blobIndIm2_cur=blobIndIm2;
    blobIndIm2_cur(1:M-1,:)=blobIndIm2(2:M,:);
    diflabel=(blobIndIm2_ref~=blobIndIm2_cur);        
   
elseif strcmp(direction,'topright')==1   

    blobIndIm2_ref=blobIndIm2;
    blobIndIm2_cur=blobIndIm2;
    blobIndIm2_cur(2:M,1:N-1)=blobIndIm2(1:M-1,2:N);
    diflabel=(blobIndIm2_ref~=blobIndIm2_cur);    

elseif strcmp(direction,'bottomright')==1    

    blobIndIm2_ref=blobIndIm2;
    blobIndIm2_cur=blobIndIm2;
    blobIndIm2_cur(1:M-1,1:N-1)=blobIndIm2(2:M,2:N);
    diflabel=(blobIndIm2_ref~=blobIndIm2_cur);    

elseif strcmp(direction,'right')==1

    blobIndIm2_ref=blobIndIm2;
    blobIndIm2_cur=blobIndIm2;
    blobIndIm2_cur(:,1:N-1)=blobIndIm2(:,2:N);
    diflabel=(blobIndIm2_ref~=blobIndIm2_cur);    

end


lvepfilt=lvep.*(blobIndIm2>0);
score=50*diflabel.*max(lvepfilt(:))+0*(max(diflabel(:))-diflabel).*lvepfilt+1;%sum(scorearray(:,:,[4]),3).*(max(diflabelfilt(:))-diflabelfilt);
