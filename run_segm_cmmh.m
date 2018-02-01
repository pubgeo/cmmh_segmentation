% This file shows how to run the algorithm introduced in:
% M.Bosch, C. Gifford, A.Dress, C. Lau, J.Skibbo, G.Christie
% "Improved image segmentation via cost minimization of multiple
% hypotheses", BMVC 17.
% If this code is used please cite above.
% This version only supports FH kernel as described in the paper.
% The software only comes pre-compiled for windows. For linux/mac mex compilation is
% needed before running.

% Copyright (c) 2017 The Johns Hopkins University/Applie?d Physics Laboratory
% This software was developed at The Johns Hopkins University/Applied Physics 
% Laboratory ("JHU/APL") that is the author thereof under the "work made for hire"
% provisions of the copyright law. Permission is hereby granted, free of charge, 
% to any person obtaining a copy of this software and associated documentation 
% (the "Software"), to use the Software without restriction, including without
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
% ARISING OUT OF THE USE OR INABILITY TO USE THE SOFTWARE.
addpath(genpath('./'))
param = 150; % tunable parameter that controls degree of over/under segmentation of each hypothesis in FH kernel
Numhyp = 30; % number of hypotheses, ideally we want a number that covers the volume space and variation in hypotheses
solver=0;%turn to 1 for Belief propagation minimizer (slower, a bit better)
im=double((imread('./BSDS300/images/train/253036.jpg')));
[segmentedmask,coloredsegment,imgsegment] = segm_cmmh(im,param,Numhyp,solver);
figure;imshow(uint8(imgsegment))