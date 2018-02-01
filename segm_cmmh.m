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
function [segmentedmask,coloredsegment,imgsegment] = segm_cmmh(im,param,Numhyp,solver)

im=double(im);
if size(im,3)==1
    'RGB color image is expected'
    imtemp=im;
    im(:,:,1)=imtemp;
    im(:,:,2)=imtemp;
    im(:,:,3)=imtemp;
end
% clear all
close all
addpath(genpath('./3rdparty/'));
pathimages='./BSDS300/images/combined/';
imgnames=dir([pathimages '*.jpg']);

 averageBoundaryError = [];
 averageRI = [];
 averageVOI = [];
 averageGCE = [];
 averageCOV = [];

 
%param = 150; % tunable parameter that controls degree of over/under segmentation of FH kernel
range = Numhyp;            
iterascale = round(exp((1:range)/range).^2*(range/2));
%iteras range should produce from highly undersegmentation to 
%oversegmented hypothesis to cover full range.
iteras = param*iterascale+1;

%solver=0;%turn to 1 for Belief propagation minimizer,    



    clearvars -except iteras solver i im averageBoundaryError averageRI averageVOI averageGCE averageCOV imgnames pathimages
            
    imorig=double(rgb2gray(uint8(im)));
    [M,N,ch]=size(im);
%     gradx=derivative5(imorig,'x');
%     grady=derivative5(imorig,'y');
%     mag=abs(sqrt(gradx.^2+grady.^2));
    
    %%%%%%%%%%%%PRE-PROCESSING%%%%%%%%%%%%%%%%%
     fprintf('pre-processing\n');
     lab_img = colorspace('Lab<-', im);
  
    iiim(:,:,1)=myrenormalize(lab_img(:,:,1));
    iiim(:,:,2)=myrenormalize(lab_img(:,:,2));
    iiim(:,:,3)=myrenormalize(lab_img(:,:,3));
    %denoiser
    iminput(:,:,1)=FGS(iiim(:,:,1),0.01,30*30,imorig,4,4);
    iminput(:,:,2)=FGS(iiim(:,:,2),0.01,30*30,imorig,4,4);
    iminput(:,:,3)=FGS(iiim(:,:,3),0.01,30*30,imorig,4,4);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%COMPUTE affinity WEIGHTS FOR the edges of the graph in the FH kernel%%%%
    [score_right,score_bottom,score_bottomright,score_topright,BW,~,boxesinit] = ...
                FH_guidescores_stal(iiim,1,1);
     scores=(score_right+score_bottom+score_bottomright+score_topright);
     score_right=scores;
     score_bottom=scores;
     score_bottomright=scores;
     score_topright=scores;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    %%%%% Multi-hypothesis pool Creation%%%%%%
    fprintf('Creation of multi-hypothesis\n');
    Niter=length(iteras);
    blobIndImaccum=zeros(M,N,Niter);
    accumean=zeros(M,N,Niter,3);
    accudifmean=zeros(M,N,Niter,3);
    accusize=zeros(M,N,Niter);
    accustd=zeros(M,N,Niter);
    accumedge=zeros(M,N,Niter);
    features=im2single(iminput);
    %right-edge affinity scale
    features(:,:,4) = score_right;%sum(score1,3);
    %bottomright-edge affinity scale
    features(:,:,5) = score_bottomright;%sum(score2,3);
    %topright-edge affinity scale
    features(:,:,6) = score_topright;%sum(score3,3);
    %bottom-edge affinity scale
    features(:,:,7) = score_bottom;%sum(score4,3);

    %if we want to enforce the superpixels to a size
    %similar to the initial regions use convarea
    %otherwise set to all 1's
    features(:,:,8) = ones(size(features,1),size(features,2));%convarea;
    features(:,:,9) = zeros(size(features,1),size(features,2));
    
    featim1=features(:,:,1);
    featim2=features(:,:,2);
    featim3=features(:,:,3);
    fprintf('hypothesis:\n');
    for iter=1:length(iteras)
        
        fprintf('%d,',iter);
        %SEGMENTER
        [blobIndIm blobBoxes sizearray meanarray stdarray] = mexFelzenSegmentIndexFeatures ...
                    (features, ...
                   iteras(iter), 100, 0, true); 

       accusize(:,:,iter) = round(sizearray);
       accumean(:,:,iter,1) = meanarray;
       accustd(:,:,iter) = round(stdarray);
       uni=unique(blobIndIm);
       pim1=zeros(M,N);
       pim2=zeros(M,N);
       pim3=zeros(M,N);
       pimm1=zeros(M,N);
       pimm2=zeros(M,N);
       pimm3=zeros(M,N);
       for ijo=1:length(uni)
            p=find(blobIndIm==uni(ijo));                    
            pim1(p)=sum(abs(featim1(p)-mean(featim1(p))))/length(p);
            pim2(p)=sum(abs(featim2(p)-mean(featim2(p))))/length(p);
            pim3(p)=sum(abs(featim3(p)-mean(featim3(p))))/length(p);
            pimm1(p)=mean(featim1(p));
            pimm2(p)=mean(featim2(p));
            pimm3(p)=mean(featim3(p));
       end
       accudifmean(:,:,iter,1) = pim1;
       accudifmean(:,:,iter,2) = pim2;
       accudifmean(:,:,iter,3) = pim3;
       accumean(:,:,iter,1) = pimm1;
       accumean(:,:,iter,2) = pimm2;
       accumean(:,:,iter,3) = pimm3;
       %gradient
   %    [Ix, Iy] = derivative5(accumean(:,:,iter),'x','y'); % Get derivatives.
   %    canim(:,:,iter) = sqrt(Ix.*Ix + Iy.*Iy);   % Gradient magnitude.

       offsetNsemgs=0;
       if iter>1
        offsetNsemgs=max(blobIndImaccum(:));
       end
       %blobIndIm = cleanupregions(reshape(blobIndIm,M,N), 100, 4);
       blobIndImaccum(:,:,iter) = blobIndIm+offsetNsemgs;

       %measure difference with scores
       [gbx,gby]=imgradientxy(blobIndIm);%derivative5(blobIndIm,'x','y');
       magg=sqrt(gbx.^2+gby.^2);
       bww=(magg>(mean(magg(:))));
       accumedge(:,:,iter) = 255*bwmorph(bww,'thin',Inf);
%                   showsegm(blobIndIm)
       clear pim1 pim2 pim3 pimm1 pimm2 pimm3
    end
    %%%%natural edges map for the edge cost%%%%
    edgeref=uint8(myrenormalize(score_right+score_bottomright+score_topright+score_bottom));
    edgeref=255*double(bwmorph(edgeref,'thin',Inf));
    
    %%%%%%%%%%select best hypothesis for each pixel%%%%%%%%%%%%
    fprintf('\nCost minimization\n');
    c = select_best_iteration_stal(Niter,iminput,accudifmean,accumean, ...
        blobIndImaccum,accumedge,edgeref,solver);
    fprintf('Post-process\n');
    %assign attributes to each segment (size and intensity average)
    blobIndImfinal1 = zeros(size(accumean,1),size(accumean,2));
    segmmeanattrib = zeros(size(accumean,1),size(accumean,2));
    segszattrib = zeros(size(accumean,1),size(accumean,2));
    contrib=zeros(Niter,1);
    for iii=1:Niter-1
        p=find(c(:)==iii);
        contrib(iii)=length(p);
        blobtemp=blobIndImaccum(:,:,iii);
        [rr,cc]=find(c==iii);
        p=sub2ind(size(blobtemp), rr, cc);
        blobIndImfinal1(p)=blobtemp(p); 
        accumeantemp=accumean(:,:,iii,1);
        accusizettemp=accusize(:,:,iii);
        segszattrib(p)=accusizettemp(p);
        segmmeanattrib(p)=accumeantemp(p);
    end
    blobIndIm=blobIndImfinal1;
    %merge segments
    oldlab=unique(blobIndIm);
    Nclusters=length(oldlab);
    neighbours=zeros(Nclusters);
    %create consecutive unique segment ids for easier management, copy segment attributes     
    newlab=1:length(oldlab);
    blobIndImtemp=blobIndIm;
    segmattribvec=zeros(length(newlab),2);
    for is=1:length(oldlab)
        p=find(blobIndImtemp==oldlab(is));                
        blobIndIm(p)=is; 

%                 [length(p) is]
        segmattribvec(is,1) = segszattrib(p(1));  
        segmattribvec(is,5) = segmmeanattrib(p(1));

        %in case we want to use other attributes for merging
        %segments
        segmattribvec(is,3)=0;
        segmattribvec(is,2)=0;
        segmattribvec(is,4)=0;
        blobIndImtemp(p)=segmattribvec(is,2);
    end    
    %merge segments with neighboring segments if similar attributes
    gradsegm=imgradient(blobIndIm);
        for ij=2:size(blobIndIm,1)-1
            for ji=2:size(blobIndIm,2)-1
                if gradsegm(ij,ji)>0
                    ne1=blobIndIm(ij+1,ji);
                    ne2=blobIndIm(ij,ji+1);
                    ne3=blobIndIm(ij-1,ji);
                    ne4=blobIndIm(ij,ji-1);
                    neighbours(blobIndIm(ij,ji),ne1)=1;
                    neighbours(blobIndIm(ij,ji),ne2)=1;
                    neighbours(blobIndIm(ij,ji),ne3)=1;
                    neighbours(blobIndIm(ij,ji),ne4)=1;
                end
            end
        end
        blobIndIm1=blobIndIm;
        blobIndIm2=zeros(M,N);
        tic;
        %merge equivalent segments from different hypothesis based on
        %attributes and update neighboring segment map
        for ih=1:Nclusters
            neighbours(ih,ih)=0;
            %check size of current box
            idx=find(neighbours(ih,:)==1);
%             cn=1;        
            for ij=1:length(idx)
                maxN=max(segmattribvec(ih,1),segmattribvec(idx(ij),1));
                minN=min(segmattribvec(ih,1),segmattribvec(idx(ij),1));
                if (maxN/minN<1.1 && minN/maxN>0.9) && ...
                        (abs(segmattribvec(ih,5)-segmattribvec(idx(ij),5))<3)
                    [p] = find(blobIndIm1==ih);
                    blobIndIm1(p) = idx(ij);
                    neighbours(ih,idx)=0;
                    neighbours(idx,ih)=0;
                    neighbours(idx,idx(ij))=1;
                    neighbours(idx(ij),idx)=1;
                   % [i ij length(idx)]
                    break;
                end
            end
        end
        %update the segment mask with consecutive IDs and neighboring map
        [Nclusters,neighbours,blobIndIm1,segmattribvec]=UpdateSegmentMap_stal(blobIndIm1,segmattribvec);
        blobIndIm=blobIndIm1; 
        
        
        
         [Nclusters,neighbours,blobIndIm1,~]=UpdateSegmentMap_stal(blobIndIm,[]);
         blobIndIm2=blobIndIm1;
        %small median filter for final refinement
        Nsgms=Nclusters;
        tic;
        for iij=1:Nclusters
            maskseg=zeros(M,N);
            ppos=find(blobIndIm1==iij);
            if isempty(ppos)
                continue;
            end
            maskseg(ppos)=1;
            CC=bwconncomp(maskseg,4);
            L = labelmatrix(CC);
            blobIndIm2(ppos)=Nsgms+double(L(ppos))-1;
            Nsgms=Nsgms+max(double(L(ppos)));
        end
        blobIndIm1=medfilt2(blobIndIm2,[5 5]);
        blobIndIm=blobIndIm1;
      
         %this is not needed (duplicated) but since it was here for the paper, we keep it
        [Nclusters,neighbours,blobIndIm1,~]=UpdateSegmentMap_stal(blobIndIm,[]);
        %initialize final mask
        blobIndIm2=blobIndIm1;        
        Nsgms=Nclusters;
        
        for iij=1:Nclusters
            maskseg=zeros(M,N);
            ppos=find(blobIndIm1==iij);
            if isempty(ppos)
                continue;
            end
            maskseg(ppos)=1;
            CC=bwconncomp(maskseg,4);
            L = labelmatrix(CC);
            blobIndIm2(ppos)=Nsgms+double(L(ppos))-1;
            Nsgms=Nsgms+max(double(L(ppos)));
        end
        %final clean up
        blobIndIm1=medfilt2(blobIndIm2,[5 5]);
        [Nclusters,neighbours,blobIndIm,~]=UpdateSegmentMap_stal(blobIndIm,[]);
        segmentedmask = cleanupregions(reshape(blobIndIm,M,N), 100, 4);
        
        %final result
        [~,coloredsegment]=showsegm2(segmentedmask,im);
        [gbx,gby]=imgradientxy(segmentedmask);
        magg=sqrt(gbx.^2+gby.^2);
        bww=(magg>(mean(magg(:))));
        imgsegment = im;
        imgsegment(:,:,1) = 255*bwmorph(bww,'thin',Inf);

