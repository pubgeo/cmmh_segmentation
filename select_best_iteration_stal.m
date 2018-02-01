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
%It takes all the segmentation hypothesis from all iterations and builds a cost function.
%for each pixel it selects the iteration hypothesis with minimum cost. Two solvers:
%solver=1 minimizes: COST+LABMDA*|label_i-label_j|
%solver~=1 minimizes: COST (=bilaccummetric)
%cost function definition:
%cost=a1*c1+a2*c2+a3*c3+a4*c4+a5*c5
%c1:L1 distance between pixel original input intensity and average
%intensity for segment pixel belongs to
%c2:label agreement among 4-neighbors
%c3:average intensity per segment difference between current iteration and
%previous and next iterations
%c4: std intensity per segment difference between current iteration and
%previous and next iterations
%c5: extra penalty for iterations that produce large segments


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

% Modified Version History:
%	10/06/2016  Marc Bosch   Initial version

function [c] = select_best_iteration_stal(Niter,im,accumean,accustd,blobIndImaccum,accumedge,edgeref,solver)
    %solver=1 minimizes: COST+LABMDA*|label_i-label_j|
    %solver~=1 minimizes: COST (=bilaccummetric)
    [M,N,ch] = size(im);        
    itercost(1,1,1:Niter-1) = 4*(1:Niter-1);
    %penalty for average intensity of segment diverge from original
    %input image
    costmean = sum(accumean,4);%costmean1+costmean2+costmean3;
    %cost for label disagreement (4-neighbors)
    temp=blobIndImaccum(:,:,2:Niter);
    temp(:,1:N-1,:)=blobIndImaccum(:,2:N,2:Niter);
    coneigl = blobIndImaccum(:,:,2:Niter)-temp(:,:,:);
    temp=blobIndImaccum(:,:,2:Niter);
    temp(:,2:N,:)=blobIndImaccum(:,1:N-1,2:Niter);
    coneigr = blobIndImaccum(:,:,2:Niter)-temp(:,:,:);
    temp=blobIndImaccum(:,:,2:Niter);
    temp(2:M,:,:)=blobIndImaccum(1:M-1,:,2:Niter);
    coneigb = blobIndImaccum(:,:,2:Niter)-temp(:,:,:);
    temp=blobIndImaccum(:,:,2:Niter);
    temp(1:M-1,:,:)=blobIndImaccum(2:M,:,2:Niter);
    coneigt = blobIndImaccum(:,:,2:Niter)-temp(:,:,:);
    costNeiglab=(coneigl~=0)+(coneigr~=0)+(coneigt~=0)+(coneigb~=0);
    %cost for average intensity of a segment changing from one iteration to another
    %(stability measure)
    %hh=reshape([-1 -1 -1 -1 8 -1 -1 -1 -1],[1 1 9]);
%     hh=reshape([-1 -1 4 -1 -1],[1 1 5]);
    hh=reshape([-1 2 -1],[1 1 3]);
    
    costdiffertaux=sum(abs(imfilter(accustd,hh,4)),4);
    costdiffert=costdiffertaux(:,:,1:end-1);
    costdiffertstdaux=abs(imfilter(accustd,hh,3));
    costdiffertstd=costdiffertstdaux(:,:,1:end-1);
   
    costedges = abs(accumedge(:,:,1:Niter-1)-repmat(edgeref,[1 1 Niter-1])).*repmat(edgeref,[1 1 Niter-1]);
    weightedge= 1*exp(abs(sum(edgeref(:)/(M*N))-(sum(sum(accumedge(:,:,2:Niter),1),2))/(M*N)));
    %     costdiffert(:,:,1) = 2*abs(accumean(:,:,2)-accumean(:,:,1));
%     costdiffert(:,:,2:Niter-1) = abs(accumean(:,:,1:end-2)-accumean(:,:,2:end-1)) + ...
%         abs(accumean(:,:,3:end)-accumean(:,:,2:end-1));
% %         costdiffert(:,:,Niter-1) = 2*abs(accumean(:,:,end-1)-accumean(:,:,end));
%     %cost for std intensity of a segment changing from one iteration to another
%     %(stability measure)
%     costdiffertstd(:,:,1) = 2*abs(accustd(:,:,2)-accustd(:,:,1));
%     costdiffertstd(:,:,2:Niter-1) = abs(accustd(:,:,1:end-2)-accustd(:,:,2:end-1)) + ...
%         abs(accustd(:,:,3:end)-accustd(:,:,2:end-1));
%         costdiffertstd(:,:,Niter-1) = 2*abs(accustd(:,:,end-1)-accustd(:,:,end));
% imshu8(edgeref)
% for iii=1:Niter-1
%     imshu8(accumedge(:,:,iii))
% end
% squeeze(weightedge)
% for iii=1:Niter-1
%     figure,subplot(1,2,1),imagesc(myrenormalize(weightedge(1,1,iii).*costdiffert(:,:,iii))),caxis([0 255]);
%     subplot(1,2,2),imagesc(blobIndImaccum(:,:,iii))
% end
%      itercost = repmat(itercost,[M N 1]);
     costedgeglob = repmat(weightedge.*sum(sum(costedges,1),2)/(255*M*N),[M N 1]);
% %         accumcostmetric = 0*costgrad+costmean+costNeiglab+30*costdiffert+30*costdiffertstd+itercost;
%     meas1=round(myrenormalize(costmean));
%     meas2=round(myrenormalize(costdiffert));
    edgeref=imdilate(edgeref,ones(5,5));
    for iii=1:Niter-1
%         w1te=0;
%         w2te=0;
%         vals=unique(meas1(:,:,iii));
%         vals2=unique(meas2(:,:,iii));
%         for ij=1:length(unique(meas1(:,:,iii)))
%             me=meas1(:,:,iii);
%             w1te=w1te+(length(find(me==vals(ij)))/(M*N))*log((M*N)/length(find(me==vals(ij))));
%         end
%         w1(1,1,iii)=(1/log(length(unique(meas1(:,:,iii)))))*w1te;
%         for ij=1:length(unique(meas2(:,:,iii)))
%             me=meas2(:,:,iii);
%             w2te=w2te+(length(find(me==vals2(ij)))/(M*N))*log((M*N)/length(find(me==vals2(ij))));
%         end
%         w2(1,1,iii)=(1/log(length(unique(meas2(:,:,iii)))))*w2te;
        oldlab=unique(blobIndImaccum(:,:,iii+1));
        blobIndIm1=blobIndImaccum(:,:,iii+1);
        costedgest=imdilate(accumedge(:,:,iii+1),ones(5,5));
        tempp=zeros(M,N);
        difp=imerode(abs(costedgest-edgeref),ones(2,2));
        for il=1:length(oldlab)            
            p=find(blobIndIm1==oldlab(il));           
            tempp(p)=sum(difp(p))/(255*max(1,length(find(edgeref(p)~=0))));
        end
        costedges(:,:,iii) = tempp;
%         figure,imagesc((costedges(:,:,iii))),caxis([min(costedges(:)) min(3,max(costedges(:)))]);
    end
     costedge=costedges;
%     w1t=w1;
%     w1=(ones(1,1,Niter-1)-w1)./((ones(1,1,Niter-1)-w1)+(ones(1,1,Niter-1)-w2));
%     w2=(ones(1,1,Niter-1)-w2)./((ones(1,1,Niter-1)-w1t)+(ones(1,1,Niter-1)-w2));
%     
    accumcostmetric = costedge;%+0.1*costdiffert+0.1*costmean;%1000*round(myrenormalize(costedge));%+repmat(w1,[M N 1]).*round(myrenormalize(costdiffert))+repmat(w2,[M N 1]).*round(myrenormalize(costdiffert));%+myrenormalize(costdiffert);%costmean+costNeiglab+1*costdiffert+1*costdiffertstd+itercost;

    ss1=exp(-((costdiffert).^2)/(2*std(costdiffert(:))));
    ss1=max(ss1(:))-ss1;
    ss=exp(-((costmean).^2)/(2*std(costmean(:))));
    ss=max(ss(:))-ss;
    ss2=exp(-((costedge).^2)/(2*std(costedge(:))));
    ss2=max(ss2(:))-ss2;
    meas1=round(myrenormalize(costdiffert));
    meas2=round(myrenormalize(costedge));
    for iii=1:Niter-1
        w1te=0;
        w2te=0;
        vals=unique(meas1(:,:,iii));
        vals2=unique(meas2(:,:,iii));
        for ij=1:length(unique(meas1(:,:,iii)))
            me=meas1(:,:,iii);
            w1te=w1te+(length(find(me==vals(ij)))/(M*N))*log((M*N)/length(find(me==vals(ij))));
        end
        w1(1,1,iii)=(1/log(length(unique(meas1(:,:,iii)))))*w1te;
        for ij=1:length(unique(meas2(:,:,iii)))
            me=meas2(:,:,iii);
            w2te=w2te+(length(find(me==vals2(ij)))/(M*N))*log((M*N)/length(find(me==vals2(ij))));
        end
        w2(1,1,iii)=(1/log(length(unique(meas2(:,:,iii)))))*w2te;
    end
    w1t=w1;
    w1=(ones(1,1,Niter-1)-w1)./((ones(1,1,Niter-1)-w1)+(ones(1,1,Niter-1)-w2));
    w2=(ones(1,1,Niter-1)-w2)./((ones(1,1,Niter-1)-w1t)+(ones(1,1,Niter-1)-w2));
%     squeeze(w1)
%     squeeze(w2)
    accumcostmetric=repmat(w2,[M N 1]).*ss1+repmat(w1,[M N 1]).*ss2;%+repmat(1,[M N Niter-1]).*costedgeglob;%ss1+ss(:,:,2:Niter)+5*ss2;
    
    bilaccummetric=accumcostmetric;
    %smooth cost metric among neighbors
    for iii=1:size(accumcostmetric,3)
        %[bilaccummetric(:,:,iii)]=bilatfilt(accumcostmetric(:,:,iii),imageToSegment(:,:,1),17,15,0);
        bilaccummetric(:,:,iii) = medfilt2(accumcostmetric(:,:,iii),[5 5]);
        bilaccummetric(:,:,iii) = FGS(accumcostmetric(:,:,iii), 0.01, 30*30, im,4,4);
       
         %   h=figure,imagesc(bilaccummetric(:,:,iii)),caxis([0 0.75]); axis off
         %   se=getframe(h);
         %   imwrite(uint8(se.cdata),['.\costtr_' num2str(iii) '.png']);
    end

    if solver==1 %minimizes:COST+LABMDA*|label_i-label_j|
        [c]= MRF_solve(bilaccummetric,0.00010,10);
    else
        [dummy,c]=min(bilaccummetric,[],3); 
    end
 
    
    
    