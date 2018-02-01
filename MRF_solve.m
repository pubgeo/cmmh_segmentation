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
%
% solves the following cost metric by belief propagation and min-sum: COST+LABMDA*|label_i-label_j| for
% smoother results
%
% datacost:              cost to minimize: MxNxNlabels
% lambda :              Smoothness weighting
% Niter:                Number of iterations (typically it converges with
% 10 or less)
%
% output:                minimized field


% Get initial segmentation, boxes, and neighbouring blobs
% Modified Version History:
%	09/19/2016  Marc Bosch   Initial modified version

function [output]= MRF_solve(datacost,lambda,Niter)
    
    
    Nlabels=size(datacost,3);    
    field.Nlabels=Nlabels;
    field.h = size(datacost,1);
    field.w = size(datacost,2);
    field.lambda = lambda;
    field.message = zeros(field.h,field.w,Nlabels,4);
    field.minlabel=zeros(field.h,field.w);
    %smoothness
    [x,y]=meshgrid(1:Nlabels,1:Nlabels);
    field.pen = lambda*min(abs(y-x),(Nlabels/2)*ones(Nlabels,Nlabels));
    prevenergy = 1000000*size(datacost,1)*size(datacost,2)*size(datacost,3);
    prevlabel = ones(size(datacost,1),size(datacost,2));
    for bp=1:Niter
        
     %   tic;               
        field=BeliefProp(field,datacost,0);%0 passing to the right
        field=BeliefProp(field,datacost,1);%1 passing to the left
        field=BeliefProp(field,datacost,2);%2 passing to up
        field=BeliefProp(field,datacost,3);%3 passing to down              
        penalty = datacost + field.message(:,:,:,1) + ...
            field.message(:,:,:,2) + ...
            field.message(:,:,:,3) + ...
            field.message(:,:,:,4);
        
        [val,pos]=sort(penalty,3);
        minlabel=pos(:,:,1);
        field.minlabel = minlabel;    
        labelleft=minlabel;labelright=minlabel;
        labelup=minlabel;labeldown=minlabel;
        labeldown(2:field.h,:)=minlabel(1:field.h-1,:);
        labelup(1:field.h-1,:)=minlabel(2:field.h,:);
        labelleft(:,1:field.w-1)=minlabel(:,2:field.w);
        labelright(:,2:field.w)=minlabel(:,1:field.w-1);
        dc = zeros(field.h,field.w);
        for iii=1:Nlabels
            blobtemp=datacost(:,:,iii);
            [rr,cc]=find(minlabel==iii);
            p=sub2ind(size(blobtemp), rr, cc);
            dc(p)=blobtemp(p); 
        end
         energy = dc + ...
             lambda*abs(field.minlabel-labelleft)+...
             lambda*abs(field.minlabel-labelright)+...
             lambda*abs(field.minlabel-labelup)+...
             lambda*abs(field.minlabel-labeldown);
         fprintf('energy %f %f\n',sum(energy(:)), sum(dc(:)));

       % toc;
        %convergence
        if (sum(energy(:))-prevenergy)/prevenergy>(-1)*(0.001)
            break;
        else
            prevenergy=sum(energy(:));
            prevlabel=field.minlabel;
        end
            
    end
    output=prevlabel;
end



function [penalty]=smoothnesspenalty(labela,labelb,lambda)
    penalty = lambda*min(abs(labela-labelb),1);    
end

function [field] = BeliefProp(fieldin,datacost,flagdir)
    
    coef=[1 1 1 1];
    if flagdir==0%right
        coef(1) = 0;
        pen=repmat(reshape(fieldin.pen,[1 fieldin.Nlabels fieldin.Nlabels]),[size(fieldin.message,1) 1 1]);
        for j=1:fieldin.w-1
            fieldin.message(:,j+1,:,2)=Messageto(fieldin,-1,j,repmat(squeeze(datacost(:,j,:)),[1 1 fieldin.Nlabels]),coef, ...
                pen);
%             squeeze(fieldin.message(:,j+1,:,2))
        end            
    elseif flagdir==1%left
        coef(2) = 0; 
        pen=repmat(reshape(fieldin.pen,[1 fieldin.Nlabels fieldin.Nlabels]),[size(fieldin.message,1) 1 1]);
        for j=fieldin.w:-1:2
            fieldin.message(:,j-1,:,1)=Messageto(fieldin,-1,j,repmat(squeeze(datacost(:,j,:)),[1 1 fieldin.Nlabels]),coef, ...
                pen);
        end        
    elseif flagdir==2%up
        coef(3) = 0;         
        pen=repmat(reshape(fieldin.pen,[1 fieldin.Nlabels fieldin.Nlabels]),[size(fieldin.message,2) 1 1]);
        for i=fieldin.h:-1:2
            fieldin.message(i-1,:,:,4)=Messageto(fieldin,i,-1,repmat(squeeze(datacost(i,:,:)),[1 1 fieldin.Nlabels]),coef, ...
                pen);
        end
    elseif flagdir==3%down
        coef(4) = 0;
        pen=repmat(reshape(fieldin.pen,[1 fieldin.Nlabels fieldin.Nlabels]),[size(fieldin.message,2) 1 1]);
        for i=1:fieldin.h-1
            fieldin.message(i+1,:,:,3)=Messageto(fieldin,i,-1,repmat(squeeze(datacost(i,:,:)),[1 1 fieldin.Nlabels]),coef, ...
                pen);
        end
    end
    field = fieldin;
end

function [update] = Messageto(fieldin,i,j,datacost,coef,pen)
    if coef(1)==0 || coef(2)==0 %accross cols <->
        cost = pen + datacost +...
            coef(1)*repmat(squeeze(fieldin.message(:,j,:,1)),[1 1 fieldin.Nlabels])+...
            coef(2)*repmat(squeeze(fieldin.message(:,j,:,2)),[1 1 fieldin.Nlabels])+...
            coef(3)*repmat(squeeze(fieldin.message(:,j,:,3)),[1 1 fieldin.Nlabels])+...
            coef(4)*repmat(squeeze(fieldin.message(:,j,:,4)),[1 1 fieldin.Nlabels]);
    else
        cost = pen + datacost +...
            coef(1)*repmat(squeeze(fieldin.message(i,:,:,1)),[1 1 fieldin.Nlabels])+...
            coef(2)*repmat(squeeze(fieldin.message(i,:,:,2)),[1 1 fieldin.Nlabels])+...
            coef(3)*repmat(squeeze(fieldin.message(i,:,:,3)),[1 1 fieldin.Nlabels])+...
            coef(4)*repmat(squeeze(fieldin.message(i,:,:,4)),[1 1 fieldin.Nlabels]);
    end
    
    update=reshape(min(cost,[],2),[size(cost,1) 1 size(cost,2)]); 
    
    if coef(3)==0 || coef(4)==0
        update=reshape(min(cost,[],2),[1 size(cost,1) size(cost,2)]);
    end

    
end