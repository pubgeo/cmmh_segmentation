function [segs,uids] = readSegs(present,iid)
% function [segs,uids] = readSegs(present,iid)
%
% Return a cell array of segmentations of an image 
% and the associated UIDs.
%
% David Martin <dmartin@eecs.berkeley.edu>
% January 2003

files = dir(fullfile(bsdsRoot,'human',present,'*'));
segs = {};
uids = {};
n = 0;
for i = 1:length(files),
  if ~files(i).isdir, continue; end
  if strcmp(files(i).name,'.'), continue; end
  if strcmp(files(i).name,'..'), continue; end
  uid = sscanf(files(i).name,'%d',1);
  file = segFilename(present,uid,iid);
  if length(dir(file))==0, continue; end
  n = n + 1;
  segs{n} = readSeg(file);
  uids{n} = uid;
end


[M,N]=size(im);
    %im=(-255)*(im-255)/255;
    gtemp=(myrenormalize(im));
    curve=bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;        
    
    %class1 = mean(gtemp(mask==1 & curve<=0));
    %class2 = mean(gtemp(mask==1 & curve>0));
    
    bbox = zeros(M,N);
    [r,c]=find(mask==1);
    if length(r)<50
        bbox(min(r(:)):max(r(:)),min(c(:)):max(c(:))) = 1;
        R1 = min(r(:));
        C1 = min(c(:));
        R2 = max(r(:));
        C2 = max(c(:));
    else
        padd=20;
        R1 = max(1,min(r(:))-padd);
        C1 = max(1,min(c(:))-padd);
        R2 = min(M,max(r(:))+padd);
        C2 = min(N,max(c(:))+padd);
        bbox(max(1,min(r(:))-padd):min(M,max(r(:))+padd),max(1,min(c(:))-padd):min(N,max(c(:))+padd)) = 1;
    end
    [r,c]=find(bbox==1);
        inimg=reshape(gtemp(bbox==1)',[max(r(:))-min(r(:))+1 max(c(:))-min(c(:))+1]);
        curvecrop=reshape(curve(bbox==1)',[max(r(:))-min(r(:))+1 max(c(:))-min(c(:))+1]);
%         inimg=zeros(M,N);
%         inimg(bbox==1)=gtemp(bbox==1);
%         curvecrop=curve;
        class1=0;class2=0;
        if R2-R1+1>1 && C2-C1+1>1            
            [u,diff]=mexTVsegmentationC(inimg',Niters,double(curvecrop)');
            u=reshape(u,[C2-C1+1 R2-R1+1])';            
            segm=zeros(M,N);
            segm(bbox==1) = u<=0;
        else
%             u = find(curve <= 1.2 & curve >= -1.2); 
%             segm=zeros(M,N);
%             segm = u>0;
            segm=mask;
        end