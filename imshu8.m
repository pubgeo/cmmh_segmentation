function handle=imshu8(im)
    im=double(im);
    handle=figure,imshow(uint8(255*im/max(im(:)))),axis on