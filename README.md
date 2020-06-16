## CMMH Segmentation 
Cost Minimization Multiple Hypotheses Segmentation 
 
This repository contains all of the necessary code, imagery, and information to reproduce JHUAPL's results published at
the British Machine Vision Conference (BMVC) in 2017.
 
M.Bosch, C. Gifford, A.Dress, C. Lau, J.Skibbo, G.Christie 
"Improved image segmentation via cost minimization of multiple hypotheses", BMVC 17 
[[pdf]](./0939.pdf)
[[mp4]](https://www.dropbox.com/s/afri5vgbhvovsn9/0939.mp4?dl=1). 
If code is used please cite our paper.

This version only supports FH kernel as described in the paper.

Results from BVMC use data in [BSDS300](./BSDS300) sourced from:
D. Martin, C. Fowlkes, D. Tal, and J. Malik. A database of human segmented natural
images and its application to evaluating segmentation algorithms and measuring eco-
logical statistics. In Proceedings of the IEEE International Conference on Computer
Vision (ICCV), pages 416–423, 2001.

### Requirements
The following Matlab resources (and their dependencies) are required:
* Matlab 2014 or later (tested 2017a)
* Toolboxes
  * image_toolbox
  * matlab

The software only comes pre-compiled for windows (most versions). For linux/mac mex compilation of the following files is needed before running code.
* [./3rdparty/mexFGS.cpp](./3rdparty/mexFGS.cpp)
* [mexFelzenSegmentIndexFeatures.cpp](./mexFelzenSegmentIndexFeatures.cpp)

### CMMH Segmentation Usage
    segm_bmvc
This command will run the experiments presented in the above reference as presented at BMVC using data in the [BSDS300](./BSDS300) directory.
    
    segm_cmmh(image, segmentation, hypotheses, bpm);
This command will run on the image provided by the user.
* **Inputs**
  * _image_- An matrix representation of an image. IE double(imread('my_image.png'))
  * _segmentation_- Tunable parameter that controls degree of over/under segmentation of each hypothesis in FH kernel
  * _hypotheses_- Number of hypotheses, ideally we want a number that covers the volume space and variation in hypotheses
  * _bpm_- Flag for belief propagation minimizer (slower, a bit better). 1 to activate, 0 to disable
* **Outputs** (mxn matrices of doubles)
  * _SegmentedMask_: each pixel is assigned the id of the segment it belongs to. It will render as greyscale
  * _ColoredSegment_: each pixel is assigned the color of the average pixel within its respective segment
  * _ImageSegment_: Input image reproduced with red boundaries denoting segment outlines 

## Geospatial Aerial Examples
The [example](./example) folder contains examples of aerial imagery from 
"S. Razakarivony and F. Jurie. Vehicle detection in aerial imagery. 
Journal on Visual Communication and Image Representation, 34(C):187–203, Jan. 2016. 
processed with the following parameters:
* segmentation = 150
* hypotheses = 30
* bpm = 0

## Troubleshooting
Some Matlab versions trigger an error claiming 'myrenormalize' cannot be found.
Either type addpath(genpath('./')) in the Matlab console or just copy 'myrenormalize.m' in the '3rdparty' folder
