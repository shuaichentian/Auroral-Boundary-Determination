# Auroral-Boundary-Determination
This project is used to extract the auroral boundaries from UVI images captured by Polar spacecraft. The corresponding work has been submitted to JGR-Space Physics, and the paper name is "Automatic Auroral Boundary Determination Algorithm with Deep Feature and Dual Level Set." 

If you intend to use the code below in your work, please cite the paper as "..." (The paper has not yet been accepted currently).

## Software we used to run the code

1. Python 3.6+
>Packages for Python:
>>OpenCV 3.4.4;   
>>tensorflow 1.12.2;    
>>keras 2.2.4;
2. Matlab2018

## Quick start
>The folder `Pre-trained CNN_PYTHON` contains the code for extracting the confidence map (see the paper for more details).  
>>The subfolder `Image_test` is used to store the UVI images.  
>>The subfolder `Confidence map` is used to store the extracted confidence map corresponding to the UVI images in `Image_test`.   
>>The `pre-trained_CNN_model.h5` contains the pre-trained model and weight parameters.  
>>Run the following command to obtain the confidence map for corresponding UVI images.     

```
python model_load.py
```
>The folder `Extract auroral boundary_MATLAB` contains the code for segmenting the auroral oval from UVI images.  
>>The subfolder `Image_test` contains the UVI images.  
>>The subfolder `Confidence map` contains the confidence map you obtained from previous step (LLSRHT result).  
>>The subfolder `ini` contains the initialization image for inner and outer zero level set curves.  
>>The subfolder `segmentation_result` is used to store the extracted auroral oval image (`.bmp` format) and the boundaries position (`.mat` format).  
>>Run the following command to extract the auroral boundaries from UVI images.  
```
matlab DemoUVI_tian1.m
```

## Demo result
The confidence map obtained from Pre-trained CNN.<br>

<div align="center"><img width="150" height="170" src="https://github.com/shuaichentian/Auroral-Boundary-Determination/blob/master/Pre-trained%20CNN_PYTHON/Confidence%20map/1997010_034012_confidence%20map.bmp"/></div>

The boundary determination process.
<div align="center"><img width="330" height="250" src="https://github.com/shuaichentian/Auroral-Boundary-Determination/blob/master/Extract%20auroral%20boundary_MATLAB/ectract_boundary.gif"/></div>

