# Auroral-Boundary-Determination
This project is used to extract the auroral boundaries from UVI images captured by Polar spacecraft.

## Requirement

1. Python 3.6+
>Packages for Python:
>>OpenCV 3.4.4;   
>>tensorflow 1.12.2;    
>>keras 2.2.4);
2. Matlab

## Quick start
Run the following command to obtain the confidence map for corresponding UVI images. (The `pre-trained_CNN_model.h5` contains the pre-trained model and weight parameters )
```
python model_load.py
```

Run the following command to extract the auroral boundaries from UVI images.
```
matlab DemoUVI_tian1.m
```

## Demo result
The confidence map obtained from Pre-trained CNN.<br>

<div align="center"><img width="150" height="170" src="https://github.com/shuaichentian/Auroral-Boundary-Determination/blob/master/Pre-trained%20CNN_PYTHON/Confidence%20map/1997010_034012_confidence%20map.bmp"/></div>

The boundary determination process.
<div align="center"><img width="330" height="250" src="https://github.com/shuaichentian/Auroral-Boundary-Determination/blob/master/Extract%20auroral%20boundary_MATLAB/ectract_boundary.gif"/></div>

## Explanation
The code for "Aurora boundary extraction" will be uploaded after the corresponding paper is published.
