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

##Demo result
The confidence map obtained from Pre-trained CNN.
![image]
