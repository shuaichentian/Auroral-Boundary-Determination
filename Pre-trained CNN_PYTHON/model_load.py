from keras.models import load_model ,Model
import numpy as np
from PIL import Image
from PIL import ImageEnhance
import cv2
import os

'''对图片进行测试'''
model=load_model('pre-trained_CNN_model.h5')

'''将图片转换为所需要的格式'''
filepath='Image_test'


for filename in os.listdir(filepath):
    imgfile=filepath+'/'+filename
    img = Image.open(imgfile)
    en_con = ImageEnhance.Contrast(img)
    contrast = 3
    img = en_con.enhance(contrast)
    img=img.convert('RGB')
    img = cv2.cvtColor(np.asarray(img),cv2.COLOR_RGB2GRAY)


    '''get the confidence map'''
    img=cv2.copyMakeBorder(img,5,5,5,5,cv2.BORDER_CONSTANT)
    img = np.array(img)
    (c, b) = np.shape(img)
    Data = []

    for i in range(5, c - 5):
        for j in range(5, b - 5):
            data = img[i - 5:i + 6, j - 5:j + 6]
            data = np.array([data])
            Data.append(data)
    Data = np.array(Data)
    Data = Data.reshape(-1, 11, 11, 1)

    Data = Data / 255.0
    file_name=filename.strip('.bmp')



    '''用来输出每个点的概率值'''
    intermediate_layer_model = Model(input=model.input, output=model.get_layer('logits').output)
    intermediate_output = intermediate_layer_model.predict(Data, batch_size=500, verbose=0)
    logits_img = [x[1] for x in intermediate_output]
    logits_img = (np.array(logits_img)) * 255
    logits_img=logits_img.astype(int)
    logits_img = logits_img.reshape(228, 200)
    cv2.imwrite('./Confidence map/'+file_name+'_confidence map.bmp',logits_img)









