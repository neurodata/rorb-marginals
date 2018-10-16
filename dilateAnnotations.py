
import numpy as np
import skimage.measure as measure
import scipy.ndimage as ndi
import tifffile as tiff
import pickle
import cv2
import os 
from slacker import Slacker

## Save data
with open('slack.ini', 'r') as fp:
    slack_token = fp.readline().strip()
    slack_user = fp.readline().strip()

slack = Slacker(slack_token)


tmp = sorted(os.listdir('annotations'))
filesA = sorted(['annotations/' + i for i in tmp])

anno = [tiff.imread(i).astype('uint16') for i in filesA[0:3]]

## make kernel for dilation
kernelD = np.ones((20 * 3, 20 * 3),np.uint8)

### dilate with OpenCV

for i in range(len(anno)):
    imD = cv2.morphologyEx(anno[i].astype('uint16'), cv2.MORPH_CLOSE, kernelD)
    tiff.imsave(file = 'dilatedAnnotations_k60/dilated_k60_' + str(i).zfill(3) + '.tif', data = imD)

slack.chat.post_message(slack_user, 'Synaptomes1 has finished dilation.')

