
import numpy as np
import json
import skimage.measure as measure
import scipy.ndimage as ndi
import tifffile as tiff
import pickle
#import cv2
import os 
from slacker import Slacker

## Set up slacker for status updates
with open('slack.ini', 'r') as fp:
    slack_token = fp.readline().strip()
    slack_user  = fp.readline().strip()

slack = Slacker(slack_token)


filesDA = sorted(['dilatedAnnotations_k64/' + i for i in os.listdir('dilatedAnnotations_k64')])
dianno = [tiff.imread(i) for i in filesDA]

Danno = np.stack(dianno)

## multiplier to go from EM to AT/IF coordinates
mult = np.array([[1,0,0], [0, 1/32, 0], [0,0,1/32]])

synapse = {}

atmp= np.transpose(np.asarray(np.where(Danno > 0)))

annoLoc = {}

slack.chat.post_message(slack_user, 'Synaptomes1 has entered the loop1 ...')

for i in range(atmp.shape[0]):
    annoLoc.setdefault(Danno[atmp[i,0], atmp[i,1], atmp[i,2]], []).append(list(atmp[i, :]))


filesSyn = sorted(['data16/synapsin_16bit/' + i for i in os.listdir('data16/synapsin_16bit')])
filesPSD = sorted(['data16/PSD95_16bit/' + i for i in os.listdir('data16/PSD95_16bit')])

syn = [tiff.imread(i) for i in filesSyn]
psd = [tiff.imread(i) for i in filesPSD]

## stacks are indexed Z, Y, X
Syn = np.stack(syn)
PSD = np.stack(psd)

slack.chat.post_message(slack_user, 'Synaptomes1 has entered the loop2 ...')

try:
    for key in annoLoc:
        ai = np.unique(np.floor(np.matmul(annoLoc[key], mult)).astype('intp'), axis=0)
    
        ## Calculate sum over annotation pixels in IF space and
        ## normalize by number of pixels
        annoID = key
        numVox = ai.shape[0]
        synsum = np.sum(Syn[ai[:, 0], ai[:, 1], ai[:, 2]]) / numVox
        psdsum = np.sum(PSD[ai[:, 0], ai[:, 1], ai[:, 2]]) / numVox
    
        synapse[str(key)] = {'annoID' : int(annoID), 
                           'synapsin': float(synsum), 
                           'psd95': float(psdsum), 
                           'voxels': int(numVox),
                           }

        with open('synapse_16bit_k64.pickle', 'wb') as foo:
            pickle.dump(synapse, foo, protocol=pickle.HIGHEST_PROTOCOL)

except Exception as e:
    print(e)
    slack.chat.post_message(slack_user, 'FAILED at loop.\n' + str(e))


## Save data
with open('results_16bit_k64.json', 'w') as fp:
    #json.dumps(synapse, fp)
    json.dump(dict(synapse), fp)



slack.chat.post_message(slack_user, 'Synaptomes1 has finished.')

