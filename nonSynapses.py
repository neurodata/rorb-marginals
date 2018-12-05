#! python3
import numpy as np
import math
import json
import skimage.measure as measure
import scipy.ndimage as ndi
import tifffile as tiff
import pickle
from slacker import Slacker
import os 
import sys

np.random.seed(1)

## multiplier to go from EM to AT/IF coordinates
mult = np.array([[1,0,0], [0, 1/32, 0], [0,0,1/32]])


## Set up slacker for status updates
with open('slack.ini', 'r') as fp:
    slack_token = fp.readline().strip()
    slack_user  = fp.readline().strip()

slack = Slacker(slack_token)


## Function to nudge annotation blobs within 
def nudgeB(B, zMax, yMax, xMax):
    ## calculate min and max for amount
    ## to nudge to stay withing bounds
    szmin = -1 * abs(np.min(B[:, 0]) - 0)
    szmax = abs(np.max(B[:, 0]) - zMax) - 1

    symin = -1 * abs(np.min(B[:, 1]) - 0)
    symax = abs(np.max(B[:, 1]) - yMax) - 1

    sxmin = -1 * abs(np.min(B[:, 2]) - 0)
    sxmax = abs(np.max(B[:, 2]) - xMax) - 1

    ## Sample random nudge in z, y, and x directions
    ## making sure to stay within bounds
    deltaz = np.random.randint(szmin, szmax, 1)[0]
    deltay = np.random.randint(symin, symax, 1)[0]
    deltax = np.random.randint(sxmin, sxmax, 1)[0]

    ## augment B to homogeneous coordinates
    Baug = np.hstack((B, np.ones(B.shape[0]).reshape(B.shape[0], 1))).astype(int)

    ## Create the affine transformation matrix
    A = np.hstack((np.vstack((np.identity(3), np.zeros(B.shape[1]))).astype(int), 
            np.array([deltaz, deltay, deltax, 2]).reshape(4,1).astype(int)))

    ## output the affine transformed data
    ## in non-homogeneous coordinates
    out = np.matmul(Baug, np.transpose(A))[:, 0:3].astype(int)
    return(out)



## Get the K = Ones(64x64) dilated annotations 
filesDA = sorted(['dilatedAnnotations_k60/' + i for i in os.listdir('dilatedAnnotations_k60')])

## remove problematic MacOS system file
for i in filesDA:
    if '.DS_Store' in i:
        filesDA.remove(i)

## build stack
dianno = [tiff.imread(i) for i in filesDA]

## Dilated annotations
Danno = np.stack(dianno)
Danno.shape

## get image bounds
zMax = Danno.shape[0]
yMax = Danno.shape[1]
xMax = Danno.shape[2]

## get non-zero locations 
atmp= np.transpose(np.asarray(np.where(Danno > 0)))

## instantiate dictionaries 
annoLoc = {}
nonSynapse = {}
synapse = {}

slack.chat.post_message(slack_user, 'Synaptomes1 has entered the loop to get annotations ...')

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
        print(key)
        sys.stdout.flush()
        ## nudge annotation i
        si = np.asarray(annoLoc[key])
        Bn = nudgeB(si, zMax, yMax, xMax)

        ## Check if the nudged annotation overlaps with an 
        ## existing annotation 
        runAgain = any(Danno[Bn[:, 0], Bn[:, 1], Bn[:, 2]] > 0)

        while(runAgain):
            Bn = nudgeB(si, zMax, yMax, xMax)
            runAgain = any(Danno[Bn[:, 0], Bn[:, 1], Bn[:, 2]] > 0)

        ai = np.unique(np.floor(np.matmul(Bn, mult)).astype('intp'), axis=0)

        ## Calculate sum over annotation pixels in IF space and
        ## normalize by number of pixels
        annoID = key
        numVox = ai.shape[0]
        synsum = np.sum(Syn[ai[:, 0], ai[:, 1], ai[:, 2]]) / numVox

        psdsum = np.sum(PSD[ai[:, 0], ai[:, 1], ai[:, 2]]) / numVox


        nonSynapse[str(key)] = {'annoID' : int(annoID), 
                                   'synapsin': float(synsum), 
                                   'psd95': float(psdsum), 
                                   'voxels': int(numVox)
                                   }



        with open('nonSynapse_16bit_k60.pickle', 'wb') as foo:
            pickle.dump(nonSynapse, foo, protocol=pickle.HIGHEST_PROTOCOL)



except Exception as e:
    print(e)
    slack.chat.post_message(slack_user, 'FAILED at loop building nonSynapses at key' +
            str(key) + '\n' + str(e))


## Save data
with open('nonSynapse_results_16bit_k64.json', 'w') as fp:
    json.dump(dict(nonSynapse), fp)



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

        with open('synapse_16bit_k60.pickle', 'wb') as foo:
            pickle.dump(synapse, foo, protocol=pickle.HIGHEST_PROTOCOL)



except Exception as e:
    print(e)
    slack.chat.post_message(slack_user, 'FAILED at loop.\n' + str(e))



## Save data
with open('synapse_results_16bit_k64.json', 'w') as fp:
    json.dump(dict(synapse), fp)



slack.chat.post_message(slack_user, 'Synaptomes1 has finished.')
