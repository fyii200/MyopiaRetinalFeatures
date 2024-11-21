"""
Script defining the 'trainedModel' class, which loads a trained DeepLabV3 Mobilenet model.

Author : Fabian Yii
Email  : fabian.yii@ed.ac.uk

2024
"""

import os
import torch
import torch.nn as nn
import cv2 as cv
import numpy as np
import matplotlib.pyplot as plt
from OD.utils import show_mask
join = os.path.join 


class trainedModel:
    def __init__(self, savedMaskDir, checkpoint="bestODmodel.pth"):
        """
        Args:
              savedMaskDir (str)             : Path to the folder in which the predicted masks are to be saved.
              checkpoint (str)               : Name of the file containing the trained weights (.pth).
        
        Return:    
              predictedFinalMask (2D array)  : Predicted optic disc mask [0, 255].
        """
        
        self.savedMaskDir = savedMaskDir
        os.makedirs(self.savedMaskDir, exist_ok = True) 
        self.checkpoint   = checkpoint
        self.device       = torch.device("cuda" if torch.cuda.is_available() else "mps")

    def load(self, dir_path):
        """
        Args:
              dir_path (str): Path to the trained weights.
        """
        
        # Load DeepLabV3 Mobilenet and its trained weights.
        self.model = torch.hub.load('pytorch/vision:v0.10.0', 
                                    'deeplabv3_mobilenet_v3_large')
        self.model.classifier[4]     = nn.LazyConv2d(1, 1)
        self.model.aux_classifier[4] = nn.LazyConv2d(1, 1)
        checkpoint_path              = join(dir_path, self.checkpoint)
        self.model.load_state_dict(torch.load(checkpoint_path, map_location = self.device))
        self.model.to(self.device)   
        self.model.eval()    
    
    def segment(self, image, image_name, blurKernelSize = (21,21), binaryThreshold = 0.95):
        """
        Args:
              image (4D ndarray)              : Input image of shape [B,C,H,W].
              image_name (str)                : Name of the image.
              blurKernelSize (int)            : Size of Gaussian blur filter used to smooth the boundary.
              binaryThreshold (float)         : Threshold above which a pixel is assigned 1,
                                                while below which = 0.
        Out:
              predictedFinalMask (2D ndarray) : Binary predicted optic disc mask (0. or 1.).
        """
        
        with torch.no_grad():
            outputs             = self.model(image.to(self.device))['out'] # (B, 1, H, W)
            output              = outputs[0,0,:,:] # (H, W)
            output              = cv.GaussianBlur(output.detach().cpu().numpy(), blurKernelSize, 0)
            output              = torch.sigmoid(torch.tensor(output)).detach().cpu().numpy()
            predictedBinaryMask = np.where(output >= binaryThreshold, 1, 0)
        
        ########################################### Quality control ###########################################
        #  If there's more than 1 predicted mask, select the mask that most probably corresponds to the optic #
        #  disc. This is done by selecting the most circular shape. If there's still more than one predicted  #
        #                     mask after this step, select the mask with the largest area.                    #
        ####################################################################################################### 
        cnts, hierarchy = cv.findContours(np.uint8(predictedBinaryMask), cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
        if len(cnts)>1:
            contour_list = []
            for contour in cnts:
                approx = cv.approxPolyDP(contour, 0.01*cv.arcLength(contour, True), True)
                area   = cv.contourArea(contour)
                if ((len(approx)>8) & (area > 30)): # only add circular contour to the list
                    contour_list.append(contour)
            if len(contour_list) != 0:
                largest_c = max(contour_list, key = cv.contourArea)
            else:
                largest_c = max(cnts, key = cv.contourArea)
            predictedFinalMask = cv.drawContours(np.zeros(image[0,0,:,:].shape), [largest_c], -1, 1, -1)
        else:
            predictedFinalMask = cv.drawContours(np.zeros(image[0,0,:,:].shape), cnts, -1, 1, -1)
   
        # Save predicted mask
        cv.imwrite(join(self.savedMaskDir, image_name), predictedFinalMask*255)
        
        # # Save original image with the mask overlay.
        # image = image[0,:,:,:].moveaxis(0,2).detach().cpu().numpy()
        # plt.imshow(image); plt.axis("off")
        # show_mask(predictedFinalMask*0.6, plt.gca())
        # plt.savefig(join(self.savedMaskDir, "overlay_" + image_name) )
        # plt.close()

        return predictedFinalMask
