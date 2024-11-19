"""
Author : Fabian Yii
Email  : fabian.yii@ed.ac.uk

2024
"""

import numpy as np
import torch
import torch.nn as nn
    
def show_mask(mask, ax, random_color=False):
    """
    Args:
          mask (2D ndarray)   : Predicted optic disc segmentation mask (binary).
          ax (str)            : plt ax object showing the original image.
          random_color (bool) : Colour of the mask overlay (random if True).
    """
    if random_color:
        color = np.concatenate([np.random.random(3), np.array([0.6])], axis=0)
    else:
        color = np.array([30/255, 144/255, 255/255, 0.6])
    h, w = mask.shape[-2:]
    mask_image = mask.reshape(h, w, 1) * color.reshape(1, 1, -1)
    ax.imshow(mask_image)   



    


