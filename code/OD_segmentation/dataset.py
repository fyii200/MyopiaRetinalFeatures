"""
Script defining the 'DatasetGenerator' class for getting and preprocessing fundus images on the UK Biobank Research Analysis Platform.

Author : Fabian Yii
Email  : fabian.yii@ed.ac.uk

2024
"""

import os
import torch
import cv2 as cv
import numpy as np
from torchvision.transforms import CenterCrop
from torch.utils.data import Dataset
join = os.path.join 


class DatasetGenerator(Dataset):
    def __init__(self, data_frame, croppedSize=(1400,1400)):
        """
        Args:
              data_frame (pandas df)       : Dataframe containing image names and their corresponding information.
              croppedSize (tuple of int)   : Desired output size of the centre crop (default to 1400 by 1400 pixels,
                                             as this removes the black/empty border around the image in UK Biobank.
        Return:
              image (3D array)             : Preprocessed colour fundus photograph. 
              image_name (str)             : Name of the retrieved fundus photograph.
        """

        self.image_names = list(data_frame.fundus)
        self.eyes        = list(data_frame.eye)
        self.crop        = CenterCrop(croppedSize)

    def __getitem__(self, index):
        """
        Out:
              image (4D tensor) : Preprocessed (cropped & normalised) image [B,C,H,W].
              image_name (str)  : Name of the image (along with the .png extension).
        """
            
        image_name = self.image_names[index]
        eye        = self.eyes[index]
        
        if eye == 'RE':
            data_dir = join('mnt', 'project', 'Bulk' 'Retinal Optical Coherence Tomography', 'Fundus (right)', str(image_name[0:2]) )
        elif eye == 'LE':
            data_dir = join('mnt', 'project' 'Bulk', 'Retinal Optical Coherence Tomography', 'Fundus (left)', str(image_name[0:2]) )
            
        image = cv.imread(join(data_dir, image_name))
        
        # Centre crop to "croppedSize"
        image        = self.crop(torch.tensor(image).moveaxis(2,0)).moveaxis(0,2).detach().numpy() 
        # Convert to RGB
        image        = cv.cvtColor(image, cv.COLOR_BGR2RGB).astype(np.float32) 
        # Resize to 560 by 560 pixels
        resizeWidth  = int(image.shape[1] * 40 / 100)  # scale down by 40%
        resizeHeight = int(image.shape[0] * 40 / 100) # scale down by 40%
        resizeDim    = (resizeWidth, resizeHeight)
        image        = cv.resize(image, resizeDim, interpolation = cv.INTER_AREA)
        # normalise to [0-1] range and reshape into [C,H,W]  
        image        = torch.tensor(image)
        image        = image.moveaxis(2,0) / torch.max(image) 
        
        return image, image_name

    def __len__(self):
        return len(self.image_names)
