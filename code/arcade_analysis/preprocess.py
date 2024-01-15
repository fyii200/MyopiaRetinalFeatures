###########################################################
#                                                         #
#                  Author: Fabian Yii                     #
#               Email: fabian.yii@ed.ac.uk                #
#                                                         #
###########################################################
###########################################################

import os
import cv2 as cv
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate, distance_transform_edt
from skimage.morphology import area_opening
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.morphology import skeletonize
join = os.path.join


class arcade_detector(object):
    def __init__(self, name, mask_path):
        self.name = name
        self.mask = cv.imread(join(mask_path, name))
        
        # Convert to 2D if image has more than 2 dimensions
        if len(self.mask.shape) > 2:
            self.mask = cv.cvtColor(self.mask, cv.COLOR_BGR2GRAY)
    
    def view(self):
        plt.imshow(self.mask, cmap="gray")
        plt.axis("off")
        plt.show()
        
    def size(self):
        return np.shape(self.mask)
        
############################## Preprocessing functions ##############################            

    def dist_transform(self, min_trigger_area=1000, min_trigger_num=2, threshold_quantile=0.99):
        areas = self.compute_areas(self.mask)
        if sum(areas > min_trigger_area) >= min_trigger_num:
            dt_mask = distance_transform_edt(self.mask)
            thres = np.quantile(dt_mask, threshold_quantile)
            self.mask = np.uint8((dt_mask > thres)*255)

    def rectangular_opening(self, min_trigger_area=1000, min_trigger_num=2, min_angle=-40, max_angle=40, num_angles=20):
        areas = self.compute_areas(self.mask)
        if sum(areas > min_trigger_area) >= min_trigger_num:
            angles = np.linspace(min_angle, max_angle, num_angles)
            for angle in angles:
                kernel = np.ones((2,5),np.uint8) 
                kernel = rotate(kernel, angle)
                self.mask = cv.erode(self.mask, kernel, iterations = 1)
                self.mask = cv.dilate(self.mask, kernel, iterations = 1)
    
    def area_opening(self, min_trigger_area=1000, min_trigger_num=2, threshold_quantile=0.8, threshold_cap=300, connectivity=1):
        areas = self.compute_areas(self.mask)
        if sum(areas > min_trigger_area) >= min_trigger_num:
            area_threshold_interim = np.quantile(areas, threshold_quantile)
            area_threshold = area_threshold_interim if area_threshold_interim <= threshold_cap else threshold_cap
            self.mask = area_opening(self.mask, area_threshold, connectivity=connectivity)
            
    def detect_parabola(self, min_trial_radius=5, max_trial_radius=15, binary_quantile=0.985):
        hough_radii = np.arange(min_trial_radius, max_trial_radius+1, 1)
        hough_results = hough_circle(self.mask, hough_radii)
        # Select the most prominent circle
        _, _, _, radius = hough_circle_peaks(hough_results, hough_radii, total_num_peaks=1)
        
        self.mask = hough_circle(self.mask, radius[0])[0]
        binary_threshold = np.quantile(self.mask, binary_quantile)
        binary_mask = (self.mask > binary_threshold)*255
        self.mask = np.uint8(binary_mask)
        
    def rectangular_closing(self, min_angle=-70, max_angle=70, num_angles=15):
        angles = np.linspace(min_angle, max_angle, num_angles)
        for angle in angles:
            kernel = np.ones((20,2),np.uint8) 
            kernel = rotate(kernel, angle)
            self.mask = cv.dilate(self.mask, kernel, iterations = 1)
            self.mask = cv.erode(self.mask, kernel, iterations = 1)
            
    def skeleton(self):
        self.mask = np.uint8(skeletonize(self.mask)*255)    
        
    def crop_around_disc(self, disc_x, disc_y, width_after_cropped=450):
        self.eye = self.which_eye(self.name, disc_x, disc_y)
        height, width = self.mask.shape
        nasal_width = 50
        temporal_width = width_after_cropped - nasal_width
        try:
            if self.eye == "LE":
                disc_x = width*0.2 if (disc_x > width/2) or np.isnan(disc_x) else disc_x
                x_min = round(disc_x - nasal_width)
                x_max = round(disc_x + temporal_width)
            elif self.eye == "RE":
                disc_x = width*0.8 if (disc_x < width/2) or np.isnan(disc_x) else disc_x
                x_min = round(disc_x - temporal_width)
                x_max = round(disc_x + nasal_width)
            x_min = 0 if x_min < 0 else x_min
            x_max = width if x_max > width else x_max
        except:
            raise ValueError("disc_x, disc_y and/or width_after_cropped not specified")
        self.disc_x = disc_x
        self.disc_y = disc_y
        self.mask = self.mask[:, x_min:x_max]
        
    def pad(self, pad_width =10):
        self.mask = np.pad(self.mask, 
                           [(0, 0), (pad_width , pad_width )], 
                           mode='constant')
        
    def crop_to_fovea(self, disc_x, disc_y, fovea_x, fovea_y, disc_temporal_width=50):
        self.eye = self.which_eye(self.name, disc_x, disc_y)
        height, width = self.mask.shape
        try:
            if self.eye == "LE":
                disc_x = width*0.2 if (disc_x > width/2) or np.isnan(disc_x) else disc_x
                x_min = round(disc_x - disc_temporal_width)
                x_max = round(fovea_x) if ~np.isnan(fovea_x) else round(width/2)
            elif self.eye == "RE":
                disc_x = width*0.8 if (disc_x < width/2) or np.isnan(disc_x) else disc_x
                x_min = round(fovea_x) if ~np.isnan(fovea_x) else round(width/2)
                x_max = round(disc_x + disc_temporal_width)
            x_min = 0 if x_min < 0 else x_min
            x_max = width if x_max > width else x_max     
        except:
            raise ValueError("disc_x, disc_y, fovea_x and/or fovea_y not specified")  
        self.mask = self.mask[:, x_min:x_max]    
        
    def horizontal_cut(self, top):
        if top:
            top_seg = self.mask[0:round(self.disc_y),:]
            top_seg_flipped = cv.flip(top_seg, 0)
            return np.concatenate((top_seg,top_seg_flipped))
        else:
            bottom_seg = self.mask[round(self.disc_y):,:]
            bottom_seg_flipped = cv.flip(bottom_seg, 0)
            return np.concatenate((bottom_seg_flipped, bottom_seg))        
        
    def rotate(self, angle):
        image_center = tuple(np.array(self.mask.shape[1::-1]) / 2)
        rot_mat = cv.getRotationMatrix2D(image_center, angle, 1.0)
        self.mask = cv.warpAffine(self.mask, rot_mat, self.mask.shape[1::-1], flags=cv.INTER_LINEAR)
        
        
############################## Internal functions ##############################
    def compute_areas(self, mask):
        # apply connected component analysis
        output = cv.connectedComponentsWithStats(mask, cv.CV_32S)
        (numLabels, labels, stats, centroids) = output
        areas = stats[:, cv.CC_STAT_AREA][1:] # ignore the first element which corresponds to background pixels
        return areas
    
    def which_eye(self, name, disc_x, disc_y):
        try:
            if (name[11:13] == "15") & (disc_x <= (560/2)): 
                eye = "LE"
            elif (name[11:13] == "16") & (disc_x >= (560/2)):
                eye = "RE"
            elif disc_x <= (560/2):
                eye = "LE"
            elif disc_x >= (560/2):
                eye = "RE"
            elif np.isnan(disc_x) or np.isnan(disc_y):
                if name[11:13] == "15":
                    eye = "LE"
                elif name[11:13] == "16":
                    eye = "RE"  
        except:
            raise ValueError("Fundus image name must conform to UK Biobank file naming convention")
        return eye
        
        

