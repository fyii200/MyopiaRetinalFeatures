#!/usr/bin/env python
"""
This is the main executable script for deriving temporal arterial/venous concavity.

Author : Fabian Yii
Email  : fabian.yii@ed.ac.uk

2024
"""

import os
import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
from preprocess import arcade_detector
from arcade_model import parabola_ls, RANSAC
join = os.path.join

homeDir = "/Users/fabianyii/Desktop/masks/"
d1 = pd.read_csv(join(homeDir,"OD_fovea_parameters.csv"))
d2 = pd.read_csv(join(homeDir,"fovea_intensity_data.csv"))
d2 = d2.rename(columns={'fundus': 'name'})
d = pd.merge(d1, d2, left_on='name', right_on='name')[0:50000]
names = list(d.name)
            
##################### Set arguments #####################  
parser = argparse.ArgumentParser(description = "Vascular Arcade Analysis")
parser.add_argument("--parallelisation",
                    action = "store_true")
parser.add_argument("--crop_to_fovea", 
                    action = "store_true")
parser.add_argument("--image_name", 
                    type = str, 
                    default = None) 
parser.add_argument("--vessel_path", 
                    type = str, 
                    default = join(homeDir, "artery_vein", "artery", "full_width")) 
parser.add_argument("--disc_centroid_img_width", 
                    type = int, 
                    default = 560)
parser.add_argument("--disc_centroid_img_height", 
                    type = int, 
                    default = 560)
parser.add_argument("--width_after_cropped", 
                    type = int, 
                    default = 450)
parser.add_argument("--pad_width", 
                    type = int, 
                    default = 25)
parser.add_argument("--morph_open_min_trigger_area", 
                    type = int, 
                    default = 1000)
parser.add_argument("--morph_open_min_trigger_num", 
                    type = int, 
                    default = 2)
parser.add_argument("--rect_open_min_angle", 
                    type = int, 
                    default = -40)
parser.add_argument("--rect_open_max_angle", 
                    type = int, 
                    default = 40)
parser.add_argument("--rect_open_num_angles", 
                    type = int, 
                    default = 15)
parser.add_argument("--dist_trans_thres", 
                    type = float, 
                    default = 0.985)
parser.add_argument("--area_open_thres_quantile", 
                    type = float, 
                    default = 0.7)
parser.add_argument("--area_open_thres_cap", 
                    type = int, 
                    default = 300)
parser.add_argument("--hough_min_trial_radius", 
                    type = int, 
                    default = 5)
parser.add_argument("--hough_max_trial_radius", 
                    type = int, 
                    default = 15)
parser.add_argument("--hough_binary_quantile", 
                    type = float, 
                    default = 0.98)
parser.add_argument("--rect_close_min_angle", 
                    type = int, 
                    default = -70)
parser.add_argument("--rect_close_max_angle", 
                    type = int, 
                    default = 70)
parser.add_argument("--rect_close_num_angles", 
                    type = int, 
                    default = 15)
parser.add_argument("--segmented_lr_fit_intercept", 
                    action = "store_true")
parser.add_argument("--show_preprocessing", 
                    action = "store_true")
parser.add_argument("--show_ransac_parabola", 
                    action = "store_true")
parser.add_argument("--save_fig", 
                    action = "store_true")
parser.add_argument("--show_ransac_lin", 
                    action = "store_true")
parser.add_argument("--show_ls_parabola", 
                    action = "store_true")
parser.add_argument("--mark_vertex", 
                    action = "store_true")
parser.add_argument("--fit_verbose", 
                    action = "store_true")
args = parser.parse_args()

            
def build_pipeline(name):
    ################ Pre-processing ################
    ## Original full-sized mask (912x912)
    vessel = arcade_detector(name, args.vessel_path)
    original_mask = vessel.mask
    
    ## If mask is not empty
    if len(np.unique(original_mask)) > 1:
        # Crop the mask around the optic disc
        resize_ratio_width = vessel.size()[1] / args.disc_centroid_img_width
        resize_ratio_height = vessel.size()[0] / args.disc_centroid_img_height
        disc_x = list(d[d.name == name].disc_x)[0] * resize_ratio_width
        disc_y = list(d[d.name == name].disc_y)[0] * resize_ratio_height
        if args.crop_to_fovea:
            fovea_x = list(d[d.name==name].fovea_x)[0] * resize_ratio_width
            fovea_y = list(d[d.name==name].fovea_y)[0] * resize_ratio_height
            vessel.crop_to_fovea(disc_x, disc_y, fovea_x, fovea_y, 100)
        else:
            vessel.crop_around_disc(disc_x, disc_y, args.width_after_cropped)
        mask_cropped = vessel.mask
        
        # Pad mask
        vessel.pad(args.pad_width)
        
#         # Rotate
#         angle = list(d[d.name==name].vertical_angle)[0]
#         if abs(angle) <= 10:
#             vessel.rotate(angle)
        
        # Remove small vessels by distance transform
        vessel.dist_transform(args.morph_open_min_trigger_area,
                              args.morph_open_min_trigger_num,
                              args.dist_trans_thres)
        dist_trans_mask = vessel.mask
    
        # Further remove small vessels by morphological area opening
        vessel.area_opening(args.morph_open_min_trigger_area, 
                            args.morph_open_min_trigger_num, 
                            args.area_open_thres_quantile, 
                            args.area_open_thres_cap, 
                            connectivity=1)
        area_opened_mask = vessel.mask
    
        # Detect parabola via circle Hough transform
        vessel.detect_parabola(args.hough_min_trial_radius, 
                               args.hough_max_trial_radius, 
                               args.hough_binary_quantile)
        raw_parabola_mask = vessel.mask
    
        # Remove small vessels by morphological area opening
        vessel.area_opening(args.morph_open_min_trigger_area, 
                            args.morph_open_min_trigger_num, 
                            args.area_open_thres_quantile, 
                            args.area_open_thres_cap)
        parabola_mask_opened = vessel.mask
    
        # Image reconstruction via morphological closing using a rectangular kernel rotated between -70 & 70
        vessel.rectangular_closing(args.rect_close_min_angle, 
                                   args.rect_close_max_angle, 
                                   args.rect_close_num_angles)
        parabola_mask_closed = vessel.mask
    
        # Skeletonize the mask
        vessel.skeleton()
        skeletonised_mask = vessel.mask
        
        ## Display mask at each preprocessing step ##
        if args.show_preprocessing:              
            fig, p = plt.subplots(2, 4, figsize=(8, 8))
            p[0,0].imshow(mask_cropped, cmap='gray'); p[0,0].axis("off")
            p[0,1].imshow(dist_trans_mask, cmap='gray'); p[0,1].axis("off")
            p[0,2].imshow(area_opened_mask, cmap='gray'); p[0,2].axis("off")
            p[0,3].imshow(raw_parabola_mask, cmap='gray'); p[0,3].axis("off")
            p[1,0].imshow(parabola_mask_opened, cmap='gray'); p[1,0].axis("off")
            p[1,1].imshow(parabola_mask_closed, cmap='gray'); p[1,1].axis("off")
            p[1,2].imshow(skeletonised_mask, cmap='gray'); p[1,2].axis("off")
            plt.subplots_adjust(wspace=0, hspace=0.01)
            
        ################ Model fitting ################
        ## RANSAC parabola
        ransac = RANSAC(vessel.mask)
        ransac.fit_parabola()
        conc_rp, med_residual_rp, top_med_residual_rp, bottom_med_residual_rp, r2_rp = ransac.compute_metrics(model="parabola", verbose=args.fit_verbose)
        if args.show_ransac_parabola:
            ransac.display_parabola(args.mark_vertex)      
        
        ## Least square parabola
        if args.show_ls_parabola:
            parabola = parabola_ls(vessel.mask)
            parabola.fit()
            conc_lsp, med_residual_lsp, r2_lsp = parabola.compute_metrics(verbose=args.fit_verbose)
            parabola.display_fit(args.mark_vertex)
            
        ## Segmented linear model
        if args.show_ransac_lin:
            ransac.fit_segmented_lr(args.segmented_lr_fit_intercept)
            # Compute standard metrics
            med_residual_lin, r2_lin = ransac.compute_metrics(model="linear", verbose=args.fit_verbose)
            # compute parabolic index
            med_residual_ratio, r2_ratio = ransac.parabola_index() 
            ransac.display_segmented_lr() 
        
        # Save figure?
        if args.save_fig:
            save_dir = join(homeDir, "outputs", "ransac_parabola")
            os.makedirs(save_dir) if not os.path.exists(save_dir) else None  
            plt.savefig(join(save_dir, name), bbox_inches="tight") 
            plt.close()
        
        # Return metrics
        if args.show_ransac_lin:
            return name, conc_rp, med_residual_rp, top_med_residual_rp, bottom_med_residual_rp, r2_rp, med_residual_lin, r2_lin, med_residual_ratio, r2_ratio
        else:    
            return name, conc_rp, med_residual_rp, top_med_residual_rp, bottom_med_residual_rp, r2_rp
        
    else:
        print("{} is empty!".format(name))  
        
        
###### Run pipeline with parallel computing ######  
if __name__ == "__main__":
    if args.parallelisation:
        print("## start running pipeline with parallelisation ##")
        with multiprocessing.Pool() as p:
            results = list(tqdm(p.imap(build_pipeline, names), total=len(names)))
            p.close()
            df = pd.DataFrame(results)
            if args.show_ransac_lin:
                column_names = ["name", "conc_rp", "med_residual_rp", "top_med_residual_rp", "bottom_med_residual_rp", "r2_rp", "med_residual_lin", "r2_lin", "med_residual_ratio", "r2_ratio"]
            else:
                column_names = ["name", "conc_rp", "med_residual_rp", "top_med_residual_rp", "bottom_med_residual_rp", "r2_rp"]
            df = df.set_axis(column_names, axis=1)
            df.to_csv("result.csv", index=False)
    
    else:
        matplotlib.use('TkAgg')
        if args.show_ransac_lin:
            name, conc_rp, med_residual_rp, top_med_residual_rp, bottom_med_residual_rp, r2_rp, med_residual_lin, r2_lin, med_residual_ratio, r2_ratio = build_pipeline(args.image_name)
        else:    
            name, conc_rp, med_residual_rp, top_med_residual_rp, bottom_med_residual_rp, r2_rp = build_pipeline(args.image_name)
        plt.show()
        
        
        
        
        
        
