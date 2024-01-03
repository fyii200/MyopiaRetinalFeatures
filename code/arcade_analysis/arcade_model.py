import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.preprocessing import PolynomialFeatures
from preprocess import arcade_detector
from sklearn.linear_model import RANSACRegressor, LinearRegression
import piecewise_regression

class parabola_ls(arcade_detector):
        def __init__(self, mask):
            self.mask = mask
            coords = vessel_pixel_coordinates(self.mask)
            self.x = coords[:,0]
            self.y = coords[:,1]
            
        def quadratic_vertex_equation(self, y, concavity, vertex_y, vertex_x):
            # Quadratic function in vertex form
            return concavity*(y - vertex_y)**2 + vertex_x
        
        def fit(self):
            self.coefficients, _ = curve_fit(self.quadratic_vertex_equation, self.y, self.x, maxfev=1000)
            self.concavity, self.vertex_y, self.vertex_x = self.coefficients
            
        def predict(self, y_input):
            return self.quadratic_vertex_equation(y_input, *self.coefficients)
    
        def display_fit(self, mark_vertex=True):
            # define a sequence of inputs between the smallest and largest known inputs
            y_input = np.arange(min(self.y), max(self.y), 1)
            # calculate the output for the range
            x_predicted = self.predict(y_input)
            # Plot
            plt.imshow(self.mask, cmap='gray')
            plt.plot(self.x, self.y, '.', markersize=1.5, color="red")              
            show_boolean = (x_predicted <= self.mask.shape[1]) & (x_predicted > 0)
            plt.plot(x_predicted[show_boolean], y_input[show_boolean], '--', color="tomato")
            if mark_vertex:
                plt.scatter(self.vertex_x, self.vertex_y, marker="x", color="tomato", s=70)
            plt.axis('off')
            
        def compute_metrics(self, decimal_place=4, verbose=True):
            residuals = self.x - self.predict(self.y)
            median_residual = np.median(abs(residuals))
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.x-np.mean(self.x))**2)
            r_squared = 1 - (ss_res / ss_tot)
            r_squared = np.round(r_squared, decimal_place)
            if verbose:
                print("#### least square parabola ####", )
                print("Concavity index: {:.4f}".format(abs(self.concavity)))
                print("Median residual: {:.4f}".format(median_residual))
                print("R2: {:.4f}".format(r_squared))
            return abs(self.concavity), median_residual, r_squared
            
            
class RANSAC(arcade_detector):
    def __init__(self, mask):
        self.mask = mask
        coords = vessel_pixel_coordinates(self.mask)
        self.x = coords[:,0]
        self.y = coords[:,1]
        self.x_reshaped = self.x.reshape((-1, 1))
        self.y_reshaped = self.y.reshape((-1, 1))
    
###################### Parabola ########################
    def fit_parabola(self):
        self.quadratic = PolynomialFeatures(degree=2)
        self.y_quadratic_trans = self.quadratic.fit_transform(self.y_reshaped)
        self.parabola = RANSACRegressor(LinearRegression(fit_intercept=False), min_samples=15)
        self.parabola.fit(self.y_quadratic_trans, self.x_reshaped)
        # get coefficients
        self.c = self.parabola.estimator_.coef_[0, 0] 
        self.b = self.parabola.estimator_.coef_[0, 1] 
        self.concavity = self.parabola.estimator_.coef_[0, 2]  
    
    def display_parabola(self, mark_vertex=True):
        y_input = np.arange(0, self.mask.shape[0])
        y_input = y_input.reshape((-1,1))
        y_input_trans = self.quadratic.fit_transform(y_input)
        self.x_predicted = self.parabola.predict(y_input_trans)
        # Plot
        plt.imshow(self.mask, cmap='gray')
        show_line_boolean = (self.x_predicted <= self.mask.shape[1]) & (self.x_predicted > 0)
        plt.plot(self.x[self.parabola.inlier_mask_], 
                 self.y[self.parabola.inlier_mask_], 
                 '.', markersize=1.5, color='red')
        plt.plot(self.x[~self.parabola.inlier_mask_], 
                 self.y[~self.parabola.inlier_mask_], 
                 '.', markersize=1.5, color='green') 
        plt.plot(self.x_predicted[show_line_boolean], 
                 y_input.reshape((-1,1))[show_line_boolean], 
                 '--', color='orange') 
        if mark_vertex:
            vertex_x, vertex_y = self.get_parabola_vertex()
            plt.scatter(vertex_x, vertex_y, marker="x", color="orange", s=70)
        plt.axis('off')  
        
    def get_parabola_vertex(self):
        # vertex point of the predicted parabola
        vertex_y = -self.b / (2*self.concavity)
        vertex_x = self.concavity*vertex_y**2 + self.b*vertex_y + self.c
        return vertex_x, vertex_y   

    
############# Segmented linear regression ###############
    def get_inliers(self):
        try:
            x_inliers = self.x[self.parabola.inlier_mask_] 
            y_inliers = self.y[self.parabola.inlier_mask_] 
            return x_inliers, y_inliers
        except:
            raise NameError("Optic disc coordinates not found; parabola must first be fitted before calling this function")
    
    def segment_mask(self):
        x_inliers, y_inliers = self.get_inliers()
        self.vertex_x, self.vertex_y = self.get_parabola_vertex()
        ## First segment of the linear regression (from top to disc centre)
        first_segment_boolean = y_inliers > self.vertex_y
        self.x_inliers_1 = x_inliers[first_segment_boolean].reshape((-1, 1))
        self.y_inliers_1 = y_inliers[first_segment_boolean].reshape((-1, 1))
        ## Second segment of the linear regression (from disc centre to bottom)
        sec_segment_boolean = y_inliers <= self.vertex_y
        self.x_inliers_2 = x_inliers[sec_segment_boolean].reshape((-1, 1))
        self.y_inliers_2 = y_inliers[sec_segment_boolean].reshape((-1, 1))
    
    def fit_segmented_lr(self, fit_intercept):
        self.segment_mask()
        self.centre_factor_x = 0 if fit_intercept else self.vertex_x
        self.centre_factor_y = 0 if fit_intercept else self.vertex_y
        ## Fit the first segment
        if len(self.y_inliers_1) > 0:
            self.linreg1 = LinearRegression(fit_intercept=fit_intercept)
            self.linreg1.fit(self.y_inliers_1 - self.centre_factor_y, self.x_inliers_1 - self.centre_factor_x)
        ## Fit the second segment    
        if len(self.y_inliers_2) > 0:
            self.linreg2 = LinearRegression(fit_intercept=fit_intercept)
            self.linreg2.fit(self.y_inliers_2 - self.centre_factor_y, self.x_inliers_2 - self.centre_factor_x)
        
    def display_segmented_lr(self):
        centre_point = self.vertex_y - self.centre_factor_y
        ## Plot
        plt.imshow(self.mask, cmap='gray')
        # First segment
        if len(self.y_inliers_1) > 0:
            top_point = np.max(self.y_inliers_1 - self.centre_factor_y) 
            y_input_1 = np.arange(centre_point, top_point) 
            y_input_1 = y_input_1.reshape((-1, 1))
            x_predicted_1 = self.linreg1.predict(y_input_1)
            show_line_boolean = abs(x_predicted_1 + self.centre_factor_x) <= self.mask.shape[1]
            plt.plot(x_predicted_1[show_line_boolean] + self.centre_factor_x, 
                     y_input_1[show_line_boolean] + self.centre_factor_y, 
                     "--", color="gray")            
        # Second segment
        if len(self.y_inliers_2) > 0:
            bottom_point = np.min(self.y_inliers_2 - self.centre_factor_y) 
            y_input_2 = np.arange(bottom_point, centre_point) 
            y_input_2 = y_input_2.reshape((-1, 1))
            x_predicted_2 = self.linreg2.predict(y_input_2)             
            show_line_boolean = abs(x_predicted_2 + self.centre_factor_x) <= self.mask.shape[1]
            plt.plot(x_predicted_2[show_line_boolean] + self.centre_factor_x, 
                     y_input_2[show_line_boolean] + self.centre_factor_y, 
                     "--", color="gray") 
        # Plot inliers & outliers
        plt.plot(self.x[self.parabola.inlier_mask_], self.y[self.parabola.inlier_mask_], '.', markersize=1.5, color='red')   
        plt.plot(self.x[~self.parabola.inlier_mask_], self.y[~self.parabola.inlier_mask_], '.', markersize=1.5, color='green')  
        plt.axis('off')   
        
################# Internal functions ###################
    def compute_metrics(self, model, verbose=True):
        inliers_boolean = self.parabola.inlier_mask_
        if model == "parabola":
            self.x_predicted = self.parabola.predict(self.y_quadratic_trans)
            # Compute residuals (include only inliers)
            residuals = self.x_reshaped[inliers_boolean] - self.x_predicted[inliers_boolean]
            ## Concavity
            concavity = abs(self.concavity)
            ## Median absolute residual (superior & inferior arcades)
            vertex_x, vertex_y = self.get_parabola_vertex()
            # Superior arcade (include only inliers)
            top_seg = self.mask[0:round(vertex_y),:]
            top_coords = vessel_pixel_coordinates(top_seg)
            top_residuals = residuals[0:len(top_coords)]
            top_median_residual = np.median(abs(top_residuals))
            # Inferior arcade  (include only inliers)
            bottom_seg = self.mask[round(vertex_y):,:]
            bottom_coords = vessel_pixel_coordinates(bottom_seg)
            bottom_residuals = residuals[-len(bottom_coords):]
            bottom_median_residual = np.median(abs(bottom_residuals))            
            
        elif model == "linear":
            residuals = []
            # First segment
            if len(self.y_inliers_1) > 0:
                x_predicted = self.linreg1.predict(self.y_inliers_1 - self.centre_factor_y)
                x_predicted = x_predicted.flatten() + self.centre_factor_x
                residuals = residuals + list(self.x_inliers_1.flatten() - x_predicted)
            # Second segment
            if len(self.y_inliers_2) > 0:
                x_predicted = self.linreg2.predict(self.y_inliers_2 - self.centre_factor_y)
                x_predicted = x_predicted.flatten() + self.centre_factor_x
                residuals = residuals + list(self.x_inliers_2.flatten() - x_predicted)
            residuals = np.array(residuals)
        else:
            raise nameError("model must be one of: 'parabola' or 'linear'")
            
        median_residual = np.median(abs(residuals))
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((self.x-np.mean(self.x))**2)
        r_squared = 1 - (ss_res / ss_tot)
                
        if verbose: 
            outliers_boolean = np.logical_not(inliers_boolean)
            print("#### RANSAC", model, ":", "Excluded", sum(outliers_boolean), "outliers out of", 
                  sum(inliers_boolean)+sum(outliers_boolean), "vessel pixels ####")
            try:
                print("Concavity index: {:.4f}".format(abs(concavity)))
                print("Top median residual: {:.4f}".format(abs(top_median_residual)))
                print("Bottom median residual: {:.4f}".format(abs(bottom_median_residual)))
            except:
                pass            
            print("Median residual: {:.4f}".format(median_residual))
            print("R2: {:.4f}".format(r_squared))

        
        try:
            return concavity, median_residual, top_median_residual, bottom_median_residual, r_squared, 
        except:
            return median_residual, r_squared

    def parabola_index(self):
        try:
            _, parabola_median_residual, parabola_r2 = self.compute_metrics(model="parabola", verbose=False)
            linear_median_residual, linear_r2 = self.compute_metrics(model="linear", verbose=False)
            median_residual_ratio = linear_median_residual / parabola_median_residual
            r2_ratio = parabola_r2 / linear_r2
            return median_residual_ratio, r2_ratio
        except:
            raise NotImplementedError("Both parabola and segmented linear models must be fitted before using this")
            
            
def vessel_pixel_coordinates(mask):
    # Get foreground (i.e. vessel) pixel coordinates 
    boolean_array = np.where(mask > 0)
    coords = np.flip(np.column_stack(boolean_array), axis=1)
    return coords                