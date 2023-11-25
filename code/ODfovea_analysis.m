%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                  Author: Fabian SL Yii                  %
%               Email: fabian.yii@ed.ac.uk                %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all, close all;

%% Get names of all masked images
% specify root directory
root = "/Users/fabianyii/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Projects/UKB_full"; % Mac
% root = "C:\Users\fslyi\OneDrive - University of Edinburgh\Projects\UKB_full";                    % Windows

% path to OD segmentation masks
OD_masks_path = fullfile(root, "outputs", "OD", "masks");   
% path to fovea segmentation masks
fovea_masks_path = fullfile(root, "outputs", "fovea", "masks");   
% specify directory where the tabular data (OD/foveal parameters) will be saved
save_csv_path = fullfile(root, "outputs", "csv");

% Get mask names
mask_names = dir(OD_masks_path);
% Only want to include filenames, i.e. exclude directories
mask_names = mask_names(~[mask_names.isdir]);

% Read full refraction dataset which contains mean SER & mean corneal radius (in mm).
refraction = readtable(fullfile(root, "data", "cleaned_data_long_PM_cohort.csv" ));
foveaData = readtable(fullfile(root, "outputs", "csv", "fovea_intensity_data.csv" ));


%% Main analysis starts here

% Cell array to store derived OD/foveal parameters of interest 
 result = {"name", "age", "SER", "meanCornealRadius", "sph", "cyl", "disc_x", "disc_y"..., 
           "od_area", "adj_od_area", "major_length", "adj_major_length"...,
           "minor_length", "adj_minor_length", "orientation"...,
           "dist", "adj_dist", "vertical_angle", "OD_seg", "fovea_seg"; 
           [], [], [], [], [], [], [], [], [], [], [] ..., 
           [], [], [], [], [], [], [], [], []};

 f = waitbar(0, 'Starting');
 for i=1:length(mask_names)

    waitbar(i/length(mask_names), f, sprintf('Progress: %d %%', floor(i/length(mask_names)*100)));

    % name of the current mask
    mask_name = mask_names(i).name;

    % Get foveal centroid coordinates
    this_fovea = foveaData(strcmp(foveaData.fundus, mask_name),:);
    fovea_centroid = [this_fovea.fovea_x, this_fovea.fovea_y];

    % Get mean SER, CR, SPH & CYL corresponding to 
    % this image from the refraction dataset 
    this_refraction = refraction(strcmp(refraction.fundus, mask_name),:);
    age = this_refraction.age;
    ser = this_refraction.SER;
    cr  = this_refraction.meanCornealRadius;
    sph = this_refraction.sph;
    cyl = this_refraction.cyl;
    % save image name, SER, CR, SPH and CYL to their respective columns in
    % the result cell array.
    result{i+1, 1} = mask_name; 
    result{i+1, 2} = age; 
    result{i+1, 3} = ser;      
    result{i+1, 4} = cr;       
    result{i+1, 5} = sph;
    result{i+1, 6} = cyl;
        
    % Read the OD mask (2D binary mask; 0 or 255)
    OD_mask_gray = imread( fullfile(OD_masks_path, mask_name) );        % 560 x 560
    
    % Read the (fovea) masked image (3D colour masked image)
    fovea_masked_img = imread( fullfile(fovea_masks_path, mask_name) ); % 560 x 560 x 3
    % Convert to 2D grayscale
    fovea_mask_gray = rgb2gray(fovea_masked_img);                       % 560 x 560
    % Convert to binary (0 or 255)
    fovea_mask_gray = uint8(imbinarize(fovea_mask_gray)*255);           % 560 x 560 (binary)


    %% Connected component analysis 
    OD_cc = bwconncomp(OD_mask_gray);
    fovea_cc = bwconncomp(fovea_mask_gray);
    
    % Only compute OD parameters if OD mask is not empty
    if OD_cc.NumObjects ~= 0
        % Compute OD centroid coordinates, area, major and minor axis length and orientation for each connected component
        OD_stats = regionprops(OD_cc, 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
        % Extract OD centroid 
        OD_centroid = [OD_stats.Centroid(1), OD_stats.Centroid(2)];
    
        % Extract OD Area
        od_area = OD_stats.Area;                           % not adjusted for magnification
        adj_od_area = littmann(od_area, ser, cr, "True");  % adjusted for magnification using Littmann's formula (see function below)
      
        % Extract OD Major Axis Length
        major_length = OD_stats.MajorAxisLength;
        adj_major_length = littmann(major_length, ser, cr, "False");
    
        % Extract OD Minor Axis Length
        minor_length = OD_stats.MinorAxisLength;
        adj_minor_length = littmann(minor_length, ser, cr, "False");

        % Extract OD orientation (angle b/w x axis and OD major axis), 
        orientation = OD_stats.Orientation;
        
        % Save OD parameters to their respective columns in the result cell array.
        result{i+1, 7} = OD_centroid(1);               % OD x coordinate
        result{i+1, 8} = OD_centroid(2);               % OD y coordinate
        result{i+1, 9} = od_area;
        result{i+1, 10} = adj_od_area;
        result{i+1, 11} = major_length;
        result{i+1, 12} = adj_major_length;
        result{i+1, 13} = minor_length;
        result{i+1, 14} = adj_minor_length;
        result{i+1, 15} = orientation;
    end

    % Only compute OD-foveal distance and angle if OD and fovea masks are not empty
    if OD_cc.NumObjects ~= 0 & fovea_cc.NumObjects ~= 0
        % Compute EUCLIDIAN DISTANCE b/w OD centroid & fovea centroid
        dist     = sqrt(sum((OD_centroid - fovea_centroid) .^2)); 
        adj_dist = littmann(dist, ser, cr, "False"); 
    
        % Vertical angle between the disc and macula, i.e. angle
        % between the horizontal midline passing through the centroid
        % of the disc & the centroid of the macula.
        vertical_angle = vert_angle(OD_centroid, fovea_centroid);      % compute using Pythagorean theorem (see function below)
            
        % Save OD-foveal parameters to their respective columns in the result cell array.
        result{i+1, 16}  = dist;      
        result{i+1, 17}  = adj_dist;  
        result{i+1, 18}  = vertical_angle;
    end
    
    % Make a note if OD mask was empty
    if OD_cc.NumObjects == 0
        result{i+1, 19} = "failed";
    end
    
    % Make a note if fovea mask was empty
    if fovea_cc.NumObjects == 0
        result{i+1, 20} = "failed";
    end
               
end

%% Save result (cell array) as csv
% write cell array to csv
path = fullfile(save_csv_path, 'OD_fovea_parameters.csv');
writecell(result, path);


%% Internal functions %%
% Function: calculates object size by accounting for the effect of 
% magnification due to ametropia. Takes image_size, spherical 
% equivalent refraction (SER) and mean (b/w strong and week meridians) 
% corneal radius (CR, in mm), and returns true (object) size. This is based
% on Littmann's original magnification formula for fundus camera.
% READ: https://link.springer.com/content/pdf/10.1007/BF00175988.pdf
function true_size = littmann(image_size, SER, CR, area)
a = 0.01 + 0.00236 * (CR - 8);
b = 0.6126 + 0.0968 * (CR - 8);
c = 30.52 + 2.57 * (CR - 8);

q = (a*SER^2 - b*SER + c) / 100;

if area == "True"
    true_size = (1.37 * q)^2 * image_size;
else
    true_size = 1.37 * q * image_size;
end

end

% Function: takes centroid1 (optic disc) coordinate and centroid2 (macula) 
% coordinate, and returns disc-fovea angle which describes the vertical
% separation between the disc and fovea. Note that we compute the absolute
% difference between disc x and macula x coordinates to disregard the
% influence of right or left eye on the sign of the vertical angle. We want
% the sign of the computed vertical angle to reflect only the spatial
% relationship between the disc and macula in the Y plane.
% NEGATIVE vertical angle means disc is HIGHER than macula (rmb origin
% [0,0] starts in the top left corner!
function vert_angle = vert_angle(centroid1, centroid2)
x_distance = abs(centroid1(1) - centroid2(1)); 
y_distance = centroid1(2) - centroid2(2);
vert_angle = atand(y_distance / x_distance);

end



