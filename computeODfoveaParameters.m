%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Relative importance of different optic disc and foveal %
%    parameters for the prediction of refractive error"    %
%                                                         %
%                  Author: Fabian SL Yii                  %
%               Email: fabian.yii@ed.ac.uk                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all, close all;

%% Get names of all masked images
% specify root directory
root = '.....'; 
% path to masked images (Fig 1 in the manuscript; bottom)
masked_img_path = fullfile(root, '.....');   
% path to original RGB images
img_path = fullfile(root, '.....');
% specify directory where the tabular data (OD/foveal parameters) will be saved
save_csv_path = fullfile(root, '.....');

% Get names of masked images
addpath masked_img_path;
img_names = dir(masked_img_path);
% Only want to include filenames, i.e. exclude directories
img_names = img_names(~[img_names.isdir]);
 
% Read full refraction dataset which contains mean SER & mean corneal radius (in mm).
refraction = readtable(fullfile(root, '.....' ));


%% Main analysis starts here

% Cell array to store derived OD/foveal parameters of interest 
 result = {'name', 'age', 'ser', 'cr', 'sph', 'cyl', 'dist', 'adj_dist', ..., 
           'vertical_angle', 'od_area', 'adj_od_area', 'major_length'...,
           'adj_major_length', 'minor_length', 'adj_minor_length'...,
           'macula_intensity', 'macula_intensity_R', 'macula_intensity_G'..., 
           'macula_intensity_B', 'scaled_macula_intensity', 'scaled_macula_intensity_R'..., 
           'scaled_macula_intensity_G', 'scaled_macula_intensity_B'..., 
           'median_intensity_R', 'median_intensity_G', 'median_intensity_B', 'median_intensity'..., 
           'orientation', 'disc_x', 'disc_y', 'macula_x', 'macula_y', 'reliability'; 
           [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []..., 
           [], [], [], [], [], [], [], [], [], [], [], [], [], []};

for i=1:length(img_names)

    % name of the current image
    img_name = img_names(i).name;

    % Get mean SER, CR, SPH & CYL corresponding to 
    % this image from the refraction dataset 
    this_refraction = refraction(strcmp(refraction.image_names, img_name),:);
    age = this_refraction.x1st_visit_age;
    ser = this_refraction.mean_ser;
    cr  = this_refraction.mean_cr;
    sph = this_refraction.mean_sph;
    cyl = this_refraction.mean_cyl;
    % save image name, SER, CR, SPH and CYL to their respective columns in
    % the result cell array.
    result{i+1, 1} = img_name; 
    result{i+1, 2} = age; 
    result{i+1, 3} = ser;      
    result{i+1, 4} = cr;       
    result{i+1, 5} = sph;
    result{i+1, 6} = cyl;
        
    % Read the image and masked image (RGB) and convert them to grayscale
    full_I = imread(fullfile(img_path, img_name));
    R = full_I(:,:,1);                % image in the red channel
    G = full_I(:,:,2);                % image in the green channel
    B = full_I(:,:,3);                % image in the blue channel
    full_gray = rgb2gray(full_I);     % convert the image from RGB to grayscale
    % Read the masked image (i.e. showing only OD and macula)
    masked_I = imread( fullfile(masked_img_path, img_name) );
    masked_lab = rgb2lab(masked_I);
    masked_R = masked_I(:,:,1);       % masked image in the red channel
    masked_G = masked_I(:,:,2);       % masked image in the green channel
    masked_B = masked_I(:,:,3);       % masked image in the blue channel
    masked_gray = rgb2gray(masked_I); % convert the masked image from RGB to grayscale
    
    % Compute the median (background) pixel intensity of the entire image, 
    % excluding zero pixels 
    gray_mask =  masked_gray==0;
    R_mask = masked_R==0;
    G_mask = masked_G==0;
    B_mask = masked_B==0;
    gray_background = median(full_gray(full_gray.*uint8(gray_mask) ~= 0));
    red_background = median(R(R.*uint8(R_mask) ~= 0));
    green_background = median(G(G.*uint8(G_mask) ~= 0));
    blue_background = median(B(B.*uint8(B_mask) ~= 0));

    %% Connected component analysis (masked image)
    % There should be two connected components, i.e. OD and macula
    cc = bwconncomp(masked_gray);
    num_cc = cc.NumObjects;                        % number of connected components
    num_pixels = cellfun(@numel, cc.PixelIdxList); % number of pixels in each component
    [largest, idx] = maxk(num_pixels, 2);          % order connected components by their size
        
    % Index for the largest component (should be OD) 
    OD_idx = idx(1); 
    % index for the second largest component (should be macula)
    macula_idx = idx(2);
    % Compute centroid coordinates, area, major and minor axis length and orientation for each connected component
    stats = regionprops(cc, 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
    % Macula centroid
    centroid2 = [stats(macula_idx).Centroid(1), stats(macula_idx).Centroid(2)];
    % OD centroid
    centroid1 = [stats(OD_idx).Centroid(1), stats(OD_idx).Centroid(2)];
    % OD Area
    od_area = stats(OD_idx).Area;                      % not adjusted for magnification
    adj_od_area = littmann(od_area, ser, cr, "True");  % adjusted for magnification using Littmann's formula (see function below)
    % OD Major Axis Length
    major_length = stats(OD_idx).MajorAxisLength;
    adj_major_length = littmann(major_length, ser, cr, "False");
    % OD Minor Axis Length
    minor_length = stats(OD_idx).MinorAxisLength;
    adj_minor_length = littmann(minor_length, ser, cr, "False");
    % EUCLIDIAN DISTANCE b/w OD centroid & macula centroid (fovea)
    dist     = sqrt(sum((centroid1 - centroid2) .^2)); 
    adj_dist = littmann(dist, ser, cr, "False"); 
    % Vertical angle between the disc and macula, i.e. angle
    % between the horizontal midline passing through the centroid
    % of the disc & the centroid of the macula.
    vertical_angle = vert_angle(centroid1, centroid2); % computed using Pythagorean theorem (see function below)
    % Median pixel intensity in the macular region
    macula_intensity = median(full_gray(createCirclesMask(full_gray, centroid2, 10))); % grayscale
    macula_intensity_R = median(R(createCirclesMask(full_gray, centroid2, 10)));       % red channel
    macula_intensity_G = median(G(createCirclesMask(full_gray, centroid2, 10)));       % green channel
    macula_intensity_B = median(B(createCirclesMask(full_gray, centroid2, 10)));       % blue channel
    % Scaled macular pixel intensity (adjusted for background median pixel intensity
    % of the entire grayscale image)
    scaled_macula_intensity = double(macula_intensity) - double(gray_background);
    scaled_macula_intensity_R = double(macula_intensity_R) - double(red_background);
    scaled_macula_intensity_G = double(macula_intensity_G) - double(green_background);
    scaled_macula_intensity_B = double(macula_intensity_B) - double(blue_background);
    % OD orientation (angle b/w x axis and OD major axis), 
    orientation = abs(stats(OD_idx).Orientation);
            
    % Save OD- and macular/foveal parameters to their respective
    % columns in the result cell array.
    result{i+1, 7}  = dist;      
    result{i+1, 8}  = adj_dist;  
    result{i+1, 9}  = vertical_angle;
    result{i+1, 10} = od_area;
    result{i+1, 11} = adj_od_area;
    result{i+1, 12} = major_length;
    result{i+1, 13} = adj_major_length;
    result{i+1, 14} = minor_length;
    result{i+1, 15} = adj_minor_length;
    result{i+1, 16} = macula_intensity;
    result{i+1, 17} = macula_intensity_R;
    result{i+1, 18} = macula_intensity_G;
    result{i+1, 19} = macula_intensity_B;
    result{i+1, 20} = scaled_macula_intensity;
    result{i+1, 21} = scaled_macula_intensity_R;
    result{i+1, 22} = scaled_macula_intensity_G;
    result{i+1, 23} = scaled_macula_intensity_B;
    result{i+1, 24} = red_background;             % median background intensity of the entire image in the red channel
    result{i+1, 25} = green_background;           % median background intensity of the entire image in the green channel
    result{i+1, 26} = blue_background;            % median background intensity of the entire image in the blue channel
    result{i+1, 27} = gray_background;            % median background intensity of the entire grayscale image
    result{i+1, 28} = orientation;
    result{i+1, 29} = centroid1(1);               % OD x coordinate
    result{i+1, 30} = centroid1(2);               % OD y coordinate
    result{i+1, 31} = centroid2(1);               % macula x coordinate
    result{i+1, 32} = centroid2(2);               % macula y coordinate
        
end

%% Save result (cell array) as csv
% write cell array to csv
path = fullfile(save_csv_path, '.....');
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



