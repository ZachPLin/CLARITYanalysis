


function [post_Macrophages, post_nanoparticle] = post_process_func(seg_Macrophages,seg_nanoparticle,save_Lymph,sample_name)

% set matlab folder to return at the end
matlab_folder = pwd;
% select folder with the data to read files
data_folder = uigetdir('','Select folder with image files'); 
cd(data_folder);
% input sample name for rename in the terminal
% sample_name = input('Enter the name of the sample: ','s');
[~,name,~]=fileparts(pwd);
sample_name = name; %strcat(cell2mat(inputdlg('Enter the name of the smple:')));
%select save directory
save_dir = uigetdir('','Select save folder');

% select and read Macrophage image files
Macrophage_file = uigetfile('','Select trained Macrophage image file');
Macrophage_info = imfinfo(Macrophage_file);
seg_Macrophages = imread(Macrophage_file,1);
for ii = 2 : size(Macrophage_info, 1)
    temp_Macrophages_tiff = imread(Macrophage_file, ii);
    seg_Macrophages = cat(3 , seg_Macrophages, temp_Macrophages_tiff);
end

num_slices = size(seg_Macrophages,3);

%display filename
shortfile = sample_name;
display(['Post-processing ' shortfile])




%%%%Load in Ilastik segmented image files

Macrophages_seg = seg_Macrophages;


tic

%%%%%Post-processing of Ilastik segmented Macrophage channel
%makes sure it is a binary image of macrophage channel
Macrophage_seg_bin = Macrophages_seg >= 1;
Macrophage_seg_crop = Macrophage_seg_bin%.*tissue_area_bin;

%creating threshold segmentation of macrophage channel
MacThresh = opening(Macrophage_seg_crop,4);

%dialate macrophages
%default parameter of dilation is 'sphere" and 1 pixels from the center,
%can be tweaked based on quality of the segmentation
se_dilate = strel('sphere',1);
Macrophage_seg_crop_dil = imdilate(MacThresh,se_dilate);

Macrophages_processed = uint32(Macrophage_seg_crop_dil);


%%%%%Final segmentation images

post_Macrophages = Macrophages_processed;

%save the post processed images
Macrophages_processed_name = strcat(shortfile,'_post_processed_Macrophages.tif');
cd(save_dir);
 imwrite(uint8(Macrophages_processed(:,:,1)),Macrophages_processed_name);
         
    for p = 2:num_slices
            imwrite(uint8(Macrophages_processed(:,:,p)),Macrophages_processed_name, 'WriteMode','append');
    end   




 toc
 
cd(matlab_folder);



display(['post-processing completed: ' shortfile])
 
end


