clear

top_level_folder = uigetdir(); %select the folder containing all the image 
%folders
cd(top_level_folder)

folders = dir('*stack*'); 
number_of_img_folders  = size(folders,1);
%looks for any folder containing the work stack
%if you are looking for something else (eg. 24hrs) you will need to change
% the key word.


%for loop to loop through image folders and process each one.
for i=1:number_of_img_folders;

    shortname = folders(i).name; %shrotname is the name of the image_folders
    %the image name should start with the folder name for this to work
    
    fprintf('working with image %s\n', shortname)
    
    img_folder_path = strcat(top_level_folder, '\', folders(i).name);
    cd(img_folder_path);
    
    %Generating image channel names
    %Change the channel number if necessary testest
    dapi_name = strcat(shortname, '_iso_ch2.tif');
    macro_name = strcat(shortname, '_iso_ch4.tif');
    ves_name = strcat(shortname, '_iso_ch2.tif');
    part_name = strcat(shortname, '_iso_ch5.tif');
    macro_seg_name = strcat(shortname, '_post_processed_macrophage.tif');
    
    %load the images, I am using loadtiff but this is inter changeable for
    %the prefered tiff loader. readMultipageTiff and imreadfast work the
    %same without telling you how long it takes to load an image.
    
    
    macro_seg_im = loadtiff(macro_seg_name);
    particle_im = loadtiff(part_name);
    macro_im = loadtiff(macro_name);
    %you can use dipiamge threshold function if desired
    %threshmac = threshold(macro_seg_im);
    
    %threshold the segmeneted macrophage
    threshmac = imbinarize(macro_seg_im);

    %mask is an image of the particle only contained in the macrophages
    
    %bsxfun does the specified operation elementwise in this case
    %multiplication. 
    
    %cast confirms that the image are of the same data type and if they are
    %different datatypes it first matches the datatypes. Matching logical
    %to uint16
    
   [pearson_coef_coloc_im, manders_overlap, k1, k2, M1, M2, mask_im,...
        particle_threshold_value, particle_threshold_image, ...
        bkg_corrected_particle_im] = correlation_coeff_fun(macro_im, ...
        particle_im,macro_seg_im);

   csv_name =  strcat(shortname, '_Correlation_Values.csv');
   
   TOV = table(pearson_coef_coloc_im, manders_overlap, k1, k2, M1, M2, particle_threshold_value);
    
   %[pearson_coef_bkg, manders_overlap_bkg,k1_bkg, k2_bkg, M1_bkg, M2_bkg] = pearsons_coeff_fun ...
   %    (bkg_macro_im, bkg_par_im, macro_seg_im)
   
    dir_csv = strcat(img_folder_path,'\', "Correlation_Coefficient");
   
   if exist(dir_csv, 'dir')~=7
        mkdir(dir_csv);
   end
   
   csv_path = strcat(dir_csv, '\' , csv_name);
   cd(dir_csv);
   writetable(TOV, csv_name)
   
   
   mask_name = strcat(dir_csv, '\', shortname, '_masked_image.tiff');
   threshold_name = strcat(dir_csv, '\', shortname, '_NP_threshold_image.tiff');
   np_bkg_name = strcat(dir_csv, '\', shortname, '_NP_background_corrected.tiff');
   
   fprintf('Writing masked image to %s\n', mask_name)
   tic
   imwrite(uint16(mask_im(:,:,1)),mask_name);
   num_slices = size(mask_im,3);
    for p = 2:num_slices
        imwrite(uint16(mask_im(:,:,p)),mask_name, 'WriteMode','append');
    end
    
   fprintf('Writing threshold image to %s\n', threshold_name)
   
   imwrite(uint16(particle_threshold_image(:,:,1)),threshold_name);
   num_slices = size(particle_threshold_image,3);
    for p = 2:num_slices
        imwrite(uint16(particle_threshold_image(:,:,p)),threshold_name, 'WriteMode','append');
    end

   fprintf('Writing background corrected image to %s\n', np_bkg_name)
   
   imwrite(uint16(bkg_corrected_particle_im(:,:,1)),np_bkg_name);
   num_slices = size(bkg_corrected_particle_im,3);
    for p = 2:num_slices
        imwrite(uint16(bkg_corrected_particle_im(:,:,p)),np_bkg_name, 'WriteMode','append');
    end
    
   toc
   fprintf('done writing image \n')
 end    
    
