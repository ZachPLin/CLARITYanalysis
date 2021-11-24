matlab_folder = pwd;

% sepcify image folder;
Main_folder = uigetdir();
cd(Main_folder);
files = dir('*-*');

%analysis start
%read all file name in the fold
for  i = 1:size(files,1)
    
    dir_timepoint = strcat(Main_folder,'\',files(i).name);
    cd(dir_timepoint)

    tic

    %Gets image name 
    [~,shortfile] = fileparts(files(i).name); %Gets image name 

    %Displays image being anaylzed and gets the file directory
    display(['Analyzing NPs around vessels for ' shortfile]) 
    main_dir_file_name = strcat(dir_timepoint,'\',shortfile);
    
    %finds post_processed_image for analyzing 
   
    post_pro_dir = strcat(dir_timepoint, '\', 'Post processing images');
    %vess_thresh_dir = strcat(main_dir_file_name,'\ves_thresh');
    %cd(vess_thresh_dir)
    
    %finds ilastik segmented image in main folder and loads image
    cd(post_pro_dir)
    vess_thresh_name = strcat(shortfile,'_vessels-filename.tif'); 
    ves_thresh_orig = loadtiff(vess_thresh_name);
    ves_thresh = changem(ves_thresh_orig,[1],[255]);
    %Gets vessel and NP image names then loads the images
    
    cd(dir_timepoint)
    %finds the original nanoparticle channel
    np_name = strcat(shortfile,'_nanoparticle_channel-filename.tif'); 

    np_ch = loadtiff(np_name);


    
    %Locates microscope metadata, reads the metadata and gets the
    %conversion from pixels into micrometer.
    
    metadata_name = strcat(shortfile,'_iso_info.csv');
    metadata = readtable(metadata_name,delimitedTextImportOptions);
    px_per_um = str2num(metadata.ExtraVar1{2})*1E6;

    %makes binary image of segmented vessels

    vess_bin = ves_thresh==1;

    %creates a distance transform of the binary vessel image
    vess_dist = bwdist(vess_bin);


    %specifies the size of the z_pixel ( this is the lightsheet thickness)
    z_physsize = px_per_um;


    %Takes the 3D vessel distance transform array and makes it into 1D and
    %rounds to the nearest whole number.
    
    %works by scanning 3D array row by row for every slice. 
    %I.e linear_dist(3)= vess_dist(1,3,1) 
    
    linear_dist = round(single(vess_dist(:)));
    % linear_lumen_dist = round(single(lumen_dist(:)));
    
    %Linear_all concatenates np array to an arrray of ones the size of
    %linear_dist
    
    linear_all = cat(2, ones(size(linear_dist)), (np_ch(:)));
    % linear_lumen_all = cat(2, ones(size(linear_lumen_dist)), (np_ch(:)));
    
    %Creates an 3D array of zeros with rows from 0 to 250 inclusive. (250
    %is arbitrary and is the pixel distance away from the vessels that we
    %are searching in). part_vs_dist_raw and pixeldistances should be
    %changed to meet the distance criteria
    
    %Column 1 in the array will display the conversion of pixel distances 
    %into micrometer distances 
    
    %Column 2:  the number of voxels in the image at a given distance
    
    %Column 3:  Sum of NP intensity at a given distance region
   
    part_vs_dist_raw = zeros(201,3);

    pixeldistances = 0:1:200;

    
    %for loop which looks at a given pixeldistance from 1:250.
    %asks logical questions where does linear distance = given
    %pixeldistance
    
    % reassigns corresponding row values in part_vs_dist_raw
    
    %Recall that linear_all = (ones, np_values) in the last line of the for 
    %loop we are reassigning the ones column to the sum of the 
    %temp_locations which is a logical array obtained from asking where
    %does linear_dist == pixeldistances. We also obtain the corresponding
    %NP values
    
    for c = 1: size(pixeldistances,2) 
        temp_locations = linear_dist==pixeldistances(c);
        part_vs_dist_raw(c,:) = [pixeldistances(c)*z_physsize sum(linear_all(temp_locations,:),1)];
    end
    
    % for c = 1: size(lumenpixeldistance,2) 
    %     temp_locations = linear_lumen_dist==lumenpixeldistance(c);
    %     part_vs_lumen_raw(c,:) = [lumenpixeldistance(c)*z_physsize sum(linear_lumen_all(temp_locations,:),1)];
    % end
    
    %Appends a 4th column to the matrix for futher storing results of
    %processing data.
    
    part_vs_dist = [part_vs_dist_raw zeros(201,1)];

    
    %Convert into mean NP intensity (divide by volume)
    part_vs_dist(:,3) = part_vs_dist(:,3)./part_vs_dist(:,2);

    
    %Remove background noise - set lowest conc to 0
    part_vs_dist(:,3) = part_vs_dist(:,3)-min(part_vs_dist(:,3));  
    %Remove backgroud noise in lumen based on the minimum noise from
    %outside the lumen

    
    %Normalize relative to max (Note, this is not useful to use across
    %different time points or images).
    
    part_vs_dist(:,4) = part_vs_dist(:,3)./max(part_vs_dist(:,3));  

    
    %Plots distance away from vessels in um (x-axis)and corresponding NP 
    %intensity at that distance on the y-axis.
    
    plot(part_vs_dist(:,1), part_vs_dist(:,3),'-r');

    %Saves results in a new subfolder called results if the results folder
    %does not already exsit.
    save_dir_name = strcat(main_dir_file_name,'\Results_name');


     if exist(save_dir_name, 'dir')~=7
            mkdir(save_dir_name);
     end
     
    cd(save_dir_name)
    
    %saves the ilatik segmentation image in the results folder. This is not
    %necessary further processing is done to the image. 
    img_skel_name = strcat(shortfile,'_vess_thresh_ilastik.tif');
    clear options;
                options.overwrite = true;
                options.compress = 'lzw';
                saveastiff(uint16(vess_bin), img_skel_name, options);

    results_table = table(part_vs_dist);
    table_name = strcat(shortfile,'_Results_Vessel_analysis_ilastik_dist100_nonorm_bkrndsub.csv');

    %writetable(results_table,table_name);
    
    labled_table = array2table(part_vs_dist, 'VariableNames', {'Distance from Vessel (um)', 'Number of Voxels', 'Mean Nanoparticle Intensity', 'Normalized Nanoparticle Intensity'})
    labled_name = strcat(shortfile, '_Distance_analysis_results.csv')
    writetable(labled_table, labled_name)
    
    
    toc
    end
%end
