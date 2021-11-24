matlab_folder = pwd;

% Met_seg_folder = 'Z:\Ben K\2018-06-13 ki67 Liver and Tumor 1hr 4 hr 24 hr';
Main_folder = uigetdir();
cd(Main_folder);
files = dir('*-*');

%dir_test = 'Z:\TumourImages\U87-GNP50nm\0.5h\MSC159-T-stack4-Nov29-2018\ves_thresh';
for  i = 1:size(files,1)
    
    dir_timepoint = strcat(Main_folder,'\',files(i).name);
    cd(dir_timepoint)

    tic

    %Gets image name 
    [~,shortfile] = fileparts(files(i).name); %Gets image name 

    %Displays image being anaylzed and gets the file directory
    display(['Post processing (segmentation, erode edge)... ' shortfile]) 
    main_dir_file_name = strcat(dir_timepoint,'\',shortfile);
    


    %Gets vessel image names then loads the images
    
    cd(dir_timepoint)

    vess_name = strcat(shortfile,'_vessel_channel_name.tif'); 
    vess_ch = loadtiff(vess_name); 
% makes sure a binary image is loaded
    vessel_seg_bin = vess_ch == 1;
    vessel_seg_crop = vessel_seg_bin;
%default parameter of erode process is 'sphere" and 1 pixels from the center,
%can be tweaked based on quality of the segmentation
    se_erode = strel('sphere',1);     
    vessel_seg_crop_er = imerode(vessel_seg_crop,se_erode);

    vessels_processed = vessel_seg_crop_er>0;
% final post processed file
    post_vessels = vessels_processed;

    
%save the post processed file
    save_dir_name = strcat(dir_timepoint,'\','Post processing images');


     if exist(save_dir_name, 'dir')~=7
            mkdir(save_dir_name);
     end
     
    cd(save_dir_name)
    
    %saves the ilatik segmentation image in the results folder. This is not
    img_vessel_name = strcat(shortfile,'_post_processed_vessels.tif');
    clear options;
                options.overwrite = true;
                options.compress = 'lzw';
                saveastiff(uint16(post_vessels), img_vessel_name, options);

    toc
    end
%end