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
    display(['Pre processing ... ' shortfile]) 
    main_dir_file_name = strcat(dir_timepoint,'\',shortfile);
    
    %finds Edge_noise_image for analyzing 
   


    %Gets vessel, macrophage and NP image names then loads the images
    
    cd(dir_timepoint)
    MQ_name = strcat(shortfile,'_macrophage_name.tif');
    np_name = strcat(shortfile,'_nanoparticle_name.tif');
    vess_name = strcat(shortfile,'_vessel_name.tif'); %_Edge_removed_vessel_iso_ch1, %_iso_ch1_Probabilities %change accordingly

    %np_ch = loadtiff(np_name);
    vess_ch = loadtiff(vess_name); 
    MQ_ch = loadtiff(MQ_name); 
    np_ch = loadtiff(np_name); 
    
    %pre-processing to normalize intensity for all channels is the same 
    %pre-processing is used to better identify vessels and macrophages in
    %Ilastik segmentation
    %note: analysis of nanoparticles shoud use the raw nanoparticle channel
    %without the pre-processing
    vessels_single = single(vess_ch);
    vessels_local = vessels_single./(0.2*max(vessels_single(:))+(im2mat(gaussf(vessels_single,size(vessels_single,1)/10))));
    vessels_loglocal = mat2im(log(vessels_local+0.1));
    vessels_processed = uint8(stretch(vessels_loglocal));
    pre_vessels = vessels_processed;

    macrophages_single = single(MQ_ch);
    macrophages_local = macrophages_single./(0.2*max(macrophages_single(:))+(im2mat(gaussf(macrophages_single,size(macrophages_single,1)/10))));
    macrophages_loglocal = mat2im(log(macrophages_local+0.1));
    macrophages_processed = uint8(stretch(macrophages_loglocal));
    pre_macrophages = macrophages_processed;
    
    nanoparticle_single = single(np_ch);
    nanopartricles_local = nanoparticle_single./(0.2*max(nanoparticle_single(:))+(im2mat(gaussf(nanoparticle_single,size(nanoparticle_single,1)/10))));
    nanoparticles_loglocal = mat2im(log(nanopartricles_local+0.1));
    nanoparticles_processed = uint8(stretch(nanoparticles_loglocal));
    pre_nanoparticles = nanoparticles_processed;
    
    

    save_dir_name = strcat(dir_timepoint,'\','Pre processing images');


     if exist(save_dir_name, 'dir')~=7
            mkdir(save_dir_name);
     end
     
    cd(save_dir_name)
    
    %saves the pre processed image 
    img_vessel_name = strcat(shortfile,'_pre_processed_vessel.tif');
    img_macrophage_name = strcat(shortfile,'_pre_processed_macrophage.tif');
    img_nanoparticle_name = strcat(shortfile,'_pre_processed_nanoparticle.tif');
    clear options;
                options.overwrite = true;
                options.compress = 'lzw';
                saveastiff(uint16(pre_vessels), img_vessel_name, options);
                saveastiff(uint16(pre_macrophages), img_macrophage_name, options);
                saveastiff(uint16(pre_nanoparticles), img_nanoparticle_name, options);;

    toc
    end
%end