function [pearson_coef_coloc_im, manders_overlap, k1, k2, M1, M2, mask_im, ...
    particle_threshold_value, particle_threshold_image, bkg_corrected_particle_im] = ...
    correlation_coeff_fun(macro_im, particle_im,...
    macro_seg_im, varargin)
%%%
%Correlation_coeff_fun

% this function calculates the Pearcon Correlation Coefficient (PCC),
% Manders Overlap Coefficient (MOC) and Manders Correlation Coefficents
% (k1,k2, M1,M2)


% 
% %Args:
% macro_im (array / unit16): Can be a 2D or a 3D image of the macrophage channel

% particle_im (array / unit16): Image of the corresponding nanoparticle channel

% macro_seg_im (array/ unit16): this is the post processed or ilastik segmented 

% image of the macrophage channel all values should be either 0 or 1 regardless of type.

% varargin (int): Optional background value for images which can be specified if desired
% these are not mandatory input arguments. Function does not currently do
% anything with these values but functionality can be added later. 

%It is recommended that you handle any background correction before feeding 
%images into function
 
% Reurns:
% 
% pearson_coef_coloc_im (double): the PCC
% 
% manders_overlap (double): the overlap coefficent of the masked nanoparticle 
% image. 
% 
% k1 (double):"colocalization split" contribution of macro intenstiy to np intensity 
% 
% k2 (doouble): "colocalization split" contribution of particle signal to macro singal
% 
% M1 (double): MCC of macrophage channel with a nanoparticle mask (not useful for zach)
% 
% M2 (double): Manders Correlation Coefficient of masked image with particle channel
% 
% 
% 
% Note: Masked image =  image of only the nanoparticle single that overlaps 
% the segmented macrophage image.

% Additional values are generated for the whole image instead of just the
% colocalized parts but these values are not returned. Change function
% output values if you are looking for something else. (Ex  replace 
% pearson_coef_coloc_im with pearson_coef_masked_im)



%%%

%% Preparing to calculate coefficients 

%fprintf('the total number of inputs %d\n', nargin);

nVargins = length(varargin);

if nVargins == 2
    bkg1 = cell2mat(varargin(1));
    bkg2 = cell2mat(varargin(2));
    
    %add additional code to handle background correction as needed
end

%thresholding and making colocalization masks
threshmac =  imbinarize(macro_seg_im);
coloc_part_w_macro= bsxfun(@times, particle_im, cast(threshmac, 'like', particle_im));
[threshpar_im, particle_threshold_value] = background_correction(particle_im);
coloc_macro_w_part= bsxfun(@times, macro_im, cast(threshpar_im, 'like', macro_im));
mask_im = coloc_part_w_macro;
particle_threshold_image = threshpar_im;
bkg_corrected_particle_im = bsxfun(@times, particle_im, cast(threshpar_im, 'like',particle_im));

% takes the loaded images and reshaps array to 1D and datatype double
%double or equivalent datatype must be used such that we don't reach the
%max value in a unit16 image and alter the calculation 

par = double(reshape(particle_im,[],1)); %particle
segmac =double((reshape(macro_seg_im,[],1))); %segmented macrophages
mac = double(reshape(macro_im,[],1)); %macrophage
par_coloc_w_mac = double(reshape(coloc_part_w_macro,[],1));
macro_coloc_w_par = double(reshape(coloc_macro_w_part,[],1));


%find average value for each channel
mac_av = mean(mac); 
par_av = mean(par);
segmac_av = mean(segmac);
par_col_av = mean(par_coloc_w_mac);
macro_col_av = mean(macro_coloc_w_par);


%subtract average from original 
mac_sub = mac -mac_av;
par_sub =  par-par_av;
segmac_sub = segmac - segmac_av;
par_col_sub = par_coloc_w_mac - par_col_av;
macro_col_sub = macro_coloc_w_par - macro_col_av;


%find dot product of two average subtracted images

%This is the same as the sumation of every element multiplied by 
%corresponding element in secon  array, this is the top half of 
%pearson coefficient 

sum_mac_par = sum(mac_sub .* par_sub);
sum_segmac_par = sum(segmac_sub .* par_sub);
sum_par_col_w_mac = sum(par_col_sub .* par_sub);

%square every value in array
mac_sub_sq = mac_sub.^2;
par_sub_sq = par_sub.^2;
segmac_sub_sq = segmac_sub.^2;
par_col_sub_sq = par_col_sub.^2;

%% Pearson Correlation coefficient 

%Note: Generally speaking R = macro, G = particle. In some cases R is also
%the masked image containing particles inside the macrophages

%bottom half of pearson coefficient fraction
mac_sqrt = sqrt(sum(mac_sub_sq) * sum(par_sub_sq));
segmac_part_sqrt = sqrt(sum(segmac_sub_sq) *sum(par_sub_sq));
par_col_sqrt = sqrt(sum(par_col_sub_sq) * sum(par_sub_sq));


%calculate pearson coefficient 
coef_mac_and_par_ims = sum_mac_par/mac_sqrt; %orriginal images

coef_segmac_par_im = sum_segmac_par/ segmac_part_sqrt; %  segmented macro image
coef_par_col_im = sum_par_col_w_mac/ par_col_sqrt; % masked image of particles in macrophage

pearson_coef_whole_im = coef_mac_and_par_ims;
pearson_coef_segmac_im = coef_segmac_par_im; 
pearson_coef_coloc_im = coef_par_col_im;

%pearsons coefficent if we use the segmented macrophage image instead of

%% Manders overlap coefficient
% overlap coefficent will change depending on background correction
% also sensitive to signal instensity

%first get dot product of macrophage and particle channel top half of eq
overlap_top_mac_par = sum(mac .*par);
overlap_top_segmac_par = sum(segmac .* par);

%square the intensity values
mac_sq = mac .^2;
segmac_sq = segmac.^2;
par_sq = par.^2;

%caluclate the bottom half of equation
overlap_bot_mac_par = sqrt((sum(mac_sq))*(sum(par_sq)));
overlap_bot_segmac_par = sqrt((sum(segmac_sq))*(sum(par_sq)));

overlap_coef_mac_par = overlap_top_mac_par/ overlap_bot_mac_par;
overlap_coef_segmac_part = overlap_top_segmac_par/ overlap_bot_segmac_par;

%return coefficients 

%manders overlap  is the overlap of all the pixels in the macrophage
%channel with all the pixels in the particle image
manders_overlap_whole_im = overlap_coef_mac_par;

%overlap_coef_segmac is the overlap of all the pixels in the segmented
%macrophage image with the particle image. Here there are no intensity
%values associated with segmented image and overlap suffers 
overlap_coef_segmac = overlap_coef_segmac_part;


%this is the overlap coefficient for just the nanoparticles contained in
%the macrophage channel with the orriginal image. Essentially the percent
%of nanoparticles contained inside the macrophages.

manders_overlap = sum(par_coloc_w_mac .* par)/ sqrt((sum(par_coloc_w_mac.^2))*(sum(par.^2)));

%% Manders Correlation Coeffient

%first we apply a mask to get just the nanoparticle singnal inside the
%macropages

%M1 = Ri_coloc = Ri if Gi>0 and Ri_coloc = 0 if Gi=0
%M2 = Gi_coloc = Gi if Ri>0 abd Gi_coloc=0 if Ri=0


%depends on the products ofe intensity of th macro and particle channels
%these are sensititve to the difference in intensity between channels.

%k1 and k2 split the calue of colocalization into a pair of separate
%parameters. Essnetially used to determine each channels contribution to
%the colocalization coefficient.

%Whole image
k1_whole = (sum(mac .* par))/ (sum(mac_sq));
k2_whole = (sum(mac .* par))/ (sum(par_sq));

%k1 and k2 represent only the colocalization image. k1 should equal 1
%because all the pixels in par_coloc_w_mac are in par

%I.e because this it is an image of the particle channel
k1 = (sum(par_coloc_w_mac .* par))/ sum(par_coloc_w_mac.^2); 
k2 = (sum(par_coloc_w_mac .* par))/ sum(par.^2);

%Describes the contribution of each channel in coloc

%M1 = requires a particle mask. A basic one has been made but M2 answers
%describes the image where only particle signal inside the macrophages is
%examined.

M1 = (sum(macro_coloc_w_par))/(sum(mac));

% By definition
% M2 = sum(Gi_coloc)/ sum(Gi) where Gi_coloc = Gi if Ri>0 
% particle_coloc =  particle_int if macro_in >0 else 0 (i.e masked img)


M2 =  sum(par_coloc_w_mac)/ sum(par);


end
