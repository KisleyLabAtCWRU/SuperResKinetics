function thd_map = LRG_SuperRes_LocalThrMap(im, e)

%LK 021914 - make function to spit out local threshold map to store the
%local threshold for identifying particles for each frame

% stolen from Bo's 'LRG_SuperRes_Particle_Identify.m' function

%% identify particles using pre-processed image data
% input: im: a single frame
% input: e.SNR_enhance = true || false, true if want to use SNR_booster
% function to increase the SNR ratio before analysis
% input: e.local_thd = true || false, true if want to use local background
% and threshold
% input: e.wide2 = 1, 1.5, 2, 2.4, 2.9, 3, ... 5, cut off distance to define a local maximum,
% defalst is 3 for wide field.
% input: e.Gauss_width = 1 to 3, Gaussian width of the PSF. default is 2.
% input: e.fitting = 'rc' || 'gs' || 'el', fitting funtion to be used, 
% radia symmetry:'rc', Gauss fitting: 'gs', Euler fitting: 'el'. default is
% 'rc'
% input: e.test = true || false, true if using test mode to generate a
% figure after particle identification. default is false.
%% initial conditions & input parameters
im = double(im);
try % use try statement to avoid input errors. If any input is missing, using default value
    if e.SNR_enhance == true
        im = SNR_booster(im);
    end
catch ME
    im = SNR_booster(im);
end
try
    local_thd = e.local_thd;% using local threshold (true) or not (false)
    wide2 = e.wide2;% cut off distance to define a local maximum
    Gauss_width = e.Gauss_width;% fitting regions will be 4*Gauss_width+1
    fitting = e.fitting;% specify the fitting algorithm to use
catch ME
    local_thd = false;
    wide2 = 3;
    Gauss_width = 2;
    fitting = 'rc';% radia symmetry:'rc', Gauss fitting: 'gs', Euler fitting: 'el'
end
try
    test = e.test;
catch ME
    test = false;
end
wide = floor(wide2);
n = 3; % how many std to add up as a threshold
% FWHM = 2.35*Gauss_width
%% calculate threshold map
w = 50;% the local region to calculate the local background and threshold
% usually 50X50 and shift by 25 is a good choice for a 512X512 image
[v h] = size(im);
count = zeros(v, h);% count store the times each pixel has contributed in local
% background calculation
bg = count; % record the local background
sd = count; % record the local standard deviation

for i = 1 : ceil(v / w * 2) - 1
    for j = 1 : ceil(h / w * 2) - 1
        % 1, select the local region
        % 2, sort the pixels based on their intensities
        % 3, find the 50% and 75% point as discussed in the paper
        % 4, calculate local background(bg), standard deviation(sd) and
        % count
        im_local = im(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1)));
        im_local = sort(im_local(:));
        n_loc = numel(im_local);
        bg(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1))) =...
            bg(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1))) +...
            im_local(round(n_loc/2));
        sd(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1))) =...
            sd(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1))) +...
            im_local(round(n_loc*0.5)) - im_local(round(n_loc*0.18));%sd = 0.82 to 0.5 of cumulative distribution
        count(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1))) =...
            count(1+w/2*(i-1):min(v,w/2*(i+1)), 1+w/2*(j-1):min(h,w/2*(j+1))) + 1;
    end % for j
end % for i
bg = bg ./ count;
sd = sd ./ count;
thd_map = bg + n * sd;% determine the local threshold

%}