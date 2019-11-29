%==========================================================================
% function load_settings()
%==========================================================================
% @descirption: Creates a global variable called simdata. Stores all the
%               constants required the SAWD algorithm and ZVEI/TREI
%               detector
% @author     : Prateek Gundannavar
% @date       : 03/11/19
%
% @input      
%             - path        location of the APDM .h5 files
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================
function load_settings(path)

global simdata;

%==========================================================================
% APMD sensor parameters
%==========================================================================

% Sampling rate (Hz)
info = h5info(path,'/Sensors');
groups = info.Groups;
fname = groups(1).Name;
simdata.fs = double(h5readatt(path,strcat(fname,'/Configuration'),...
    'Sample Rate'));

% Sampling interval
simdata.ts = 1/simdata.fs;

%==========================================================================
% Location based constants
%==========================================================================

% Rough altitude [m] 
simdata.altitude=100;

% Rough latitude [degrees]
simdata.latitude=38;

% Magnitude of the local gravity vector [m/s^2]
simdata.g=gravity(simdata.latitude,simdata.altitude);

%==========================================================================
% Detector parameters
%==========================================================================

% Detector type to be used. You can chose between: 
% ZUPT     - zero-velocity detector
% FOG-ZUPT - freezing of gait and zero-velocity detector
% FOG      - Only freezing of gait (no splitting of gravity and inertial)
simdata.detector_type='FOG-ZUPT';

% Standard deviation of the acceleromter noise [m/s^2]. This is used to 
% control the zero-velocity detectors trust in the accelerometer data.
simdata.sigma_a=1; 

% Standard deviation of the gyroscope noise [rad/s]. This is used to 
% control the zero-velocity detectors trust in the gyroscope data.
simdata.sigma_g=0.8;

% Window size of the zero-velocity detector [samples] 
simdata.Window_size=simdata.fs/8;

% Threshold used in the zero-velocity detector. If the test statistics are 
% below this value the zero-velocity hypothesis is chosen.  
simdata.gamma=2.00; 

% Threshold used in the trembling detector. If the test statistics are 
% below this value the trembling hypothesis is chosen.  
simdata.kappa=34.53; 

%==========================================================================
% Sparsity-assisted gait segmentation
%==========================================================================
simdata.NFFT = 512;               % FFT length to compute spectrum
simdata.spe = 0.95;               % Spectral edge frequency

simdata.plot = false;              % Plot signal decompostion
simdata.train = false;             % Generate template DWT coefficients
end



%% SUBFUNCTIONS

%==========================================================================
%  funtion g=gravity(lambda,h) 
%==========================================================================
%
% @description: Function for calculating the magnitude of the local gravity. 
%
% @input      
%             - lambda     latitude [degrees] 
%             - h          altitude [m]
% @output
%             - g          gravit at given location (m/sec^2)
%==========================================================================
function g=gravity(lambda,h)

lambda=pi/180*lambda;
gamma=9.780327*(1+0.0053024*sin(lambda)^2-0.0000058*sin(2*lambda)^2);
g=gamma-((3.0877e-6)-(0.004e-6)*sin(lambda)^2)*h+(0.072e-12)*h^2;

end
