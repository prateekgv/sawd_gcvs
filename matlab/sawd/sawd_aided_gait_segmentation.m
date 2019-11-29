%==========================================================================
% function [toe_off,heel_stk] = sawd_aided_gait_segmentation(u, zupt, foot)
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB script is run the SAWD
%                application
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019, Released under the 3-Clause BSD license
%==========================================================================
function [toe_off,heel_stk] = sawd_aided_gait_segmentation(u, zupt, foot)

toe_off = nan(1,length(u(5,:))); 
heel_stk = nan(1,length(u(3,:))); 

global simdata;
fs = simdata.fs;
NFFT = simdata.NFFT;

% Determine the start and stop time of a gait cycle
n_zupt= zeros(1,length(zupt));
n_zupt(find(zupt==0))= 1;
[start_t,stop_t]= edge_detection(n_zupt);

% Remove those segments which are less than 0.25 seconds
ind = find((stop_t - start_t)/fs < 0.25);
start_t(ind) = []; stop_t(ind) = [];

% Init variables
Hf_gx_avg = [];
gx_avg = [];
k_app = [];
rmse = [];
run_time = [];

% Load template
load dwt_l.mat
load dwt_r.mat

% figure(1);
for ii= 1:length(start_t)
    
    t1 = start_t(ii);
    t2 = stop_t(ii);
    gx= u(5,t1:t2)';
    
    % Preprocessing the signal
    gx = -gx;                                                       % sign change
    gx = ((2*(gx-min(gx)))./(max(gx)-min(gx)))-1;                   % scaling [-1,1]
    y = interp1(linspace(0,1,length(gx)),gx,linspace(0,1,fs));      % linear interpolation
    gx_avg = [gx_avg; y];
    
    % Compute the Fourier Spectrum
    [Hf_gx,~] = freqz(y, NFFT, 'whole');
    Hf_gx_avg = [Hf_gx_avg, abs(Hf_gx).^2/NFFT];
    
    % Sparsity-assisted gait segmentation
    wc = 0.025; lam = 0.05; mu = 0.1;
    tic;
    [~, ~, k_gx, k_app] = sawd_L1(y, 2, wc, lam, mu, k_app);
    t_toc = toc;
    run_time = [run_time, t_toc]; 
    
    % Compute RMSE
    if strcmp(foot,'left')
        rmse_i = sqrt(mean((k_gx(:,2)-dwt_l).^2));
        rmse = [rmse, rmse_i];
    elseif strcmp(foot,'right')
        rmse_i = sqrt(mean((k_gx(:,2)-dwt_r).^2));
        rmse = [rmse, rmse_i];
    end
    
    % Check validity
    if rmse_i < 0.46
        binX = determine_gait_segments(gx);
        [toe_off,heel_stk,~] = label_gait_segments(binX,u(5,:),t1,toe_off,heel_stk);
    end

end

if (~simdata.train)
    fprintf('GCVS RMSE: Mean = %f \t Std. Dev = %f \n',mean(rmse),std(rmse));
    fprintf('GCVS TIME: Mean = %f \t Std. Dev = %f \n',mean(run_time),std(run_time));
end

end


%==========================================================================
% function [toe_off,heel_stk] = label_gait_segments(iGS,gx,t1,toe_off,heel_stk)
%==========================================================================
function [toe_off,heel_stk,iGS] = label_gait_segments(binX,gx,t1,toe_off,heel_stk)

% Gyro data
[ms_start,ms_stop] = edge_detection(binX);
iGS = [ms_start,ms_stop-1];

if ~isempty(iGS)
    if (size(iGS,1) >= 2)
        % first local maxima (toe-off event)
        [~,ind] = max(gx(t1+iGS(1,1):t1+iGS(1,2)));
        toe_off(1,t1+iGS(1,1)+ind-1) = gx(t1+iGS(1,1)+ind-1);
        
        % second local maxima (heel-strike event)
        max_val = 0; ss = 2;
        for ii=2:size(iGS,1)
            [max_v, max_i] = max(gx(t1+iGS(ii,1):t1+iGS(ii,2)));
            if max_v > max_val
                max_val = max_v; ind = max_i; ss = ii;
            end
        end
        heel_stk(1,t1+iGS(ss,1)+ind-1) = gx(t1+iGS(ss,1)+ind-1);
    end
end

end

%==========================================================================
% function [binX] = determine_gait_segments(en_gx)
%==========================================================================
function [binX] = determine_gait_segments(en_gx)

th_gx = 0.5;
binX = zeros(1,length(en_gx));
binX(en_gx < th_gx) = 1;

end

