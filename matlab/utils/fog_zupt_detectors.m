%==========================================================================
% function [zupt,fog]=fog_zupt_detectors(u)
%==========================================================================
% @descirption: Detects zero-velocity and termor instances based on the
%               accelerometer and gyroscope data. 
% @author     : Prateek Gundannavar
% @date       : 06/07/18
%
% @input      
%             - u        3xN accelerometer data and 3xN gyroscope data
%
% @output
%             - zupt     1xN vector consisting of zero-velocity instances
%             - fog      1xN vector consisting of tremor instances
%
% @references 
%             - I. Skog, P. Handel, J. O. Nilsson and J. Rantakokko, 
%               "Zero-Velocity Detection—An Algorithm Evaluation," in 
%               IEEE Transactions on Biomedical Engineering, vol. 57, 
%               no. 11, pp. 2657-2666, Nov. 2010.
%
%             - G. V. Prateek, I. Skog, M. E. McNeely, R. P. Duncan, 
%               G. M. Earhart and A. Nehorai, "Modeling, Detecting, and 
%               Tracking Freezing of Gait in Parkinson Disease using 
%               Inertial Sensors," in IEEE Transactions on Biomedical 
%               Engineering.
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================
function [zupt,fog,ff]=fog_zupt_detectors(u)

global simdata;

zupt=zeros(1,length(u));
fog=zeros(1,length(u));
fog_zupt=zeros(1,length(u));
ff=nan(1,length(u));
ts_array = [];
W=simdata.Window_size;
fs=simdata.fs;

switch simdata.detector_type
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                    With Alt. Minimization                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'FOG-ZUPT'
        [T_zupt,T_fog,T_zupt_fog]=FOG_ZUPT(u);
        
        % Fix the edges of the detector statistics
        T_zupt_fog=[T_zupt_fog min(T_zupt_fog)*ones(1,floor(W))];
        T_zupt=[T_zupt min(T_zupt)*ones(1,floor(W))];
        T_fog=[T_fog min(T_fog)*ones(1,floor(W))];
        
        T_zupt_fog = T_zupt_fog(1:length(u));
        T_zupt = T_zupt(1:length(u));
        T_fog = T_fog(1:length(u));
        
        for k=1:length(T_zupt)
            % determine fog and zupt vs not stationary instancess
            if(T_zupt_fog(k)<=simdata.kappa)
                fog_zupt(k)=1;
            end 
            % determine if zupt or fog
            if(((T_fog(k)>simdata.gamma)) ...
                    && fog_zupt(k)==1)
                fog(k)=1;
                % determine flat-foot stage
                if ~isempty(ts_array) && (length(ts_array)/fs > 0.1)
                    [~,ind] = min(ts_array);
                    ff(1,k-length(ts_array)+ind) = u(5,k-length(ts_array)+ind);
                    ts_array = [];
                end
            elseif(((T_fog(k)<=simdata.gamma))...
                    && fog_zupt(k)==1)
                zupt(k)=1;
                ts_array = [ts_array,T_zupt(k)];
            else
                % determine flat-foot stage
                if ~isempty(ts_array) && (length(ts_array)/fs > 0.1)
                    [~,ind] = min(ts_array);
                    ff(1,k-length(ts_array)+ind) = u(5,k-length(ts_array)+ind);
                    ts_array = [];
                end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                  Zero-velocity detection                        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 'ZUPT'
        [T_zupt,T_fog,T_zupt_fog,T_theta]=ZUPT(u);
        
        for k=1:length(T_zupt)
            if (T_zupt(k)<simdata.gamma)
                zupt(k:k+W-1)=ones(1,W);
            end
        end
        
        % Fix the edges of the detector statistics
        T_zupt=[min(T_zupt)*ones(1,floor(W/2)) T_zupt min(T_zupt)*ones(1,floor(W/2))];
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                 Without Alt. Minimization                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'FOG'
        [T_zupt,T_fog,T_zupt_fog,T_theta]=FOG(u);
        
        for k=1:length(T_zupt)
            % determine fog and zupt vs not stationary instancess
            if(T_zupt_fog(k)<=simdata.kappa)
                fog_zupt(k)=1;
            end
            % determine if zupt or fog
            if(((T_fog(k)>simdata.gamma)) ...
                    && fog_zupt(k)==1)
                fog(k)=1;
            elseif(((T_fog(k)<=simdata.gamma))...
                    && fog_zupt(k)==1)
                zupt(k)=1;
            end
        end
        
        % Fix the edges of the detector statistics
        T_zupt_fog=[min(T_zupt_fog)*ones(1,floor(W/2)) T_zupt_fog min(T_zupt_fog)*ones(1,floor(W/2))];
        T_zupt=[min(T_zupt)*ones(1,floor(W/2)) T_zupt min(T_zupt)*ones(1,floor(W/2))];
        T_fog=[min(T_fog)*ones(1,floor(W/2)) T_fog min(T_fog)*ones(1,floor(W/2))];
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [T_zupt,T_fog,T_zupt_fog]= FOG_ZUPT(u)
%
%> @brief Proposed method
%>
%> @param[out]  T_fog      The test statistics of the fog_zupt detector
%> @param[out]  T_zupt     The test statistics of the zupt detector
%> @param[in]   u          The IMU data vector.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T_zupt,T_fog,T_zupt_fog]= FOG_ZUPT(u)

global simdata;

g=simdata.g;
sigma2_a=simdata.sigma_a^2;
sigma2_g=simdata.sigma_g^2;
W=simdata.Window_size;


N=length(u);
T_zupt=zeros(1,N-W+1);
T_fog=zeros(1,N-W+1);
T_zupt_fog=zeros(1,N-W+1);

for k=1:N-W+1
    
    % zupt detector
    ya_m= mean(u(1:3,k:k+W-1),2);
    
    % matrix of all readings 3xW
    Ya= u(1:3,k:k+W-1);
    Yw= u(4:6,k:k+W-1);
    
    % initial vector
    va_i= ya_m/norm(ya_m);
    
    % iterate alternative minimization
    iter=2;
    max_iter= 10;
    cost_val= zeros(max_iter,1);
    epsilon= 1e-3;
    %cost_temp= zeros(max_iter,1);
    
    while (iter <= max_iter)
        Ga_new= (bsxfun(@minus, Ya, g*va_i))*(bsxfun(@minus, Ya,g*va_i))';
        [V1,D1] = eigs(Ga_new);
        [~,index] = max(abs(diag(D1)));  % The maximum eigenvalue and its index
        ua_i= V1(:,index);               % The associated eigenvector in V
        %cost_temp(iter,1)= ua_i'*Ga_new*ua_i;
        %fprintf('%e\t', cost_temp(iter,1));
        
        va_i = (1/(g*W))*(sum(Ya,2)-ua_i);
        va_i= va_i/norm(va_i);
        Ga_new= (bsxfun(@minus, Ya, g*va_i))*(bsxfun(@minus, Ya,g*va_i))';
        cost_val(iter,1)= ua_i'*Ga_new*ua_i;
        %fprintf('%e\t %d\n', cost_val(iter,1),k);
        if (cost_val(iter,1)-cost_val(iter-1,1))<epsilon
            break;
        end
        iter=iter+1;
    end
     
    %to determine if the foot is stationary/foging vs moving.
    Ga_new= (bsxfun(@minus,Ya,g*va_i))*(bsxfun(@minus,Ya,g*va_i))';
    T_zupt_fog(k)= (1/sigma2_a)*ua_i'*Ga_new*ua_i;
    
    temp= bsxfun(@minus,Ya,g*ya_m/norm(ya_m));
    T_zupt(k)= (1/sigma2_g)*norm(Yw(:))^2 + (1/sigma2_a)*norm(temp(:))^2;
    
    %to determine fog, we use both accelerometer and gyroscope values
    T_fog(k)= T_zupt(k)-T_zupt_fog(k);
    

end


T_zupt= T_zupt./W;
T_fog= T_fog./W;
T_zupt_fog= T_zupt_fog./W;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion T = ZUPT(u)
%
%> @brief Function that runs the generalized likelihood test (SHOE detector). 
%>
%> @param[out]  T          The test statistics of the detector 
%> @param[in]   u          The IMU data vector.     
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T_zupt,T_fog,T_zupt_fog,T_theta]=ZUPT(u)

global simdata;

g=simdata.g;
sigma2_a=simdata.sigma_a^2;
sigma2_g=simdata.sigma_g^2;
W=simdata.Window_size;


N=length(u);
T_zupt=zeros(1,N-W+1);
T_fog=T_zupt;
T_zupt_fog=T_zupt;
T_theta=T_zupt;

for k=1:N-W+1
    % zupt detector
    ya_m= mean(u(1:3,k:k+W-1),2);
    
    % matrix of all readings 3xW
    Ya= u(1:3,k:k+W-1);
    Yw= u(4:6,k:k+W-1);
    
    temp= bsxfun(@minus,Ya,g*ya_m/norm(ya_m));
    T_zupt(k)= (1/sigma2_g)*norm(Yw(:))^2 + (1/sigma2_a)*norm(temp(:))^2;
    
end

T_zupt=T_zupt./W;

end



