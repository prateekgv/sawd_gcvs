%==========================================================================
% main.m
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB script is to detect gait
%                events, namely toe-off, heel-strike, and midstance.
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================

%==========================================================================
% Clear workspace
%==========================================================================
close all; clc; clear all;

%==========================================================================
% Add utilities path
%==========================================================================
addpath('./utils/');
addpath('./sawd/');

%==========================================================================
% Select dataset
%==========================================================================
path = './data/';
[ h5names, ss ] = select_dataset( path );
fpath = strcat(path,h5names{ss});

%==========================================================================
% Load datasets 
%==========================================================================
txt = ['Extract data from ', path];
disp(txt)
[u_l,u_r,q_proc_l,q_proc_r,t_l,t_r] = load_data(fpath);

%==========================================================================
% Load settings 
%==========================================================================
load_settings(fpath);

%==========================================================================
% Run the zero-velocity and tremor detector
%==========================================================================
disp('Run the detector for left foot data')
[zupt_l,fog_l,ff_l]=fog_zupt_detectors(u_l(1:6,:));

disp('Run the detector for right foot data')
[zupt_r,fog_r,ff_r]=fog_zupt_detectors(u_r(1:6,:));

%==========================================================================
% SAWD-aided gait segmentation
%==========================================================================
disp('Run sparsity-aided gait segmentation for left foot data')
[to_l,hs_l] = sawd_aided_gait_segmentation(u_l,zupt_l,'left');
plot_fftohs(-u_l(5,:),-ff_l,-to_l,-hs_l,zupt_l,'left',ss);

disp('Run sparsity-aided gait segmentation for right foot data')
[to_r,hs_r] = sawd_aided_gait_segmentation(u_r,zupt_r,'right');
plot_fftohs(-u_r(5,:),-ff_r,-to_r,-hs_r,zupt_r,'right',ss);

