%% This code implements an optical flow algorithm that is commonly used to 
%% solve 2D image registration problems. 
%% 
%% The goal of image registration is to match different sets of data into 
%% one coordinate system. An optical flow algorithm is an image registration 
%% technique which aims at estimating a dense velocity field assuming an 
%% apparent grey level intensity conservation during displacement. In other 
%% words, each spatio-temporal variation of grey level intensity is attributed 
%% to "motion". However, with a problem stated like that, an infinite number 
%% of velocity fields may potentially be a solution: any voxel-matching is 
%% a solution as long as the paired voxels have similar grey level intensities. 
%% An additional constraint is thus required. The method originally proposed 
%% by Horn and Schunck introduced an additional spatial regularization 
%% constraint by assuming that the motion field is smooth in the neighborhood 
%% of each estimation point.
%% 
%% This code contains 1 tutorial "test_Optical_Flow_L2L2_2D.m" designed to:
%% - load an abdominal MR-image series,
%% - estimate dense velocity field between 2 opposite time of the breathing 
%% cycle,
%% - diplay the image registration results.
%% 
%% This code has been written by Baudouin Denis de Senneville (Institut de 
%% Mathématiques de Bordeaux, UMR 5231 CNRS, Université de Bordeaux 351, 
%% cours de la Libération - F 33 405 TALENCE).
%% 
%% This code was developed under the commercial software Matlab ©1994-2021 
%% The MathWorks, Inc.
%% 
%% Please cite the following paper if you are using this code:
%% 
%% Zachiu C., Papadakis N., Ries M., Moonen C. T. W., Denis de Senneville B., 
%% An improved optical flow tracking technique for real-time MR-guided beam 
%% therapies in moving organs, Physics in Medicine and Biology, 2015; 60(23), 9003.

clc;
close all;
clear all;

addpath(pwd);

%%=========================  Read data file ================================================

[data,dimx,dimy,dimz,no_dyn] = load_dat('./data/abdomen2D.dat');

%%========================= Input parameters for the image registration algorithm ==========

%% Dynamic image used as the reference position
reference_dynamic      = 4;

%% Dynamic image used as the current position
current_dynamic        = 7;

%% Weighting factor between the data fidelity term and the regularisation term
%% Low: motion is highly sensitive to grey level intensity variation. 
%% High: the estimated motion is very regular along the space.
alpha                  = 0.1;

%% Number of resolution level in the multi-resolution scheme
nb_resolution_level    = 3;

%%========================= Adjustement of grey level intensities =========================

%% Get the reference image for the registration
Iref = flipud(data(:, :, :, reference_dynamic)');

%% Get the current image to register
Icur = flipud(data(:, :, :, current_dynamic)');

%% Normalize the reference image
Iref = (Iref - min(Iref(:)))/(max(Iref(:)) - min(Iref(:)));

%% Normalize the current image by adjusting the mean of the images
Icur = Icur * (mean(Iref(:))/mean(Icur(:)));

%% Estimate motion field
tic;
[Ireg, u, v] = MultiRes_Optical_Flow_2D(Iref, Icur, alpha, nb_resolution_level);
toc;

%% Display registered images & estimated motion field
display_result2D(Iref,Icur,Ireg,complex(u, v));
