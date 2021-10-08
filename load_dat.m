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

function [image,dimx,dimy,dimz,no_dyn] = load_dat(file_name)

disp('Loading data file');

%% Open the data file
f = fopen(file_name, 'r');

%% Extract the header size
header_size = fread(f, 1, 'int');

%% Extract the header data
header = fread(f, header_size, 'int');

%% Store image dimension parameters
dimx = header(1);
dimy = header(2);
if (header_size == 3)
  dimz   = 1;
  no_dyn = header(3);
end
if (header_size == 4)
  dimz   = header(3);
  no_dyn = header(4);
end

%% Extract data from file
image = fread(f, 'float');
image = reshape(image,dimx,dimy,dimz,no_dyn);

%% Close data file
fclose(f);
