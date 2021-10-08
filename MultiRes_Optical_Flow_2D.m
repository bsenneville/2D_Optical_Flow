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

function [Ireg, u, v, nb_iter] = MultiRes_Optical_Flow_2D(Iref, Icur, alpha, nb_resolution_level)

[dimx,dimy] = size(Iref);

u = double(zeros(dimx,dimy));
v = double(zeros(dimx,dimy));

%% Loop on multi-resolution levels
for k=nb_resolution_level-1:-1:0
    
    scale=2^k;
    dimxk=dimx/scale;
    dimyk=dimy/scale;
    
    if (nb_resolution_level>=0)
        
        %% Resize motion field to the current resolution
        u=double(imresize(double(u), [dimxk dimyk], 'cubic'))*2.;
        v=double(imresize(double(v), [dimxk dimyk], 'cubic'))*2.;
        
        %% Resize images to the curent resolution
        Irefk = imresize(Iref, [dimxk dimyk], 'cubic');
        Icurk = imresize(Icur, [dimxk dimyk], 'cubic');
        
    else
        Irefk = Iref;
        Icurk = Icur;
    end
    
    %% Perform the motion estimation process
    [uk, vk] = OF_HS_2D(Irefk,Icurk,alpha,u,v);
    
    %% Update motion estimates
    u = u + uk;
    v = v + vk;
    
end

%% Register the current image using the final motion estimate
Ireg = imageRegistration_2D(Icur,u,v);

end







