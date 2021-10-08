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

function [u, v] = OF_HS_2D(Iref,Icur,alpha,uprec,vprec) %% Algo L2L2

[dimx,dimy] = size(Iref);

fprintf('Compute 2D Horn&Schunck algorithm on image resolution [ %d x %d ]\n',dimx,dimy);

nb_iter_FP_max = 1000;   % Max number of iterations in the fix point scheme
epsilon_pf = 0.0001;

u = zeros(dimx,dimy);
v = zeros(dimx,dimy);

%% Move the current image with the current velocity field estimate
Ireg = imageRegistration_2D(Icur,uprec+u,vprec+v);

uk = zeros(dimx,dimy);
vk = zeros(dimx,dimy);

error2 = zeros(nb_iter_FP_max);

kernelM = ones(3,3)/9.;

Ix = zeros(dimx,dimy);
Iy = zeros(dimx,dimy);

%% Compute the Laplacian of the current motion estimate
lapl_uprec = conv2(uprec + u, kernelM, 'same') - uprec - u;
lapl_vprec = conv2(vprec + v, kernelM, 'same') - vprec - v;

%% Compute the spatio-temporal gradient
Ix(2:end-1,:) = (Ireg(3:end,:) - Ireg(1:end-2,:))/2.;
Iy(:,2:end-1) = (Ireg(:,3:end) - Ireg(:,1:end-2))/2.;
It = Ireg - Iref;

%% Set boundary conditions on the spatial gradients
lapl_uprec(1,:) = lapl_uprec(2,:);
lapl_uprec(end,:) = lapl_uprec(end-1,:);
lapl_uprec(:,1) = lapl_uprec(:,2);
lapl_uprec(:,end) = lapl_uprec(:,end-1);
lapl_vprec(1,:) = lapl_vprec(2,:);
lapl_vprec(end,:) = lapl_vprec(end-1,:);
lapl_vprec(:,1) = lapl_vprec(:,2);
lapl_vprec(:,end) = lapl_vprec(:,end-1);
Ix(1,:) = Ix(2,:);
Ix(end,:) = Ix(end-1,:);
Ix(:,1) = Ix(:,2);
Ix(:,end) = Ix(:,end-1);
Ix(1,1) = Ix(2,2);
Ix(end,1) = Ix(end-1,2);
Ix(1,end) = Ix(2,end-1);
Ix(end,end) = Ix(end-1,end-1);
Iy(1,:) = Iy(2,:);
Iy(end,:) = Iy(end-1,:);
Iy(:,1) = Iy(:,2);
Iy(:,end) = Iy(:,end-1);
Iy(1,1) = Iy(2,2);
Iy(end,1) = Iy(end-1,2);
Iy(1,end) = Iy(2,end-1);
Iy(end,end) = Iy(end-1,end-1);

%% Compute IxIx IyIy
IxIx = Ix .^ 2;
IyIy = Iy .^ 2;
denom = alpha + IxIx + IyIy;

%% Loop for the fix point scheme
for it_pf = 1 : nb_iter_FP_max   % Algo L2L2
    
    utmp = uk;
    vtmp = vk;
    
    mean_x = conv2(uk, kernelM, 'same') + lapl_uprec;
    mean_y = conv2(vk, kernelM, 'same') + lapl_vprec;
    
    %% Correct zero padding of Matlab conv2
    mean_x(1,:) = mean_x(2,:);
    mean_x(end,:) = mean_x(end-1,:);
    mean_x(:,1) = mean_x(:,2);
    mean_x(:,end) = mean_x(:,end-1);
    mean_x(1,1) = mean_x(2,2);
    mean_x(end,1) = mean_x(end-1,2);
    mean_x(1,end) = mean_x(2,end-1);
    mean_x(end,end) = mean_x(end-1,end-1);
    mean_y(1,:) = mean_y(2,:);
    mean_y(end,:) = mean_y(end-1,:);
    mean_y(:,1) = mean_y(:,2);
    mean_y(:,end) = mean_y(:,end-1);
    mean_y(1,1) = mean_y(2,2);
    mean_y(end,1) = mean_y(end-1,2);
    mean_y(1,end) = mean_y(2,end-1);
    mean_y(end,end) = mean_y(end-1,end-1);
    
    phi = (mean_x .* Ix + mean_y .* Iy + It) ./ denom;
    
    uk = mean_x - Ix .* phi;
    vk = mean_y - Iy .* phi;
    
    residu = mean(sqrt((uk(:) - utmp(:)) .^ 2 + (vk(:) - vtmp(:)) .^2));
    
    error2(it_pf) = residu;
    
    if (residu < epsilon_pf)
        break;
    end
    
end

%% Update u and v
u = u + uk;
v = v + vk;

end
