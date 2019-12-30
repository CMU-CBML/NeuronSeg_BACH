function parameters = setparameters_neuron04_seg()

% rho for THB-spline refinement 
rho(1,1) = 0.01;
rho(2,1) = 0.05;
rho(3,1) = 1;

%gaussian integration order
orderGauss = 6;

%max level of resolution
maxlevel = 3;

%lambda in regularization
lambda1 = 0.01;
lambda2 = 0.01;
lambda3 = 0.01;

%max iterations
maxiteration = 50;

%time step for each level
timestep(1,1) = 0.0005;
timestep(2,1) = 0.0003;
timestep(3,1) = 0.0001;

%number of control grid elements for first level
nelemx = 40;
nelemy = 40;

%polynomial order for THB-spline
pU = 3;
pV = 3;


sigma =5; %G_sigma: the key parameter which needs to be tuned properly for local image fitting
sigma_phi = 0.5; %G_sigma_phi: the variance of regularized Gaussian kernel
epsilon = 1; %epsilon: for width of level set function
parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,...
    'orderGauss',orderGauss,'rho',rho,'timestep',timestep,'lambda1',lambda1,'lambda2',lambda2,'lambda3',lambda3,'maxiteration',maxiteration,...
    'sigma',sigma,'sigma_phi',sigma_phi,'epsilon',epsilon);
end