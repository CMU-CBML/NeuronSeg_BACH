function parameters = setparameters_neuron01_seg()

rho(1,1) = 0.01;
rho(2,1) = 0.05;
rho(3,1) = 1;

orderGauss = 6;

maxlevel = 3;

lambda1 = 0.01;
lambda2 = 0.01;
lambda3 = 0.01;

maxiteration = 10;

timestep(1,1) = 0.0003;
timestep(2,1) = 0.0001;
timestep(3,1) = 0.00005;

nelemx = 40;
nelemy = 40;

pU = 3;
pV = 3;

sigma =5;% the key parameter which needs to be tuned properly.
sigma_phi = 0.5;% the variance of regularized Gaussian kernel
epsilon = 1;
parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,...
    'orderGauss',orderGauss,'rho',rho,'timestep',timestep,'lambda1',lambda1,'lambda2',lambda2,'lambda3',lambda3,'maxiteration',maxiteration,...
    'sigma',sigma,'sigma_phi',sigma_phi,'epsilon',epsilon);
end