function param = setparameters_helix()

%This function sets the parameters for the neuron segmentation

%--Input Variables

%--Output Variables

%param: struct variable storing the parameters for the neuron
%segmentation consisting of the following fields
%maxiteration: number of maximum iterations for each refinement level
%pU, pV, pW: B-spline degree in each parametric direction
%maxlevel: number of image resolution levels
%maxlevelg: number of refinement levels
%orderGauss: Gaussian quadrature order
%m_var, n_var, o_var: number of elements in each direction at the first
%refinement level
%lambda1, lambda2, lambda3: regularization parameters
% nelemU, nelemV, nelemW: number of elements at each refinement levels
% nobU, nobV, nobW: number of splines in each direction at each level
%kU, kV, kW: number of knot vectors in each refinement level
%rho: threshold parameter for the refinement of B-splines
%timestep: time step value for each refinement level
%sigma: G_sigma: the key parameter which needs to be tuned properly for local image fitting
%sigma_phi:%G_sigma_phi: the variance of regularized Gaussian kernel
%epsilon:width of level set function

%Maximum number of iterations
maxiteration = 1;

%Degree of B-splines
pU = 2;
pV = 2;
pW = 2;

%Maximum number of resolution levels
maxlevel = 3;
maxlevelg = 3;

%Gaussian quadrature order
orderGauss = 4;

%initial grid: number of elements in x,y,z direction for level 1
m_var = 10;
n_var = 10;
o_var = 10;

%regularization parameters
lambda_1 = 0.01;
lambda_2 = 0.01;
lambda_3 = 0.01;

%number of elements in each direction at each level
nelemU = zeros(maxlevelg,1);
nelemV = zeros(maxlevelg,1);
nelemW = zeros(maxlevelg,1);

%number of splines in each direction at each level
nobU = zeros(maxlevelg,1);
nobV = zeros(maxlevelg,1);
nobW = zeros(maxlevelg,1);

%number of knot vectors in each direction at each level
kU = zeros(maxlevelg,1);
kV = zeros(maxlevelg,1);
kW = zeros(maxlevelg,1);

rho = zeros(maxlevelg,1); %refinement parameter
timestep = zeros(maxlevelg,1); %timestep for each refinement level

%number of elements in each direction
for level = 1:maxlevelg
    nelemU(level,1) = m_var*2^(level-1);
    nelemV(level,1) = n_var*2^(level-1);
    nelemW(level,1) = o_var*2^(level-1);
    
    kU(level,1) = nelemU(level,1)+2*pU+1;
    kV(level,1) = nelemV(level,1)+2*pU+1;
    kW(level,1) = nelemW(level,1)+2*pU+1;
    
    nobU(level,1) = kU(level,1) - pU - 1;
    nobV(level,1) = kV(level,1) - pV - 1;
    nobW(level,1) = kW(level,1) - pW - 1;
end

rho(1,1) = 0.5; %level 2 refinement
rho(2,1) = 1.0; %level 3 refinement
rho(3,1) = 1;

timestep(1,1) = 0.0004;
timestep(2,1) = 0.00008;
timestep(3,1) = 0.00001;

sigma = 20;% the key parameter which needs to be tuned properly.
sigma_phi = 20;% the variance of regularized Gaussian kernel
epsilon = 0.25;

%make a struct variable param, with all the parameters
param = struct('pU',pU,'pV',pV,'pW',pW,'maxiteration',maxiteration,'maxlevel',maxlevel,'nelemU',nelemU,'nelemV',nelemV,'nelemW',nelemW,...
    'orderGauss',orderGauss,'kU',kU,'kV',kV,'kW',kW,'nobU',nobU,...
    'nobV',nobV,'nobW',nobW,'rho',rho,'timestep',timestep,...
    'lambda_1',lambda_1,'lambda_2',lambda_2,'lambda_3',lambda_3,...
    'sigma',sigma,'sigma_phi',sigma_phi,'epsilon',epsilon);
end