function param = setparameters_olfactoryfiber()
%This function sets the parameters for the neuron segmentation

%Maximum number of iterations
maxiteration = 10;

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
m_var = 30;
n_var = 30;
o_var = 30;

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

rho(1,1) = 1; %level 2 refinement
rho(2,1) = 3; %level 3 refinement
rho(3,1) = 1;

timestep(1,1) = 0.001; 
timestep(2,1) = 0.0005;
timestep(3,1) = 0.0001;

sigma = 20;% the key parameter which needs to be tuned properly.
sigma_phi = 100;% the variance of regularized Gaussian kernel
epsilon = 10;

%make a struct variable param, with all the parameters
param = struct('pU',pU,'pV',pV,'pW',pW,'maxiteration',maxiteration,'maxlevel',maxlevel,'nelemU',nelemU,'nelemV',nelemV,'nelemW',nelemW,...
    'orderGauss',orderGauss,'kU',kU,'kV',kV,'kW',kW,'nobU',nobU,...
    'nobV',nobV,'nobW',nobW,'rho',rho,'timestep',timestep,...
    'lambda_1',lambda_1,'lambda_2',lambda_2,'lambda_3',lambda_3,...
    'sigma',sigma,'sigma_phi',sigma_phi,'epsilon',epsilon);
end