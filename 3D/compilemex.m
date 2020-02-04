path = pwd;
addpath(strcat(path,'/data'));
addpath(strcat(path,'/BsplineComposition'));
addpath(strcat(path,'/IterationLoopFunctions'));
addpath(strcat(path,'/LevelSetFunctions'));
addpath(strcat(path,'/PostProcessing'));
addpath(strcat(path,'/Initialization'));
addpath(strcat(path,'/SetParameters'));
addpath(strcat(path,'/THBSpline'));
disp('Added paths....');

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
%cfg.ReportPotentialDifferences = false;
cfg.GlobalDataSyncMethod = 'NoSync';

%% Define argument types for entry-point 'img2coef3D'.
ARGS1 = cell(1,1);
ARGS1{1} = cell(1,4);
ARGS1{1}{1} = coder.typeof(0,[1 Inf],[0 1]);
ARGS1{1}{2} = coder.typeof(0);
ARGS1{1}{3} = coder.typeof(0);
ARGS1{1}{4} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg img2coef3D -args ARGS1{1} -o img2coef3D

%% Define argument types for entry-point 'BsplineComposeImage3D'.
ARGS2 = cell(1,1);
ARGS2{1} = cell(7,1);
ARGS2{1}{1} = coder.typeof(0);
ARGS2{1}{2} = coder.typeof(0);
ARGS2{1}{3} = coder.typeof(0);
ARGS2{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS2{1}{5} = coder.typeof(0,[1 Inf],[0 1]);
ARGS2{1}{6} = coder.typeof(0,[1 Inf],[0 1]);
ARGS2{1}{7} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg BsplineComposeImage3D -args ARGS2{1} -o BsplineComposeImage3D

%% Define argument types for entry-point 'BsplineComposeImage3D_single'.
ARGS3 = cell(1,1);
ARGS3{1} = cell(7,1);
ARGS3{1}{1} = coder.typeof(single(0), [Inf,4,4],[0,1]);
ARGS3{1}{2} = coder.typeof(single(0), [Inf,4,4],[0,1]);
ARGS3{1}{3} = coder.typeof(single(0), [Inf,4,4],[0,1]);
ARGS3{1}{4} = coder.typeof(0, [Inf,Inf,Inf],[0,1]);
ARGS3{1}{5} = coder.typeof(0);
ARGS3{1}{6} = coder.typeof(0);
ARGS3{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg BsplineComposeImage3D_single -args ARGS3{1} -o BsplineComposeImage3D_single

%% Define argument types for entry-point 'computenewPoints'.
ARGS4 = cell(1,1);
ARGS4{1} = cell(7,1);

ARGS4{1}{1} = struct;
ARGS4{1}{1}.nzsplines = coder.typeof(int64(0), [Inf,1],[1,0]);
ARGS4{1}{1} = coder.typeof(ARGS4{1}{1}, [Inf,1],[1,0]);

ARGS4{1}{2} = coder.typeof(double(0), [Inf,4],[1,0]);

ARGS4{1}{3} = struct;
ARGS4{1}{3}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{3} = coder.typeof(ARGS4{1}{3}, [Inf,1],[1,0]);

ARGS4{1}{4} = struct;
ARGS4{1}{4}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{4} = coder.typeof(ARGS4{1}{4}, [Inf,1],[1,0]);

ARGS4{1}{5} = struct;
ARGS4{1}{5}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{5} = coder.typeof(ARGS4{1}{5}, [Inf,1],[1,0]);

ARGS4{1}{6} = struct;
ARGS4{1}{6}.mat = coder.typeof(single(0), [Inf,4,4,4],[1,0,0,0]);
ARGS4{1}{6} = coder.typeof(ARGS4{1}{6}, [Inf,1],[1,0]);

ARGS4{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg computenewPoints -args ARGS4{1} -o computenewPoints

%% Define argument types for entry-point 'tripleIterLoop'.
ARGS5 = cell(1,1);
ARGS5{1} = cell(4,1);
ARGS5{1}{1} = coder.typeof(0,[1 3]);

ARGS5{1}{2} = struct;
ARGS5{1}{2}.active_cell = coder.typeof(0);
ARGS5{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS5{1}{2} = coder.typeof(ARGS5{1}{2},[Inf  1],[1 0]);

ARGS5{1}{3} = struct;
ARGS5{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS5{1}{3} = coder.typeof(ARGS5{1}{3},[Inf  1],[1 0]);
ARGS5{1}{4} = coder.typeof(0,[Inf  4],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop -args ARGS5{1}

%% Define argument types for entry-point 'compute_Integ_Domain_fidelity'.
ARGS6 = cell(1,1);
ARGS6{1} = cell(11,1);
ARGS6{1}{1} = struct;
ARGS6{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS6{1}{1} = coder.typeof(ARGS6{1}{1},[Inf  1],[1 0]);

ARGS6{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS6{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS6{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS6{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);

ARGS6{1}{6} = coder.typeof(0,[Inf  4],[1 0]);

ARGS6{1}{7} = struct;
ARGS6{1}{7}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{7} = coder.typeof(ARGS6{1}{7},[Inf  1],[1 0]);

ARGS6{1}{8} = coder.typeof(0,[Inf  1],[1 0]);
ARGS6{1}{9} = coder.typeof(0,[Inf  1],[1 0]);
ARGS6{1}{10} = coder.typeof(0,[Inf  1],[1 0]);

ARGS6{1}{11} = coder.typeof(single(0),[Inf  3],[1 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_fidelity -args ARGS6{1}

%% Define argument types for entry-point 'compute_Integ_Domain_regularization'.
ARGS6 = cell(1,1);
ARGS6{1} = cell(21,1);

ARGS6{1}{1} = struct;
ARGS6{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS6{1}{1} = coder.typeof(ARGS6{1}{1},[Inf  1],[1 0]);

ARGS6{1}{2} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{3} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{4} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{5} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS6{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);

ARGS6{1}{11} = coder.typeof(0,[Inf  4],[1 0]);

ARGS6{1}{12} = struct;
ARGS6{1}{12}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{12} = coder.typeof(ARGS6{1}{12},[Inf  1],[1 0]);

ARGS6{1}{13} = struct;
ARGS6{1}{13}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{13} = coder.typeof(ARGS6{1}{13},[Inf  1],[1 0]);

ARGS6{1}{14} = struct;
ARGS6{1}{14}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS6{1}{14} = coder.typeof(ARGS6{1}{14},[Inf  1],[1 0]);

ARGS6{1}{15} = coder.typeof(0);
ARGS6{1}{16} = coder.typeof(0);
ARGS6{1}{17} = coder.typeof(0);

ARGS6{1}{18} = coder.typeof(0,[Inf  1],[1 0]);
ARGS6{1}{19} = coder.typeof(0,[Inf  1],[1 0]);
ARGS6{1}{20} = coder.typeof(0,[Inf  1],[1 0]);

ARGS6{1}{21} = coder.typeof(single(0),[Inf  3],[1 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_regularization -args ARGS6{1}

%% Define argument types for entry-point 'BsplineCompose3D'.
ARGS8 = cell(1,1);
ARGS8{1} = cell(9,1);

ARGS8{1}{1} = coder.typeof(0);
ARGS8{1}{2} = coder.typeof(0);
ARGS8{1}{3} = coder.typeof(0);

ARGS8{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{5} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{6} = coder.typeof(0,[1 Inf],[0 1]);

ARGS8{1}{7} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{8} = coder.typeof(0,[1 Inf],[0 1]);
ARGS8{1}{9} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg BsplineCompose3D -args ARGS8{1}

%% Define argument types for entry-point 'storePixelPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(int64(0));
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[Inf  3],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{6} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{7}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{7}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{7}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{7}.node = coder.typeof(0);
ARGS{1}{7}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nodes = coder.typeof(0,[Inf 1],[1 0]);
ARGS{1}{7} = coder.typeof(ARGS{1}{7},[Inf  1],[1 0]);
ARGS{1}{8} = struct;
ARGS{1}{8}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{8} = coder.typeof(ARGS{1}{8},[Inf  1],[1 0]);
ARGS{1}{9} = struct;
ARGS{1}{9}.pU = coder.typeof(0);
ARGS{1}{9}.pV = coder.typeof(0);
ARGS{1}{9}.pW = coder.typeof(0);
ARGS{1}{9}.maxiteration = coder.typeof(0);
ARGS{1}{9}.maxlevel = coder.typeof(0);
ARGS{1}{9}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.orderGauss = coder.typeof(0);
ARGS{1}{9}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.rho = coder.typeof(0,[3 1]);
ARGS{1}{9}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.lambda_1 = coder.typeof(0);
ARGS{1}{9}.lambda_2 = coder.typeof(0);
ARGS{1}{9}.lambda_3 = coder.typeof(0);
ARGS{1}{9}.sigma= coder.typeof(0);
ARGS{1}{9}.sigma_phi= coder.typeof(0);
ARGS{1}{9}.epsilon= coder.typeof(0);
ARGS{1}{9} = coder.typeof(ARGS{1}{9});

%% Invoke MATLAB Coder.
codegen -config cfg storePixelPhi -args ARGS{1}

%% Define argument types for entry-point 'GaussPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(8,1);
ARGS{1}{1} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = struct;
ARGS{1}{2}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{2}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{2}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{2}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{2}.node = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{6} = struct;
ARGS{1}{6}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6},[Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.pU = coder.typeof(0);
ARGS{1}{7}.pV = coder.typeof(0);
ARGS{1}{7}.pW = coder.typeof(0);
ARGS{1}{7}.maxiteration = coder.typeof(0);
ARGS{1}{7}.maxlevel = coder.typeof(0);
ARGS{1}{7}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.orderGauss = coder.typeof(0);
ARGS{1}{7}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.rho = coder.typeof(0,[3 1]);
ARGS{1}{7}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.lambda_1 = coder.typeof(0);
ARGS{1}{7}.lambda_2 = coder.typeof(0);
ARGS{1}{7}.lambda_3 = coder.typeof(0);
ARGS{1}{7}.sigma= coder.typeof(0);
ARGS{1}{7}.sigma_phi= coder.typeof(0);
ARGS{1}{7}.epsilon= coder.typeof(0);
ARGS{1}{7} = coder.typeof(ARGS{1}{7});
ARGS{1}{8} = coder.typeof(0);
%% Invoke MATLAB Coder.
codegen -config cfg GaussPhi -args ARGS{1}

