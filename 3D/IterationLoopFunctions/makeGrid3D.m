function [ gridX, gridY, gridZ, pointsX, pointsY,pointsZ ] = makeGrid3D(ngx, ngy, ngz, nx, ny ,nz )

%makes a 3D grid overlaying the image

%--Input Variable:
%ngx: number of grid points in X-direction
%ngy: number of grid points in Y-direction
%ngz: number of grid points in Z-direction
%nx: image size in X-direction
%ny: image size in Y-direction
%nz: image size in Z-direction

%--Output Variable:

%gridX: grid x-vertices
%gridY: grid y-vertices
%gridZ: grid z-vertices
%pointsX: pixel grid points in X-direction
%pointsY: pixel grid points in Y-direction
%pointsZ: pixel grid points in Z-direction

pointsX = round(linspace(1,nx,ngx));
pointsY = round(linspace(1,ny,ngy));
pointsZ = round(linspace(1,nz,ngz));

[gridX, gridY, gridZ] = meshgrid(pointsX, pointsY, pointsZ);

end