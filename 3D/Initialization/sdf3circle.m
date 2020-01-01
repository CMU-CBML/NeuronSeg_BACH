function f = sdf3circle(nrow,ncol,nslice, ic,jc,kc, r)
%   sdf2circle(nrow,ncol, ic,jc,r) computes the signed distance to a circle
%   input: 
%       nrow: number of rows
%       ncol: number of columns
%       (ic,jc): center of the circle
%       r: radius of the circle
%   output: 
%       f: signed distance to the circle
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li


[X,Y,Z] = meshgrid(1:ncol, 1:nrow, 1:nslice);

f = sqrt((X-jc).^2 + (Y-ic).^2 + (Z-kc).^2)-r;

