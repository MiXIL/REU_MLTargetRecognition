function [R,dt] = ray_range_2medium(p_ant,p_scat,h_gnd,er1,er2)
%ray_range_2medium - compute the free space range bin of a point reflector due to travel through 2 media. 
% Useful for calculating apparent free space range of subsurface targets
% for a UAV GPR.
%
% Coordinates:
%                 ^ z
%                /
%               /
%       tx/rx   ----> x
%              |
%              |
%              v
%              y
%
% -------------------------- ground
%
%                   o target
%
% Syntax:  R = ray_range_2medium(ptx,pscat,h_gnd,er1,er2)
%
% Inputs:
%    p_ant - [x,y,z] coordinate of tx/rx antenna (uav)
%    p_scat - [x,y,z] coordinate of target scatterer
%    h_gnd - height of antenna/ UAV above surface
%    er1 - relative permittivity of upper media (nominally air -> er1 = 1)
%    er2 - relative permittivity of lower media (nominally ground)
%
% Outputs:
%    R - free space range bin
%    dt - time delay (dt = 2*R/c)
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: refraction_path.m

% Author: Samuel Prager
% Microwave Systems, Sensors and Imaging Lab (MiXiL) 
% University of Southern California
% Email: sprager@usc.edu
% Created: 2018/03/28 10:47:39; Last Revised: 2018/03/28 10:47:39
%
% Copyright 2012-2018 University of Southern California
%------------- BEGIN CODE --------------

% coordinates of target relative to antenna
p_rel = p_scat-p_ant;

% length of x-z plane triangle base
rloc = sqrt(p_rel(1)^2+p_rel(3)^2);

R = sqrt(er1)*sqrt(h_gnd^2 + rloc^2*(h_gnd/p_rel(2))^2)+sqrt(er2)*sqrt((p_rel(2)-h_gnd)^2+rloc^2*((p_rel(2)-h_gnd)/p_rel(2))^2);
dt = 2*R/3e8;
%------------- END OF CODE --------------
