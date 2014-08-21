%
% A script for adding 
%
% Author : Atiyah Elsheikh, Austrian Institute of Technology GmbH © 2014 
% 


%
% Copyright (C) Atiyah Elsheikh 
% (Atiyah.Elsheikh@ait.ac.at,a.m.g.Elsheikh@gmail.com) 2014, 
%  AIT Austrian Institute of Technology GmbH
% This file is part of the software dymmat toolbox
% dymmat is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% dymmat is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% You should have received a copy of the GNU Lesser General Public License
% along with Foobar. If not, see <http://www.gnu.org/licenses/>
%

dymhome = 'C:\Program Files (x86)\Dymola 2015\'; 
currdir = pwd; 
base = strcat(dymhome,'\mfiles')
traj = strcat(dymhome,'\mfiles\traj')
dymtools = strcat(dymhome,'\mfiles\dymtools')

% Advanced Matlab Interface to Dymola 
dymmat = strcat(currdir,'\mymfiles\dymmat')

addpath(base,traj,dymtools,dymmat);
