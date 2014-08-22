%
% This directory contains the following code:
% 
% Routines for setting and loading dymola simulations within 
%   setSimulation.m   : setting a dymola simulation record 
%   loadSimulation.m  : simulating a dymola simulation record 
% 
% Help routines 
%   listPars.m      : listing model parameters of a simulation structure
%   listParsReg.m   :
%
%   getSimPar.m     : geting the value of a named parameter
%   setSimPar.m     : setting the value of a named parameter
%   setSimParReg.m  : setting the value of parameters matching a regular
%                     expression
%   setSimParRegUniform.m  : setting the value of parameters to random value
%                            according to a uniform distribution, the 
%                            parameters can be specified with reg. express.
%  
%
%  Author : Atiyah Elsheikh, Austrian Institute of Technology GmbH 
%  Years  : 2011-2014
%


%
% Copyright (C) Atiyah Elsheikh 
% Atiyah.Elsheikh@ait.ac.at,a.m.g.Elsheikh@gmail.com) 2011-2014, 
% AIT Austrian Institute of Technology GmbH
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
% along with dymmat. If not, see <http://www.gnu.org/licenses/>
%