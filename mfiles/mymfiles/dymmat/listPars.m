%
% list parameters of a simulation structure created by @setSimulation  
% call: 
%     listPars(simopt[,ind])
%
%  simopt : A simulation structure set by the call setSimulation 
%  ind: optional, index of the parameters to be listed, e.g.
%       simopt.activeparind
%
%
%  Author : Atiyah Elsheikh, Austrian Institute of Technology GmbH 
%  Years  : 2011-2014
%


%
% Copyright (C) Atiyah Elsheikh 
% Atiyah.Elsheikh@ait.ac.at,a.m.g.Elsheikh@gmail.com) 2014, 
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

function listPars(simopt,ind)

    npars = length(simopt.plabels); 
    
    if(nargin < 2)
        for i=1:npars 
            fprintf('%5d  %50s  = %f\n',i,simopt.plabels{i},simopt.pvalues(i)); 
        end
    else
        for i=1:npars 
            if(sum(i == ind) > 0)
                fprintf('%5d  %50s  = %f\n',i,simopt.plabels{i},simopt.pvalues(i));
            end
        end
    end 
    
end 
