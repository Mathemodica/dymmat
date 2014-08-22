%
% Set a specific model parameter to a specific value 
% call: 
%     simopt = setSimPar(simopt,pname,pval)
%
%  simopt   : A simulation structure set by the call setSimulation
%  pname    : parameter name 
%  pval     : parameter value 
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

function simopt = setSimPar(simopt,pname,pval)

    npars = length(simopt.plabels); 
    
    for i=1:npars 
        if(strcmp(simopt.plabels{i},pname))
            simopt.pvalues(i) = pval;
            return; 
        end
    end

    error('input parameter name not found'); 
end