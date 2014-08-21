function [simpleSimOpt] = setSimulation(varargin)
% Create a structure for simulating a Dymola program, one for
% simulating the original program
% 
% [simpleSimOp] = setSimulations([simpleSimOp,]'Name1',value1,'Name2',value2,...) 
% creates a structure in which the named parameters have the specified values
% any unspecified name gets a default value
% 
% Simulation Properities:
% =======================
% 
% * SimPar: the simulation parameters used by dymosim, should be given,
% otherwise program-defined simulation is provided. 
% >> simPar = [0,1,0,100,1.e-4,8]; 
% or
% >> simPar = [0,1,0.01,0,1.e-4,8];
%
% 
% * x0values: The start values of variables, if not given, default is
% loaded from dsin.txt (or .mat)
%
%
% * pvalues: The values of parameters for simulation, if not given default 
% is loaded from dsin.txt (or .mat) 
% 
% 
% * varName: a cell array of (many) variable names as in [s,names]=dymosim, 
% each element is a string/regular expression. It represents the desired 
% variables, whose results will be returned from a simulation. Regular 
% expressions are instantiated within the next simulation. Default is all 
% existing varnames. Setting $varName$ erases $varInd$. $varInd$ are 
% generated within the next simulation.% 
% 
% 
% * acitvePar: a cell array of the active parameters, which can be changed
% e.g. during an optimization. Each cell is a string/regular expression.
% regular expressions will be inistiated. If not set, default will be all 
% existing active parameters present in the code, loaded in $dplabels$ 
% in the order of declaration 
% Setting $activePar$ erases $activeParInd$, which will be generated within the
% next simulation. $activePar$ is a field only in $gradSimOpt$ but not in
% $simpleSimOpt$.
%
%
% * newSim: If a new simulation should be performed true or false, eg.
% for another parameter set. default is true. In case false, $AllRes$ 
% and $AllLab$ should be not empty. [AllRes,AllLab] are the 
% results of an older simulation dymosim()
%
%
% * saveAllRes: Whether the hole simulation resuls [sim,lab]=dymosim()
% should be stored in $AllRes$ and $simLabels$ default is 
% 'false
%
% 
% * SimDir: the directory where dymosim.exe is. Default is current
% directory for $simpleSimOpt$ and $SimDir$/gradient for $gradSimOpt$
%
% 
% 
% Properities automatically set after setSimulation(...):
% =======================================================
% 
% 
% * plabels: The names of the parameters 
%
%
% * x0labels: The name of start variables 
%
% 
%
% The following are properties only in $gradSimOpt$
% * activeParInd: the indeces of required active parameters within
% dplabels, length(activeParInd) = length(activePar), otherwise error.
% should not be set by the user. It is generated automatically if
% $activePar$ ist set. Otherwise, only $dpind$ is used. 
%
% 
%
%
% Properties automatically within the next simulation:
% ====================================================
%
%
% * succeed: Whetehr the simulation terminated successifully.
%
%
% * varInd: a cell array of the indeces of the results of required
% variables in s where [s]=dymosim(SimPar), each cell represent the indeces 
% of a group of variables. default $varInd$ will be determined itself, if 
% $varName$ is set. the property $varInd$ has a higher priority than 
% $varName$. That is if both exists, $varInd$ will be directly used within
% the next simulation.
% Setting $varName$ erases $varInd$. $varName$ is generated within the next
% simulation.
%
% 
% * AllRes: Results of older simulation, needed in case other results
% for more refVar are needed, as in [s]  = dymosim(); $simResults$ will be
% stored in case $saveAllRes$ is on 
% 
%
% * AllLab: The labels of older simulations, as in [s,n] = dymosim();
% needed in case results for other refVar are sought. $AllLab$ will be
% stored in case $saveAllRes$ is on
%
% 
% * x: The results of the required variables corresponding to $varName$                                                                                                                
%
%
% Author : Atiyah Elsheikh, Austrian Institute of Technology GmbH
% Years  : 2014
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

%%
% Menu
%

    basicFields = {
      'SimPar'    ,   'vector of 6/7 values, see help dymosim'
      'x0values'  ,   'vector of length = length(x0), dymosim(simpar,x0,p)'
      'pvalues'   ,   'vector of length = length(p), dymosim(simpar,x0,p)'
      'varName'   ,   'cell array of regular expressions/strings'
      'activePar' ,   'cell array of regular expressions/strings' 
      'newSim'    ,   'true or false'
      'saveAllRes',   'true or false'
      'SimDir'    ,   'a path to the directory of the original model as a string'
    };
    if (nargin == 0) && (nargout == 0)
        fprintf('arguments : \n');  
        for i=1:size(basicFields,1) 
            fprintf('%s : %s \n',basicFields{i,1},basicFields{i,2});
        end 
        fprintf('\n');  
        return;
    end;

%%
% Fields of the simulation structure 
% 
    % basic simulation fields  
    Names = lower(basicFields(:,1)); 
   
    % essential active parameters (to be set e.g. via an optimization 
    parNames = lower({ 'activePar'});
        
    % fields to be computed / instantiated 
    outNames = lower({ 
                'varInd'
                'x0labels' 
                'plabels'
                'AllRes'
                'AllLab'
                'x'
                'time'
                'succeed'
            });
    
    % fields to be computed related to active parameters 
    outParNames = lower({  
                    'activeParInd'
                });
        
%%
%  initialization of propoereties and simulation structure 
% 

    firstArgIsStruct = nargin > 0 && isstruct(varargin{1}); 

    if(firstArgIsStruct)
        simpleSimOpt = varargin{1}; 
    end; 
    
    for i=1:length(Names)
        control.(Names{i}) = 0; %not set yet 
    end;
    for i=1:length(parNames)
        control.(parNames{i}) = 0; %not set yet
    end;
    
    if(~firstArgIsStruct)
        for i=1:length(Names)
            simpleSimOpt.(Names{i}) = [];
        end;
    
        for i=1:length(parNames)
            simpleSimOpt.(parNames{i}) = [];
        end;
        
        for i=1:length(outNames)
            simpleSimOpt.(outNames{i}) = [];
        end; 
        
         for i=1:length(outParNames)
             simpleSimOpt.(outParNames{i}) = [];
        end;     
        
    end
    clear parNames;
    clear Names; 
    
%%
% setting properties 
% 

    argind = firstArgIsStruct + 1;
   
    while argind <= nargin 
        arg = lower(varargin{argind}); 
        control.(arg) = 1;  
        % Setting Properties
        if(isfield(simpleSimOpt,arg))
            if(strcmp(arg,'varname')) % adding more variables 
                len = length(simpleSimOpt.varname);
                for k=1:length(varargin{argind+1}) 
                    simpleSimOpt.varname{len+k}{1} = ...
                                    varargin{argind+1}{k};
                end
             elseif(strcmp(arg,'activepar'))
                len = length(simpleSimOpt.activepar);
                for k=1:length(varargin{argind+1}) 
                    simpleSimOpt.activepar{len+k} = ...
                       varargin{argind+1}{k};
                end
            else        
                simpleSimOpt.(arg) = varargin{argind+1};
            end;
        else 
            error('setSimulation: argument %d is expected to be a property name\n',i); 
        end;
        argind = argind+2; 
    end; % while
    
    clear argind; 
       
    if(control.activepar) 
        simpleSimOpt.activeparind = [];
    end;     

%%
% default values for common property 
% 

    if(~firstArgIsStruct)
        if(~control.newsim) 
            simpleSimOpt.newsim = true;    
        end
        
        if(isempty(simpleSimOpt.simpar))
            simpleSimOpt.simpar = [0,1,0.01,0,1.e-4,8];
        end;   
        
        if(isempty(simpleSimOpt.saveallres)) 
            simpleSimOpt.saveallres = false; 
        end;    
        
        if(isempty(simpleSimOpt.simdir)) 
            simpleSimOpt.simdir = pwd;
        end;   
    end;
    
%%
%  Instantiating parameter names and their indices 
% 

    if(isfield(simpleSimOpt,'activepar') ... 
            && ~isempty(simpleSimOpt.activepar) ... 
            && isempty(simpleSimOpt.activeparind))
        [simpleSimOpt.activepar,simpleSimOpt.activeparind] = ...
            gradActInd(simpleSimOpt.activepar,simpleSimOpt.plabels);
    end

    
%%
% Output Propererties of simpleSimOpt
% 
    if(~firstArgIsStruct)
        currentDir = pwd;         
        cd(simpleSimOpt.simdir);
    
        % check if dsin.txt exists ...
        if(~exist('dsin.txt','file'))
            error(strcat('dsin.txt does not exist in ',  simpleSimOpt.simdir)); 
        end 
    
        [p1,x01,pnames1,x0names1]=loaddsin('dsin.txt');
        cd(currentDir);        
    
        if(isempty(simpleSimOpt.pvalues)) 
            simpleSimOpt.pvalues = p1; 
        end;    
    
        if(isempty(simpleSimOpt.x0values)) 
            simpleSimOpt.x0values = x01; 
        end;
    
        if(isempty(simpleSimOpt.plabels))
            simpleSimOpt.plabels = vecStr2cellArr(pnames1); 
        end;
    
        if(isempty(simpleSimOpt.x0labels))
            simpleSimOpt.x0labels = vecStr2cellArr(x0names1);
        end;
        
         clear pnames1 x0labels1 p1 x01 currentDir;
    end
    
%%    
% instantiate activePar 
% generate activeParInd
% 

    if(~isempty(simpleSimOpt.activepar) && isempty(simpleSimOpt.activeparind))
        [simpleSimOpt.activepar,gradSimOpt.activeparind] = ...
           gradActInd(simpleSimOpt.activepar,simpleSimOpt.plabels);
    end
    
%%
% Finalize 
%
   checkFieldsConsistency(simpleSimOpt,control,outNames,outParNames); 
   simpleSimOpt = orderfields(simpleSimOpt);
   clear outNames outParNames control;
   

%%
% function that matches a set of regexp $reg$ to a set of labels $labels$ 
% instantiate reg and generate the indeces of $reg{i}$ within $labels$
% 
    
    function [ret , ind] = gradActInd(reg,labels)
          
        
        cnt = 1; 
        m = length(reg);
        n = size(labels,1);
        flag = zeros(n,1);             
        for i=1:m
            r = reg{i};
            matchfound = 0;
            for j=1:n
                if(flag(j)) 
                    continue;
                elseif(regexp(labels{j},r) == 1)
                    flag(j) = 1;
                    matchfound = 1; 
                end;
            end;
            if(~matchfound)
                warning(strcat('regular expression : ',r, ' not matchable'));
            end;
        end;
               
        N = sum(flag);
        ind = zeros(N,1);
        ret = cell(N,1);
        for i=1:n
            if flag(i)
                ret{cnt,1} = labels{i};
                ind(cnt,1)= i;
                cnt = cnt + 1;
            end;
        end;
        
        clear cnt m n flag r N 
    
%%
% whether user properties make sense
%
    function checkFieldsConsistency(simpleSimOpt,control,outNames,outParNames) 
    
        for i=1:length(outNames) 
            if(isfield(control,outNames{i})) 
                if(control.(outNames{i}))
                    error(outNames{i},' should not be set by the user'); 
                end; 
            end;
        end; 

        for i=1:length(outParNames) 
            if(isfield(control,outParNames{i}))
                if (control.(outParNames{i}))
                    error(outParNames{i},' should not be set by the user'); 
                end;
            end;
        end; 
        
        if(~simpleSimOpt.newsim)
            if(isempty(simpleSimOpt.allres) || isempty(simpleSimOpt.alllab))
                error('setSimulation: in case newSim false, simulation results and labels should be set in simResuls and simLabels\n');
            end;
        end;
    
%%
% the indices of a set of strings (activePar) in a list of strings (parList) 
%
    function [ind] = indexIn(parList,activePar) 
    
    n = size(parList,1);
    m = size(activePar,1);
    cnt = 0; 
    ind = zeros(m,1);
    for j=1:m
        for i=cnt+1:n
            if strcmp(parList{i},activePar{j}) 
                ind(j) = i;
                cnt = i;
                break;
            end;
        end;
    end;
   
    clear n m cnt;
 
    
%%
% tramsfom a vector of strings to a cell array without leading blanks 
% 
    function x = vecStr2cellArr(x) 
        N = size(x,1);
        t = cell(N,1); 
        for i = 1:N 
            t{i,1} = deblank(x(i,:));
        end;
        x = t;
        clear t N;     
