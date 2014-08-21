% Call:
%   function simopt = loadSimulation(simopt)
%
% Describtion: 
%   Rund and Load the simulation specified by simopt
%
% Input:
% - simopt: Simulation Handle created by @setSimulation
% - simroutine : the name of simulation routine
% 
% OutPut: 
% - simopt: Results are loaded within the simulation Handle
%    - simroutine : @dymosim by default
%
% Persistent:
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
% function simopt = loadSimulation(simopt)
% 
function simopt = loadSimulation(simopt,simroutine)


    if(nargin < 2) 
        simroutine = @dymosim
    end

%%
% simulate
%

if(simopt.newsim)
    currentdir = pwd; 
    cd(simopt.simdir);

    if(isempty(simopt.alllab) && ( ...
        isempty(simopt.varind) || ...
            (isfield(simopt,'dpind') && isempty(simopt.dpind)) ) )  
        [simopt.allres,simopt.alllab] = ...
            simroutine(simopt.simpar,simopt.x0values,simopt.pvalues,1);
        simopt.alllab = vecStr2cellArr(simopt.alllab);
    else 
        simopt.allres = ...
            simroutine(simopt.simpar,simopt.x0values,simopt.pvalues,1);
    end;     
    
    cd(currentdir);
    clear currentdir;
else 
   
    % check whetehr simopt.allres & alllab is here 
    if(~simpleSimOpt.newsim)
        if(isempty(simpleSimOpt.allres) || isempty(simpleSimOpt.alllab))
            error('setSimulation: in case newSim false, simulation results and labels should be set in AllRes and AllLab\n');
        end;
    end;
      
end; 

%%
%see simopt.varind 
% if empty then instantiate simopt.varname
%
    
for i = 1:length(simopt.varname)
    if(length(simopt.varind) < i || isempty(simopt.varind{i}))
        [t1,t2] = ...
            instantiateReg(simopt.varname{i},simopt.alllab);
        simopt.varname{i} = t1;
        simopt.varind{i} = t2;
        clear t1 t2;
    end;
end;

%%
% if varnames were not regular 
% 

len_ = length(simopt.varname); 
flag_ = 1; 
for i = 1:len_
    if(length(simopt.varname{i}) > 1) 
        flag_ = 0;
        break; 
    end
end

if(flag_) 
    %varind_{1} = []; 
    for i=1:len_
        varname_{1}(i,1) = simopt.varname{i}(1); 
        varind_{1}(i,1) = simopt.varind{i}; 
    end
    clear simopt.varname simopt.varind; 
    simopt.varname = varname_; 
    simopt.varind  = varind_; 
end

clear len_ varname_ varind_ flag_; 


%% 
% load results; 
% 
    simopt.time = simopt.allres(:,1);
    for i=1:length(simopt.varind)
        simopt.x{i} = simopt.allres(:,simopt.varind{i});    
    end;


%%
% finalize
% 
if(~simopt.saveallres)
        simopt.allres = [];
        simopt.alllab = [];
end;

if(simopt.time(end) + 1e-9 < simopt.simpar(2))
    simopt.succeed = 0;
else 
    simopt.succeed = 1;
end; 


%%
% 
% 

function x = vecStr2cellArr(x) 

    N = size(x,1);
    t = cell(N,1); 
    for i = 1:N 
        t{i,1} = deblank(x(i,:));
    end;
    x = t;
    clear t N;


%%
% function that matches a set of regexp $reg$ to a set of labels $labels$ 
% instantiate reg and generate the indeces of $reg{i}$ within $labels$
% $ind$ is the indeces of variables within labels, that matches $reg$ 
% 
    
function [reg , ind] = instantiateReg(reg,labels)
         
    cnt = 1; 
    m = size(reg,1);
    n = size(labels,1);
    flag = zeros(n,1);             
    for i=1:m
        r = reg{i};
        isreg = ~(isletter(r) | (r >= '0' & r <= '9') | r == '_'); 
        if(any(isreg))
            for j=1:n
                %labels{j,1}
                %regexp(labels{j,1},r)
                if(flag(j)) 
                    continue;
                %elseif(regexp(labels{j,1},r))    
                elseif(regexp(labels{j,1},r)==1)
                    flag(j) = 1;
                end;
            end;
        else
            for j=1:n
                if(flag(j)) 
                    continue; 
                elseif(strcmp(labels{j,1},r))
                    flag(j) = 1; 
                end;
            end;
        end;
    end;
               
    N = sum(flag);
    ind = zeros(N,1);
    reg = cell(N,1);
    for i=1:n
        if flag(i)
            reg{cnt,1} = labels{i,1};
            ind(cnt,1)= i;
            cnt = cnt + 1;
        end;
    end;
    
    if(cnt <= N)
        display(reg)
        display(ind)
        error('not all variables were found');
    end;
    
    [ind,I] = sort(ind); 
    reg = reg(I);
    %sort reg accoarding to ind
    
    clear cnt m n flag r N I 
        
%%
% find the indeces of the g_$desVar$[1] in $varList$
% NAP is the number of active Parameters
% 

function [ind] = gradInd(varList,desVar,NAP,prefix) 
    
    n = size(desVar,1);
    gVar = cell(n,1);

    for i=1:n 
       desVarW{i} = strrep(desVar{i},'.','_');  
    end
    
    for i=1:n
        gVar{i,1} = strcat(strcat(prefix,desVar{i}),'[1]');
        gVarW{i,1} = strcat(strcat(prefix,desVarW{i}),'[1]');
    end;
    
    m = size(varList,1);
  %  cnt = 1; 
    ind = zeros(n,1);
    for i=1:n
        for j=1:m
            if (strcmp(gVar{i},varList{j}) || strcmp(gVarW{i},varList{j})) 
                ind(i) = j;
   %             cnt = cnt + NAP;
                break;
            end;
    %        cnt = cnt + 1; 
        end;
        
        if(size(ind,1) < i) 
           fprintf('Warning a derivative for %s does not exist\n',gVar{i});
           ind(i) = -1;
        end;  
    end;
   
    clear n m cnt;
    
    
%%
% 
        
