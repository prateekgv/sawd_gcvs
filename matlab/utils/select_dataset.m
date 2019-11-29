%==========================================================================
% function [ h5names, ss ] = select_dataset( path )
%==========================================================================
% @descirption: Selects the dataset from which we extract the inertial
%               sensor measurements.
% @author     : Prateek Gundannavar
% @date       : 06/07/18
%
% @input      
%             - path        location of the APDM .h5 files
% @output
%             - h5names     A cell that contains all the .h5 files
%             - ss          Input selected
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================
function [ h5names, ss ] = select_dataset( path )

temp = dir(path);
h5names= {}; k = 1;
disp('List of datasets available');
for i=3:length(temp)
    fname = temp(i).name;
    if ~isempty(strfind(fname,'.h5'))
        h5names{k} = fname;
        fprintf('%d \t %s \n', k, h5names{k});
        k = k + 1;
    end
end

prompt = 'Select dataset: ';
ss = input(prompt);
end

