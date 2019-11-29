%==========================================================================
% function [ start_index, stop_index ] = edge_detection( fog )
%==========================================================================
% @descirption: Detect the edges of the detector
% @author     : Prateek Gundannavar
% @date       : 03/11/19
%
% @input      
%             - fog       array containing boolean enteries
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================
function [ start_index, stop_index ] = edge_detection( fog )

if (fog(1,1)==1) && (fog(1,end)==1)
    ind=[1 abs(diff(fog)) 1 >0]>0;
elseif (fog(1,1)==1) && ~(fog(1,end)==1)
    ind=[1 abs(diff(fog)) >0]>0;
elseif ~(fog(1,1)==1) && (fog(1,end)==1)
    ind=[abs(diff(fog)) 1 >0]>0;
elseif ~(fog(1,1)==1) && ~(fog(1,end)==1)
    ind=[abs(diff(fog)) >0]>0;
end

ind=find(ind>0);
n=length(ind);
if (n~=0 && mod(n,2)==0)
    M = reshape(ind,2,(length(ind)/2))';
    start_index=M(:,1);
    stop_index=M(:,2);
elseif (n~=0 && mod(n,2)==1)
    ind=ind(1:end-1);
    M = reshape(ind,2,(length(ind)/2))';
    start_index=M(:,1);
    stop_index=M(:,2);
else
    start_index = [];
    stop_index = [];
end

end

