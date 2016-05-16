function nc_transitions = get_nc_transitions(can_exit_oocyte)

%created 10/5/16
%last edit 10/5/16
%based on figure 1 a in caceres 2005 development
%can_exit_oocyte is a logical input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
    can_exit_oocyte=0;
end

connections{1} = [2,3,5,9];
connections{2} = [1,4,6,10];
connections{3} = [1,7,11];
connections{4} = [2,8,12];
connections{5} = [1,13];
connections{6} = [2,14];
connections{7} = [3,15];
connections{8} = [4,16];
connections{9} = 1;
connections{10} = 2;
connections{11} = 3;
connections{12} = 4;
connections{13} = 5;
connections{14} = 6;
connections{15} = 7;
connections{16} = 8;

nc_transitions = zeros(16);
for j=1:16
    nc_transitions(j,connections{j}) = 1; %transition rates
    nc_transitions(j,j) = -numel(connections{j}); %rate leaving
end

if ~can_exit_oocyte
    nc_transitions(1,:)=0; %no transitions to leave oocyte (cell 1)
end