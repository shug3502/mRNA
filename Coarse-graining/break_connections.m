function Q = break_connections(Q,index1,index2)
%pick a cell to block a given RC both IN and OUT
Q(index1,index2)=0;
Q(index2,index1)=0;
r = 1:16; r(index1)=[];
Q(index1,index1)=-sum(Q(index1,r));
r = 1:16; r(index2)=[];
Q(index2,index2)=-sum(Q(index2,r));
