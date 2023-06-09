function mol=CARTcell_ccm_C(free_net,free_xch,in)

kernel_net=[ 1 -0 -0 -0 -0 -0 -0 -3 -1 1 1 -3 -0 -0 1 -1 -0 -1 -1 -1;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 0.5 0.5 0.5 0.5 1 -0 0.5 1 0.5 -0.5 0.5 -0 0.5 -0.5 0.5 -0 -0.5 0.5 0.5;
 -0 0.5 0.5 0.5 0.5 1 -0 0.5 1 0.5 -0.5 0.5 -0 0.5 -0.5 0.5 -0 -0.5 0.5 0.5;
 -0 0.5 0.5 0.5 0.5 -0 -0 -2.5 -0 0.5 -0.5 0.5 -0 0.5 -0.5 0.5 -0 -0.5 0.5 0.5;
 -0 0.5 0.5 0.5 0.5 -0 -0 -0.5 -0 0.5 -0.5 0.5 -0 0.5 -0.5 0.5 -0 -0.5 0.5 0.5;
 -0 0.5 0.5 0.5 0.5 -0 -0 -0.5 -0 0.5 -0.5 0.5 -0 0.5 -0.5 0.5 -0 -0.5 0.5 0.5;
 -0 -0.5 0.5 0.5 0.5 -0 -0 -0.5 -0 0.5 -0.5 0.5 -0 0.5 -0.5 0.5 -0 -0.5 0.5 0.5;
 -0 -0 1 1 1 -0 -0 -0 -0 1 -1 1 -0 1 -1 1 -0 -1 1 1;
 -0 -0 1 1 1 -0 -0 -0 -0 1 -1 1 -0 1 -1 1 -0 -1 1 1;
 -0 -0 -0 1 1 -0 -0 -0 -0 1 -1 1 -0 1 -1 1 -0 -1 1 1;
 -0 -0 -0 1 1 -0 1 -0 -0 -0 -1 -0 1 1 -0 -0 -0 -0 -0 1;
 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 1 -0 -0 -1 -0 -1 1 -0 1 -1 -0 1 -1 -0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -1 -0 -0 -0 -0 -0 -0 -0;
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 -0 3 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 3 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 1 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 2 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -1 1 -0 -0 -1 1 -0 -0 1 1;
 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -1 1 -0 -0 1 1;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -1 1 -0 -0 1 -0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0;
 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 1 -1 -0 -0 -1 -0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 1 -1 -0 -0 -1 -0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -1 1 -0 -0 -0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];

kernel_xch=[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

net=sparse(kernel_net)*free_net;
xch=sparse(kernel_xch)*free_xch;

f=xch+max(0,net);
b=xch-min(0,net);

%%
A1=zeros(126);
A1(1,1)=-b(1)-b(19)-b(20)-b(23)-b(31)-b(35)-b(36)-f(2)-f(21);
A1(1,19)=f(19)+b(21);
A1(1,20)=f(20);
A1(1,21)=f(23);
A1(1,22)=f(31);
A1(1,23)=f(35);
A1(1,24)=f(36);
A1(2,2)=-b(5)-b(26)-b(28)-f(6);
A1(2,4)=f(28);
A1(2,8)=f(26);
A1(2,25)=f(5);
A1(2,26)=b(6);
A1(3,3)=-b(25)-f(26)-f(27);
A1(3,6)=b(26);
A1(3,7)=b(27);
A1(3,27)=f(25);
A1(4,2)=b(28);
A1(4,4)=-b(27)-f(28);
A1(4,9)=f(27);
A1(5,5)=-b(7)-f(8)-f(14);
A1(5,26)=f(7);
A1(5,28)=b(8);
A1(6,3)=f(26);
A1(6,6)=-b(5)-b(26)-b(28)-f(6);
A1(6,7)=f(28);
A1(6,29)=f(5);
A1(6,30)=b(6);
A1(7,3)=f(27);
A1(7,6)=b(28);
A1(7,7)=-b(27)-f(28);
A1(8,2)=b(26);
A1(8,8)=-b(28)-f(26);
A1(8,31)=f(28);
A1(9,4)=b(27);
A1(9,9)=-b(24)-f(27)-f(29);
A1(9,13)=f(24);
A1(10,10)=-b(31)-b(32)-f(33);
A1(10,14)=f(31);
A1(10,32)=b(33);
A1(11,11)=-b(34)-f(50);
A1(11,32)=f(34);
A1(12,12)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A1(12,33)=f(21);
A1(12,34)=f(34)+b(33);
A1(12,35)=f(40);
A1(12,36)=b(19);
A1(13,9)=b(24);
A1(13,13)=-b(23)-f(24)-f(25);
A1(13,37)=f(23);
A1(13,38)=b(25);
A1(14,10)=b(31);
A1(14,14)=-b(12)-b(20)-f(13)-f(16)-f(21)-f(31)-f(41);
A1(14,16)=f(12);
A1(14,39)=f(20);
A1(14,40)=b(13);
A1(14,41)=b(21);
A1(14,42)=b(41);
A1(15,15)=-b(31)-b(32)-f(33);
A1(15,33)=f(31);
A1(15,43)=b(33);
A1(16,14)=b(12);
A1(16,16)=-b(11)-b(19)-f(12);
A1(16,41)=f(19);
A1(16,44)=f(11);
A1(17,17)=-b(35)-b(41)-f(36)-f(43);
A1(17,34)=f(35);
A1(17,45)=f(41)+b(43);
A1(17,46)=b(36);
A1(18,18)=-b(36)-f(37);
A1(18,47)=f(36);
A1(18,48)=b(37);
A1(19,1)=f(21)+b(19);
A1(19,19)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A1(19,20)=f(40);
A1(19,49)=f(34)+b(33);
A1(20,1)=b(20);
A1(20,19)=b(40);
A1(20,20)=-b(39)-f(20)-f(40);
A1(20,50)=f(39);
A1(21,1)=b(23);
A1(21,21)=-b(22)-f(23);
A1(21,51)=f(22);
A1(22,1)=b(31);
A1(22,22)=-b(12)-b(20)-f(13)-f(16)-f(21)-f(31)-f(41);
A1(22,52)=f(12);
A1(22,53)=f(20);
A1(22,54)=b(13);
A1(22,55)=b(21);
A1(22,56)=b(41);
A1(23,1)=b(35);
A1(23,23)=-b(33)-f(34)-f(35);
A1(23,55)=f(33)+b(34);
A1(24,1)=b(36);
A1(24,24)=-b(35)-b(41)-f(36)-f(43);
A1(24,49)=f(35);
A1(24,57)=f(41)+b(43);
A1(25,2)=b(5);
A1(25,25)=-b(4)-f(5)-f(18)-f(22);
A1(25,58)=f(4);
A1(25,59)=b(22);
A1(26,2)=f(6);
A1(26,5)=b(7);
A1(26,26)=-b(6)-f(7);
A1(27,3)=b(25);
A1(27,27)=-b(23)-f(24)-f(25);
A1(27,59)=f(23);
A1(27,60)=b(24);
A1(28,5)=f(8);
A1(28,28)=-b(7)-b(8)-b(26)-b(27)-f(9)-f(28);
A1(28,61)=f(7);
A1(28,62)=f(26)+f(27);
A1(28,63)=b(9);
A1(28,64)=b(28);
A1(29,6)=b(5);
A1(29,29)=-b(4)-f(5)-f(18)-f(22);
A1(29,37)=b(22);
A1(29,65)=f(4);
A1(30,6)=f(6);
A1(30,30)=-b(6)-f(7);
A1(30,66)=b(7);
A1(31,8)=b(28);
A1(31,31)=-b(27)-f(28);
A1(31,60)=f(27);
A1(32,10)=f(33);
A1(32,11)=b(34);
A1(32,32)=-b(33)-f(34)-f(35);
A1(32,67)=b(35);
A1(33,12)=b(21);
A1(33,15)=b(31);
A1(33,33)=-b(12)-b(20)-f(13)-f(16)-f(21)-f(31)-f(41);
A1(33,35)=f(20);
A1(33,36)=f(12);
A1(33,68)=b(13);
A1(33,69)=b(41);
A1(34,12)=f(33)+b(34);
A1(34,17)=b(35);
A1(34,34)=-b(33)-f(34)-f(35);
A1(35,12)=b(40);
A1(35,33)=b(20);
A1(35,35)=-b(39)-f(20)-f(40);
A1(35,70)=f(39);
A1(36,12)=f(19);
A1(36,33)=b(12);
A1(36,36)=-b(11)-b(19)-f(12);
A1(36,71)=f(11);
A1(37,13)=b(23);
A1(37,29)=f(22);
A1(37,37)=-b(22)-f(23);
A1(38,13)=f(25);
A1(38,38)=-b(25)-f(26)-f(27);
A1(38,72)=b(26);
A1(38,73)=b(27);
A1(39,14)=b(20);
A1(39,39)=-b(39)-f(20)-f(40);
A1(39,41)=b(40);
A1(39,74)=f(39);
A1(40,14)=f(13);
A1(40,40)=-b(13)-f(17);
A1(41,14)=f(21);
A1(41,16)=b(19);
A1(41,39)=f(40);
A1(41,41)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A1(41,75)=f(34)+b(33);
A1(42,14)=f(41);
A1(42,42)=-b(41)-f(42);
A1(43,15)=f(33);
A1(43,43)=-b(33)-f(34)-f(35);
A1(43,76)=b(34);
A1(43,77)=b(35);
A1(44,16)=b(11);
A1(44,44)=-b(10)-f(11)-f(15);
A1(44,78)=f(10);
A1(45,17)=f(43)+b(41);
A1(45,45)=-b(43)-b(45)-f(41)-f(46)-f(49);
A1(45,79)=b(46);
A1(46,17)=f(36);
A1(46,46)=-b(36)-f(37);
A1(46,80)=b(37);
A1(47,18)=b(36);
A1(47,47)=-b(35)-b(41)-f(36)-f(43);
A1(47,75)=f(35);
A1(47,81)=f(41)+b(43);
A1(48,18)=f(37);
A1(48,48)=-b(37)-f(38);
A1(48,83)=b(38);
A1(49,19)=f(33)+b(34);
A1(49,24)=b(35);
A1(49,49)=-b(33)-f(34)-f(35);
A1(50,20)=b(39);
A1(50,50)=-b(38)-f(39);
A1(50,82)=f(38);
A1(51,21)=b(22);
A1(51,51)=-b(4)-f(5)-f(18)-f(22);
A1(51,72)=b(5);
A1(51,84)=f(4);
A1(52,22)=b(12);
A1(52,52)=-b(11)-b(19)-f(12);
A1(52,55)=f(19);
A1(52,85)=f(11);
A1(53,22)=b(20);
A1(53,53)=-b(39)-f(20)-f(40);
A1(53,55)=b(40);
A1(53,83)=f(39);
A1(54,22)=f(13);
A1(54,54)=-b(13)-f(17);
A1(55,22)=f(21);
A1(55,23)=f(34)+b(33);
A1(55,52)=b(19);
A1(55,53)=f(40);
A1(55,55)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A1(56,22)=f(41);
A1(56,56)=-b(41)-f(42);
A1(57,24)=f(43)+b(41);
A1(57,57)=-b(43)-b(45)-f(41)-f(46)-f(49);
A1(57,86)=b(46);
A1(58,25)=b(4);
A1(58,58)=-b(3)-f(4);
A1(59,25)=f(22);
A1(59,27)=b(23);
A1(59,59)=-b(22)-f(23);
A1(60,27)=f(24);
A1(60,31)=b(27);
A1(60,60)=-b(24)-f(27)-f(29);
A1(61,28)=b(7);
A1(61,61)=-b(6)-f(7);
A1(61,64)=f(6);
A1(62,28)=b(26)+b(27);
A1(62,62)=-b(25)-f(26)-f(27);
A1(62,87)=f(25);
A1(63,28)=f(9);
A1(63,63)=-b(9)-f(10);
A1(63,85)=b(10);
A1(64,28)=f(28);
A1(64,61)=b(6);
A1(64,64)=-b(5)-b(26)-b(28)-f(6);
A1(64,88)=f(5);
A1(64,89)=f(26);
A1(65,29)=b(4);
A1(65,65)=-b(3)-f(4);
A1(66,30)=f(7);
A1(66,66)=-b(7)-f(8)-f(14);
A1(66,90)=b(8);
A1(67,32)=f(35);
A1(67,67)=-b(35)-b(41)-f(36)-f(43);
A1(67,91)=f(41)+b(43);
A1(67,92)=b(36);
A1(68,33)=f(13);
A1(68,68)=-b(13)-f(17);
A1(69,33)=f(41);
A1(69,69)=-b(41)-f(42);
A1(70,35)=b(39);
A1(70,70)=-b(38)-f(39);
A1(70,80)=f(38);
A1(71,36)=b(11);
A1(71,71)=-b(10)-f(11)-f(15);
A1(71,93)=f(10);
A1(72,38)=f(26);
A1(72,51)=f(5);
A1(72,72)=-b(5)-b(26)-b(28)-f(6);
A1(72,73)=f(28);
A1(72,94)=b(6);
A1(73,38)=f(27);
A1(73,72)=b(28);
A1(73,73)=-b(27)-f(28);
A1(74,39)=b(39);
A1(74,74)=-b(38)-f(39);
A1(74,95)=f(38);
A1(75,41)=f(33)+b(34);
A1(75,47)=b(35);
A1(75,75)=-b(33)-f(34)-f(35);
A1(76,43)=f(34);
A1(76,76)=-b(34)-f(50);
A1(77,43)=f(35);
A1(77,77)=-b(35)-b(41)-f(36)-f(43);
A1(77,96)=f(41)+b(43);
A1(77,97)=b(36);
A1(78,44)=b(10);
A1(78,78)=-b(9)-f(10);
A1(78,98)=f(9);
A1(79,45)=f(46);
A1(79,79)=-b(44)-b(46)-f(47);
A1(80,46)=f(37);
A1(80,70)=b(38);
A1(80,80)=-b(37)-f(38);
A1(81,47)=f(43)+b(41);
A1(81,81)=-b(43)-b(45)-f(41)-f(46)-f(49);
A1(81,99)=b(46);
A1(82,50)=b(38);
A1(82,82)=-b(37)-f(38);
A1(82,97)=f(37);
A1(83,48)=f(38);
A1(83,53)=b(39);
A1(83,83)=-b(38)-f(39);
A1(84,51)=b(4);
A1(84,84)=-b(3)-f(4);
A1(85,52)=b(11);
A1(85,63)=f(10);
A1(85,85)=-b(10)-f(11)-f(15);
A1(86,57)=f(46);
A1(86,86)=-b(44)-b(46)-f(47);
A1(87,62)=b(25);
A1(87,87)=-b(23)-f(24)-f(25);
A1(87,100)=f(23);
A1(87,101)=b(24);
A1(88,64)=b(5);
A1(88,88)=-b(4)-f(5)-f(18)-f(22);
A1(88,100)=b(22);
A1(88,102)=f(4);
A1(89,64)=b(26);
A1(89,89)=-b(28)-f(26);
A1(89,103)=f(28);
A1(90,66)=f(8);
A1(90,90)=-b(7)-b(8)-b(26)-b(27)-f(9)-f(28);
A1(90,93)=b(9);
A1(90,104)=f(7);
A1(90,105)=f(26)+f(27);
A1(90,106)=b(28);
A1(91,67)=f(43)+b(41);
A1(91,91)=-b(43)-b(45)-f(41)-f(46)-f(49);
A1(91,107)=b(46);
A1(92,67)=f(36);
A1(92,92)=-b(36)-f(37);
A1(92,95)=b(37);
A1(93,71)=b(10);
A1(93,90)=f(9);
A1(93,93)=-b(9)-f(10);
A1(94,72)=f(6);
A1(94,94)=-b(6)-f(7);
A1(94,108)=b(7);
A1(95,74)=b(38);
A1(95,92)=f(37);
A1(95,95)=-b(37)-f(38);
A1(96,77)=f(43)+b(41);
A1(96,96)=-b(43)-b(45)-f(41)-f(46)-f(49);
A1(96,109)=b(46);
A1(97,77)=f(36);
A1(97,82)=b(37);
A1(97,97)=-b(36)-f(37);
A1(98,78)=b(9);
A1(98,98)=-b(7)-b(8)-b(26)-b(27)-f(9)-f(28);
A1(98,108)=f(8);
A1(98,110)=f(7);
A1(98,111)=f(26)+f(27);
A1(98,112)=b(28);
A1(99,81)=f(46);
A1(99,99)=-b(44)-b(46)-f(47);
A1(100,87)=b(23);
A1(100,88)=f(22);
A1(100,100)=-b(22)-f(23);
A1(101,87)=f(24);
A1(101,101)=-b(24)-f(27)-f(29);
A1(101,103)=b(27);
A1(102,88)=b(4);
A1(102,102)=-b(3)-f(4);
A1(103,89)=b(28);
A1(103,101)=f(27);
A1(103,103)=-b(27)-f(28);
A1(104,90)=b(7);
A1(104,104)=-b(6)-f(7);
A1(104,106)=f(6);
A1(105,90)=b(26)+b(27);
A1(105,105)=-b(25)-f(26)-f(27);
A1(105,113)=f(25);
A1(106,90)=f(28);
A1(106,104)=b(6);
A1(106,106)=-b(5)-b(26)-b(28)-f(6);
A1(106,114)=f(5);
A1(106,115)=f(26);
A1(107,91)=f(46);
A1(107,107)=-b(44)-b(46)-f(47);
A1(108,94)=f(7);
A1(108,98)=b(8);
A1(108,108)=-b(7)-f(8)-f(14);
A1(109,96)=f(46);
A1(109,109)=-b(44)-b(46)-f(47);
A1(110,98)=b(7);
A1(110,110)=-b(6)-f(7);
A1(110,112)=f(6);
A1(111,98)=b(26)+b(27);
A1(111,111)=-b(25)-f(26)-f(27);
A1(111,116)=f(25);
A1(112,98)=f(28);
A1(112,110)=b(6);
A1(112,112)=-b(5)-b(26)-b(28)-f(6);
A1(112,117)=f(5);
A1(112,118)=f(26);
A1(113,105)=b(25);
A1(113,113)=-b(23)-f(24)-f(25);
A1(113,119)=f(23);
A1(113,120)=b(24);
A1(114,106)=b(5);
A1(114,114)=-b(4)-f(5)-f(18)-f(22);
A1(114,119)=b(22);
A1(114,121)=f(4);
A1(115,106)=b(26);
A1(115,115)=-b(28)-f(26);
A1(115,122)=f(28);
A1(116,111)=b(25);
A1(116,116)=-b(23)-f(24)-f(25);
A1(116,123)=f(23);
A1(116,124)=b(24);
A1(117,112)=b(5);
A1(117,117)=-b(4)-f(5)-f(18)-f(22);
A1(117,123)=b(22);
A1(117,125)=f(4);
A1(118,112)=b(26);
A1(118,118)=-b(28)-f(26);
A1(118,126)=f(28);
A1(119,113)=b(23);
A1(119,114)=f(22);
A1(119,119)=-b(22)-f(23);
A1(120,113)=f(24);
A1(120,120)=-b(24)-f(27)-f(29);
A1(120,122)=b(27);
A1(121,114)=b(4);
A1(121,121)=-b(3)-f(4);
A1(122,115)=b(28);
A1(122,120)=f(27);
A1(122,122)=-b(27)-f(28);
A1(123,116)=b(23);
A1(123,117)=f(22);
A1(123,123)=-b(22)-f(23);
A1(124,116)=f(24);
A1(124,124)=-b(24)-f(27)-f(29);
A1(124,126)=b(27);
A1(125,117)=b(4);
A1(125,125)=-b(3)-f(4);
A1(126,118)=b(28);
A1(126,124)=f(27);
A1(126,126)=-b(27)-f(28);

B1=zeros(126,23);
B1(1,1)=-f(1);
B1(10,2)=-f(32);
B1(12,3)=-f(48);
B1(15,4)=-f(32);
B1(19,5)=-f(48);
B1(41,6)=-f(48);
B1(45,7)=-f(45);
B1(55,8)=-f(48);
B1(57,9)=-f(45);
B1(58,10)=-f(3);
B1(65,11)=-f(3);
B1(79,12)=-f(44);
B1(81,13)=-f(45);
B1(84,14)=-f(3);
B1(86,15)=-f(44);
B1(91,16)=-f(45);
B1(96,17)=-f(45);
B1(99,18)=-f(44);
B1(102,19)=-f(3);
B1(107,20)=-f(44);
B1(109,21)=-f(44);
B1(121,22)=-f(3);
B1(125,23)=-f(3);

y1=[ in.CO2_IN_1; in.AC_IN_2; in.OAA_IN_2; in.AC_IN_1; in.OAA_IN_4; in.OAA_IN_3; in.Glu_IN_3; in.OAA_IN_1; in.Glu_IN_1; in.GLC_IN_3; in.GLC_IN_2; in.Gln_IN_3; in.Glu_IN_2; in.GLC_IN_1; in.Gln_IN_1; in.Glu_IN_4; in.Glu_IN_5; in.Gln_IN_2; in.GLC_IN_4; in.Gln_IN_4; in.Gln_IN_5; in.GLC_IN_5; in.GLC_IN_6];
x1=sparse(A1)\(B1*sparse(y1));

%%
A2=zeros(86);
A2(1,1)=-b(31)-b(32)-f(33);
A2(1,9)=f(31);
A2(1,14)=b(33);
A2(2,2)=-b(25)-f(26)-f(27);
A2(2,4)=b(26);
A2(2,5)=b(27);
A2(2,13)=f(25);
A2(3,3)=-b(34)-f(50);
A2(3,14)=f(34);
A2(4,2)=f(26);
A2(4,4)=-b(5)-b(26)-b(28)-f(6);
A2(4,5)=f(28);
A2(4,15)=f(5);
A2(4,16)=b(6);
A2(5,2)=f(27);
A2(5,4)=b(28);
A2(5,5)=-b(27)-f(28);
A2(6,6)=-b(27)-f(28);
A2(6,17)=b(28);
A2(7,7)=-b(7)-f(8)-f(14);
A2(7,18)=f(7);
A2(7,19)=b(8);
A2(8,8)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A2(8,9)=f(21);
A2(8,10)=b(19);
A2(8,20)=f(34)+b(33);
A2(8,21)=f(40);
A2(9,1)=b(31);
A2(9,8)=b(21);
A2(9,9)=-b(12)-b(20)-f(13)-f(16)-f(21)-f(31)-f(41);
A2(9,10)=f(12);
A2(9,21)=f(20);
A2(9,22)=b(13);
A2(9,23)=b(41);
A2(10,8)=f(19);
A2(10,9)=b(12);
A2(10,10)=-b(11)-b(19)-f(12);
A2(10,24)=f(11);
A2(11,11)=-b(36)-f(37);
A2(11,12)=f(36);
A2(11,25)=b(37);
A2(12,11)=b(36);
A2(12,12)=-b(35)-b(41)-f(36)-f(43);
A2(12,20)=f(35);
A2(12,26)=f(41)+b(43);
A2(13,2)=b(25);
A2(13,13)=-b(23)-f(24)-f(25);
A2(13,27)=f(23);
A2(13,28)=b(24);
A2(14,1)=f(33);
A2(14,3)=b(34);
A2(14,14)=-b(33)-f(34)-f(35);
A2(14,29)=b(35);
A2(15,4)=b(5);
A2(15,15)=-b(4)-f(5)-f(18)-f(22);
A2(15,30)=f(4);
A2(15,31)=b(22);
A2(16,4)=f(6);
A2(16,16)=-b(6)-f(7);
A2(16,32)=b(7);
A2(17,6)=f(28);
A2(17,17)=-b(5)-b(26)-b(28)-f(6);
A2(17,18)=b(6);
A2(17,33)=f(5);
A2(18,7)=b(7);
A2(18,17)=f(6);
A2(18,18)=-b(6)-f(7);
A2(19,7)=f(8);
A2(19,19)=-b(7)-b(8)-b(26)-b(27)-f(9)-f(28);
A2(19,34)=f(7);
A2(19,35)=f(26)+f(27);
A2(19,36)=b(9);
A2(19,37)=b(28);
A2(20,8)=f(33)+b(34);
A2(20,12)=b(35);
A2(20,20)=-b(33)-f(34)-f(35);
A2(21,8)=b(40);
A2(21,9)=b(20);
A2(21,21)=-b(39)-f(20)-f(40);
A2(21,38)=f(39);
A2(22,9)=f(13);
A2(22,22)=-b(13)-f(17);
A2(23,9)=f(41);
A2(23,23)=-b(41)-f(42);
A2(24,10)=b(11);
A2(24,24)=-b(10)-f(11)-f(15);
A2(24,39)=f(10);
A2(25,11)=f(37);
A2(25,25)=-b(37)-f(38);
A2(25,41)=b(38);
A2(26,12)=f(43)+b(41);
A2(26,26)=-b(43)-b(45)-f(41)-f(46)-f(49);
A2(26,42)=b(46);
A2(27,13)=b(23);
A2(27,27)=-b(22)-f(23);
A2(27,33)=f(22);
A2(28,13)=f(24);
A2(28,28)=-b(24)-f(27)-f(29);
A2(28,43)=b(27);
A2(29,14)=f(35);
A2(29,29)=-b(35)-b(41)-f(36)-f(43);
A2(29,44)=f(41)+b(43);
A2(29,45)=b(36);
A2(30,15)=b(4);
A2(30,30)=-b(3)-f(4);
A2(31,15)=f(22);
A2(31,31)=-b(22)-f(23);
A2(32,16)=f(7);
A2(32,32)=-b(7)-f(8)-f(14);
A2(32,46)=b(8);
A2(33,17)=b(5);
A2(33,27)=b(22);
A2(33,33)=-b(4)-f(5)-f(18)-f(22);
A2(33,47)=f(4);
A2(34,19)=b(7);
A2(34,34)=-b(6)-f(7);
A2(34,37)=f(6);
A2(35,19)=b(26)+b(27);
A2(35,35)=-b(25)-f(26)-f(27);
A2(35,48)=f(25);
A2(36,19)=f(9);
A2(36,36)=-b(9)-f(10);
A2(36,49)=b(10);
A2(37,19)=f(28);
A2(37,34)=b(6);
A2(37,37)=-b(5)-b(26)-b(28)-f(6);
A2(37,50)=f(5);
A2(37,51)=f(26);
A2(38,21)=b(39);
A2(38,38)=-b(38)-f(39);
A2(38,52)=f(38);
A2(39,24)=b(10);
A2(39,39)=-b(9)-f(10);
A2(39,46)=f(9);
A2(40,40)=-b(37)-f(38);
A2(40,45)=f(37);
A2(40,53)=b(38);
A2(41,25)=f(38);
A2(41,41)=-b(38)-f(39);
A2(41,54)=b(39);
A2(42,26)=f(46);
A2(42,42)=-b(44)-b(46)-f(47);
A2(43,28)=f(27);
A2(43,43)=-b(27)-f(28);
A2(44,29)=f(43)+b(41);
A2(44,44)=-b(43)-b(45)-f(41)-f(46)-f(49);
A2(44,55)=b(46);
A2(45,29)=f(36);
A2(45,40)=b(37);
A2(45,45)=-b(36)-f(37);
A2(46,32)=f(8);
A2(46,39)=b(9);
A2(46,46)=-b(7)-b(8)-b(26)-b(27)-f(9)-f(28);
A2(46,56)=f(7);
A2(46,57)=f(26)+f(27);
A2(46,58)=b(28);
A2(47,33)=b(4);
A2(47,47)=-b(3)-f(4);
A2(48,35)=b(25);
A2(48,48)=-b(23)-f(24)-f(25);
A2(48,59)=f(23);
A2(48,60)=b(24);
A2(49,36)=f(10);
A2(49,49)=-b(10)-f(11)-f(15);
A2(49,61)=b(11);
A2(50,37)=b(5);
A2(50,50)=-b(4)-f(5)-f(18)-f(22);
A2(50,59)=b(22);
A2(50,62)=f(4);
A2(51,37)=b(26);
A2(51,51)=-b(28)-f(26);
A2(51,63)=f(28);
A2(52,38)=b(38);
A2(52,52)=-b(37)-f(38);
A2(52,64)=f(37);
A2(53,40)=f(38);
A2(53,53)=-b(38)-f(39);
A2(53,65)=b(39);
A2(54,41)=f(39);
A2(54,54)=-b(39)-f(20)-f(40);
A2(54,66)=b(20);
A2(54,67)=b(40);
A2(55,44)=f(46);
A2(55,55)=-b(44)-b(46)-f(47);
A2(56,46)=b(7);
A2(56,56)=-b(6)-f(7);
A2(56,58)=f(6);
A2(57,46)=b(26)+b(27);
A2(57,57)=-b(25)-f(26)-f(27);
A2(57,68)=f(25);
A2(58,46)=f(28);
A2(58,56)=b(6);
A2(58,58)=-b(5)-b(26)-b(28)-f(6);
A2(58,69)=f(5);
A2(58,70)=f(26);
A2(59,48)=b(23);
A2(59,50)=f(22);
A2(59,59)=-b(22)-f(23);
A2(60,48)=f(24);
A2(60,60)=-b(24)-f(27)-f(29);
A2(60,63)=b(27);
A2(61,49)=f(11);
A2(61,61)=-b(11)-b(19)-f(12);
A2(61,66)=b(12);
A2(61,67)=f(19);
A2(62,50)=b(4);
A2(62,62)=-b(3)-f(4);
A2(63,51)=b(28);
A2(63,60)=f(27);
A2(63,63)=-b(27)-f(28);
A2(64,52)=b(37);
A2(64,64)=-b(36)-f(37);
A2(64,71)=f(36);
A2(65,53)=f(39);
A2(65,65)=-b(39)-f(20)-f(40);
A2(65,72)=b(40);
A2(66,54)=f(20);
A2(66,61)=f(12);
A2(66,66)=-b(12)-b(20)-f(13)-f(16)-f(21)-f(31)-f(41);
A2(66,67)=b(21);
A2(66,73)=b(13);
A2(66,74)=b(41);
A2(67,54)=f(40);
A2(67,61)=b(19);
A2(67,66)=f(21);
A2(67,67)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A2(67,75)=f(34)+b(33);
A2(68,57)=b(25);
A2(68,68)=-b(23)-f(24)-f(25);
A2(68,76)=f(23);
A2(68,77)=b(24);
A2(69,58)=b(5);
A2(69,69)=-b(4)-f(5)-f(18)-f(22);
A2(69,76)=b(22);
A2(69,78)=f(4);
A2(70,58)=b(26);
A2(70,70)=-b(28)-f(26);
A2(70,79)=f(28);
A2(71,64)=b(36);
A2(71,71)=-b(35)-b(41)-f(36)-f(43);
A2(71,80)=f(35);
A2(71,81)=f(41)+b(43);
A2(72,65)=f(40);
A2(72,72)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A2(72,82)=f(34)+b(33);
A2(73,66)=f(13);
A2(73,73)=-b(13)-f(17);
A2(74,66)=f(41);
A2(74,74)=-b(41)-f(42);
A2(75,67)=f(33)+b(34);
A2(75,75)=-b(33)-f(34)-f(35);
A2(76,68)=b(23);
A2(76,69)=f(22);
A2(76,76)=-b(22)-f(23);
A2(77,68)=f(24);
A2(77,77)=-b(24)-f(27)-f(29);
A2(77,79)=b(27);
A2(78,69)=b(4);
A2(78,78)=-b(3)-f(4);
A2(79,70)=b(28);
A2(79,77)=f(27);
A2(79,79)=-b(27)-f(28);
A2(80,71)=b(35);
A2(80,80)=-b(33)-f(34)-f(35);
A2(81,71)=f(43)+b(41);
A2(81,81)=-b(43)-b(45)-f(41)-f(46)-f(49);
A2(81,83)=b(46);
A2(82,72)=f(33)+b(34);
A2(82,82)=-b(33)-f(34)-f(35);
A2(82,84)=b(35);
A2(83,81)=f(46);
A2(83,83)=-b(44)-b(46)-f(47);
A2(84,82)=f(35);
A2(84,84)=-b(35)-b(41)-f(36)-f(43);
A2(84,85)=f(41)+b(43);
A2(85,84)=f(43)+b(41);
A2(85,85)=-b(43)-b(45)-f(41)-f(46)-f(49);
A2(85,86)=b(46);
A2(86,85)=f(46);
A2(86,86)=-b(44)-b(46)-f(47);

B2=zeros(86,28);
B2(1,2)=-f(32);
B2(6,1)=-f(27);
B2(8,4)=-f(48);
B2(17,3)=-f(26);
B2(26,6)=-f(45);
B2(30,8)=-f(3);
B2(31,5)=-b(23);
B2(42,9)=-f(44);
B2(43,7)=-b(28);
B2(44,10)=-f(45);
B2(47,11)=-f(3);
B2(55,14)=-f(44);
B2(62,15)=-f(3);
B2(65,12)=-b(20);
B2(66,13)=-b(31);
B2(67,19)=-f(48);
B2(72,16)=-f(21);
B2(72,17)=-b(19);
B2(72,22)=-f(48);
B2(75,18)=-b(35);
B2(78,23)=-f(3);
B2(80,20)=-f(33);
B2(80,21)=-b(34);
B2(81,24)=-f(45);
B2(83,26)=-f(44);
B2(84,25)=-b(36);
B2(85,27)=-f(45);
B2(86,28)=-f(44);

y2=[ cauchy(x1(3,:),x1(9,:)); in.AC_IN_1_2; cauchy(x1(3,:),x1(8,:)); in.OAA_IN_2_3; cauchy(x1(1,:),x1(13,:)); in.Glu_IN_2_3; cauchy(x1(2,:),x1(8,:)); in.GLC_IN_1_2; in.Gln_IN_2_3; in.Glu_IN_4_5; in.GLC_IN_2_3; cauchy(x1(1,:),x1(14,:)); cauchy(x1(1,:),x1(15,:)); in.Gln_IN_4_5; in.GLC_IN_4_5; cauchy(x1(1,:),x1(14,:)); cauchy(x1(1,:),x1(16,:)); cauchy(x1(1,:),x1(17,:)); in.OAA_IN_1_2; cauchy(x1(10,:),x1(12,:)); cauchy(x1(11,:),x1(12,:)); in.OAA_IN_3_4; in.GLC_IN_5_6; in.Glu_IN_3_4; cauchy(x1(1,:),x1(18,:)); in.Gln_IN_3_4; in.Glu_IN_1_2; in.Gln_IN_1_2];
x2=sparse(A2)\(B2*sparse(y2));

%%
A3=zeros(34);
A3(1,1)=-b(7)-b(8)-b(26)-b(27)-f(9)-f(28);
A3(1,2)=f(8);
A3(1,3)=b(9);
A3(1,11)=f(7);
A3(1,12)=f(26)+f(27);
A3(1,13)=b(28);
A3(2,1)=b(8);
A3(2,2)=-b(7)-f(8)-f(14);
A3(2,14)=f(7);
A3(3,1)=f(9);
A3(3,3)=-b(9)-f(10);
A3(3,4)=b(10);
A3(4,3)=f(10);
A3(4,4)=-b(10)-f(11)-f(15);
A3(4,5)=b(11);
A3(5,4)=f(11);
A3(5,5)=-b(11)-b(19)-f(12);
A3(5,6)=b(12);
A3(5,15)=f(19);
A3(6,5)=f(12);
A3(6,6)=-b(12)-b(20)-f(13)-f(16)-f(21)-f(31)-f(41);
A3(6,7)=b(13);
A3(6,8)=b(41);
A3(6,15)=b(21);
A3(6,16)=f(20);
A3(7,6)=f(13);
A3(7,7)=-b(13)-f(17);
A3(8,6)=f(41);
A3(8,8)=-b(41)-f(42);
A3(9,9)=-b(5)-b(26)-b(28)-f(6);
A3(9,10)=f(28);
A3(9,14)=b(6);
A3(9,17)=f(5);
A3(10,9)=b(28);
A3(10,10)=-b(27)-f(28);
A3(11,1)=b(7);
A3(11,11)=-b(6)-f(7);
A3(11,13)=f(6);
A3(12,1)=b(26)+b(27);
A3(12,12)=-b(25)-f(26)-f(27);
A3(12,18)=f(25);
A3(13,1)=f(28);
A3(13,11)=b(6);
A3(13,13)=-b(5)-b(26)-b(28)-f(6);
A3(13,19)=f(5);
A3(13,20)=f(26);
A3(14,2)=b(7);
A3(14,9)=f(6);
A3(14,14)=-b(6)-f(7);
A3(15,5)=b(19);
A3(15,6)=f(21);
A3(15,15)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A3(15,16)=f(40);
A3(15,21)=f(34)+b(33);
A3(16,6)=b(20);
A3(16,15)=b(40);
A3(16,16)=-b(39)-f(20)-f(40);
A3(16,22)=f(39);
A3(17,9)=b(5);
A3(17,17)=-b(4)-f(5)-f(18)-f(22);
A3(17,23)=f(4);
A3(17,24)=b(22);
A3(18,12)=b(25);
A3(18,18)=-b(23)-f(24)-f(25);
A3(18,25)=f(23);
A3(18,26)=b(24);
A3(19,13)=b(5);
A3(19,19)=-b(4)-f(5)-f(18)-f(22);
A3(19,25)=b(22);
A3(19,27)=f(4);
A3(20,13)=b(26);
A3(20,20)=-b(28)-f(26);
A3(20,28)=f(28);
A3(21,15)=f(33)+b(34);
A3(21,21)=-b(33)-f(34)-f(35);
A3(22,16)=b(39);
A3(22,22)=-b(38)-f(39);
A3(22,29)=f(38);
A3(23,17)=b(4);
A3(23,23)=-b(3)-f(4);
A3(24,17)=f(22);
A3(24,24)=-b(22)-f(23);
A3(25,18)=b(23);
A3(25,19)=f(22);
A3(25,25)=-b(22)-f(23);
A3(26,18)=f(24);
A3(26,26)=-b(24)-f(27)-f(29);
A3(26,28)=b(27);
A3(27,19)=b(4);
A3(27,27)=-b(3)-f(4);
A3(28,20)=b(28);
A3(28,26)=f(27);
A3(28,28)=-b(27)-f(28);
A3(29,22)=b(38);
A3(29,29)=-b(37)-f(38);
A3(29,30)=f(37);
A3(30,29)=b(37);
A3(30,30)=-b(36)-f(37);
A3(30,31)=f(36);
A3(31,30)=b(36);
A3(31,31)=-b(35)-b(41)-f(36)-f(43);
A3(31,32)=f(35);
A3(31,33)=f(41)+b(43);
A3(32,31)=b(35);
A3(32,32)=-b(33)-f(34)-f(35);
A3(33,31)=f(43)+b(41);
A3(33,33)=-b(43)-b(45)-f(41)-f(46)-f(49);
A3(33,34)=b(46);
A3(34,33)=f(46);
A3(34,34)=-b(44)-b(46)-f(47);

B3=zeros(34,12);
B3(6,1)=-b(31);
B3(9,2)=-f(26);
B3(10,3)=-f(27);
B3(15,5)=-f(48);
B3(21,4)=-b(35);
B3(23,7)=-f(3);
B3(24,6)=-b(23);
B3(27,8)=-f(3);
B3(32,9)=-f(33);
B3(32,10)=-b(34);
B3(33,11)=-f(45);
B3(34,12)=-f(44);

y3=[ cauchy(x1(1,:),x2(1,:)); cauchy(x2(2,:),x1(8,:)); cauchy(x2(2,:),x1(9,:)); cauchy(x1(1,:),x2(12,:)); in.OAA_IN_1_2_3; cauchy(x1(1,:),x2(13,:)); in.GLC_IN_1_2_3; in.GLC_IN_4_5_6; cauchy(x1(10,:),x2(8,:)); cauchy(x1(11,:),x2(8,:)); in.Glu_IN_2_3_4; in.Gln_IN_2_3_4];
x3=sparse(A3)\(B3*sparse(y3));

%%
A4=zeros(13);
A4(1,1)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A4(1,3)=f(34)+b(33);
A4(1,4)=f(40);
A4(2,2)=-b(35)-b(41)-f(36)-f(43);
A4(2,3)=f(35);
A4(2,5)=f(41)+b(43);
A4(3,1)=f(33)+b(34);
A4(3,2)=b(35);
A4(3,3)=-b(33)-f(34)-f(35);
A4(4,1)=b(40);
A4(4,4)=-b(39)-f(20)-f(40);
A4(4,6)=f(39);
A4(5,2)=f(43)+b(41);
A4(5,5)=-b(43)-b(45)-f(41)-f(46)-f(49);
A4(5,7)=b(46);
A4(6,4)=b(39);
A4(6,6)=-b(38)-f(39);
A4(6,8)=f(38);
A4(7,5)=f(46);
A4(7,7)=-b(44)-b(46)-f(47);
A4(8,6)=b(38);
A4(8,8)=-b(37)-f(38);
A4(8,9)=f(37);
A4(9,8)=b(37);
A4(9,9)=-b(36)-f(37);
A4(9,10)=f(36);
A4(10,9)=b(36);
A4(10,10)=-b(35)-b(41)-f(36)-f(43);
A4(10,11)=f(35);
A4(10,12)=f(41)+b(43);
A4(11,10)=b(35);
A4(11,11)=-b(33)-f(34)-f(35);
A4(12,10)=f(43)+b(41);
A4(12,12)=-b(43)-b(45)-f(41)-f(46)-f(49);
A4(12,13)=b(46);
A4(13,12)=f(46);
A4(13,13)=-b(44)-b(46)-f(47);

B4=zeros(13,11);
B4(1,1)=-f(21);
B4(1,2)=-b(19);
B4(1,5)=-f(48);
B4(2,3)=-b(36);
B4(4,4)=-b(20);
B4(5,6)=-f(45);
B4(7,7)=-f(44);
B4(11,8)=-f(33);
B4(11,9)=-b(34);
B4(12,10)=-f(45);
B4(13,11)=-f(44);

y4=[ cauchy(x1(1,:),x2(9,:)); cauchy(x1(1,:),x2(10,:)); cauchy(x1(1,:),x2(11,:)); cauchy(x1(1,:),x2(9,:)); in.OAA_IN_2_3_4; in.Glu_IN_1_2_3; in.Gln_IN_1_2_3; cauchy(x2(1,:),x1(12,:)); cauchy(x2(3,:),x1(12,:)); in.Glu_IN_3_4_5; in.Gln_IN_3_4_5];
x4=sparse(A4)\(B4*sparse(y4));

%%
A5=zeros(10);
A5(1,1)=-b(21)-b(34)-b(40)-b(48)-f(19)-f(30)-f(33);
A5(1,2)=f(40);
A5(1,6)=f(34)+b(33);
A5(2,1)=b(40);
A5(2,2)=-b(39)-f(20)-f(40);
A5(2,5)=f(39);
A5(3,3)=-b(36)-f(37);
A5(3,4)=b(37);
A5(3,7)=f(36);
A5(4,3)=f(37);
A5(4,4)=-b(37)-f(38);
A5(4,5)=b(38);
A5(5,2)=b(39);
A5(5,4)=f(38);
A5(5,5)=-b(38)-f(39);
A5(6,1)=f(33)+b(34);
A5(6,6)=-b(33)-f(34)-f(35);
A5(7,3)=b(36);
A5(7,7)=-b(35)-b(41)-f(36)-f(43);
A5(7,8)=f(35);
A5(7,9)=f(41)+b(43);
A5(8,7)=b(35);
A5(8,8)=-b(33)-f(34)-f(35);
A5(9,7)=f(43)+b(41);
A5(9,9)=-b(43)-b(45)-f(41)-f(46)-f(49);
A5(9,10)=b(46);
A5(10,9)=f(46);
A5(10,10)=-b(44)-b(46)-f(47);

B5=zeros(10,9);
B5(1,1)=-f(21);
B5(1,2)=-b(19);
B5(1,5)=-f(48);
B5(2,3)=-b(20);
B5(6,4)=-b(35);
B5(8,6)=-f(33);
B5(8,7)=-b(34);
B5(9,8)=-f(45);
B5(10,9)=-f(44);

y5=[ cauchy(x1(1,:),x3(6,:)); cauchy(x1(1,:),x3(5,:)); cauchy(x1(1,:),x3(6,:)); cauchy(x1(1,:),x4(2,:)); in.OAA_IN_1_2_3_4; cauchy(x2(1,:),x2(8,:)); cauchy(x2(3,:),x2(8,:)); in.Glu_IN_2_3_4_5; in.Gln_IN_2_3_4_5];
x5=A5\(B5*sparse(y5));

%%
A6=zeros(10);
A6(1,1)=-b(28)-f(26);
A6(1,2)=f(28);
A6(1,3)=b(26);
A6(2,1)=b(28);
A6(2,2)=-b(27)-f(28);
A6(2,4)=f(27);
A6(3,1)=f(26);
A6(3,3)=-b(5)-b(26)-b(28)-f(6);
A6(3,5)=f(5);
A6(3,6)=b(6);
A6(4,2)=b(27);
A6(4,4)=-b(24)-f(27)-f(29);
A6(4,7)=f(24);
A6(5,3)=b(5);
A6(5,5)=-b(4)-f(5)-f(18)-f(22);
A6(5,8)=f(4);
A6(5,9)=b(22);
A6(6,3)=f(6);
A6(6,6)=-b(6)-f(7);
A6(7,4)=b(24);
A6(7,7)=-b(23)-f(24)-f(25);
A6(7,9)=f(23);
A6(7,10)=b(25);
A6(8,5)=b(4);
A6(8,8)=-b(3)-f(4);
A6(9,5)=f(22);
A6(9,7)=b(23);
A6(9,9)=-b(22)-f(23);
A6(10,7)=f(25);
A6(10,10)=-b(25)-f(26)-f(27);

B6=zeros(10,5);
B6(3,1)=-f(28);
B6(6,2)=-b(7);
B6(8,5)=-f(3);
B6(10,3)=-b(26);
B6(10,4)=-b(27);

y6=[ cauchy(x1(4,:),x3(1,:)); cauchy(x1(5,:),x3(1,:)); cauchy(x1(6,:),x3(1,:)); cauchy(x1(7,:),x3(1,:)); in.GLC_IN_3_4_5_6];
x6=A6\(B6*sparse(y6));

%%
A7=zeros(9);
A7(1,1)=-b(23)-f(24)-f(25);
A7(1,2)=b(24);
A7(1,3)=b(25);
A7(1,4)=f(23);
A7(2,1)=f(24);
A7(2,2)=-b(24)-f(27)-f(29);
A7(2,5)=b(27);
A7(3,1)=f(25);
A7(3,3)=-b(25)-f(26)-f(27);
A7(4,1)=b(23);
A7(4,4)=-b(22)-f(23);
A7(4,6)=f(22);
A7(5,2)=f(27);
A7(5,5)=-b(27)-f(28);
A7(6,4)=b(22);
A7(6,6)=-b(4)-f(5)-f(18)-f(22);
A7(6,7)=f(4);
A7(6,8)=b(5);
A7(7,6)=b(4);
A7(7,7)=-b(3)-f(4);
A7(8,6)=f(5);
A7(8,8)=-b(5)-b(26)-b(28)-f(6);
A7(8,9)=b(6);
A7(9,8)=f(6);
A7(9,9)=-b(6)-f(7);

B7=zeros(9,7);
B7(3,1)=-b(26);
B7(3,2)=-b(27);
B7(5,3)=-b(28);
B7(7,6)=-f(3);
B7(8,4)=-f(26);
B7(8,5)=-f(28);
B7(9,7)=-b(7);

y7=[ cauchy(x2(4,:),x3(1,:)); cauchy(x2(5,:),x3(1,:)); cauchy(x1(2,:),x6(1,:)); cauchy(x1(3,:),x6(1,:)); cauchy(x2(6,:),x3(1,:)); in.GLC_IN_2_3_4_5_6; cauchy(x2(7,:),x3(1,:))];
x7=A7\(B7*sparse(y7));

%%
A8=zeros(4);
A8(1,1)=-b(35)-b(41)-f(36)-f(43);
A8(1,2)=f(41)+b(43);
A8(1,4)=f(35);
A8(2,1)=f(43)+b(41);
A8(2,2)=-b(43)-b(45)-f(41)-f(46)-f(49);
A8(2,3)=b(46);
A8(3,2)=f(46);
A8(3,3)=-b(44)-b(46)-f(47);
A8(4,1)=b(35);
A8(4,4)=-b(33)-f(34)-f(35);

B8=zeros(4,5);
B8(1,1)=-b(36);
B8(2,4)=-f(45);
B8(3,5)=-f(44);
B8(4,2)=-f(33);
B8(4,3)=-b(34);

y8=[ cauchy(x1(1,:),x5(3,:)); cauchy(x2(1,:),x4(1,:)); cauchy(x2(3,:),x4(1,:)); in.Glu_IN_1_2_3_4_5; in.Gln_IN_1_2_3_4_5];
x8=A8\(B8*sparse(y8));

%%
A9=zeros(5);
A9(1,1)=-b(3)-f(4);
A9(1,2)=b(4);
A9(2,1)=f(4);
A9(2,2)=-b(4)-f(5)-f(18)-f(22);
A9(2,3)=b(5);
A9(2,5)=b(22);
A9(3,2)=f(5);
A9(3,3)=-b(5)-b(26)-b(28)-f(6);
A9(3,4)=b(6);
A9(4,3)=f(6);
A9(4,4)=-b(6)-f(7);
A9(5,2)=f(22);
A9(5,5)=-b(22)-f(23);

B9=zeros(5,5);
B9(1,5)=-f(3);
B9(3,1)=-f(26);
B9(3,2)=-f(28);
B9(4,3)=-b(7);
B9(5,4)=-b(23);

y9=[ cauchy(x2(2,:),x6(1,:)); cauchy(x3(10,:),x3(1,:)); cauchy(x3(2,:),x3(1,:)); cauchy(x1(1,:),x7(1,:)); in.GLC_IN_1_2_3_4_5_6];
x9=A9\(B9*sparse(y9));

%%
A10=zeros(1);
A10(1,1)=-b(33)-f(34)-f(35);

B10=zeros(1,3);
B10(1,1)=-f(33);
B10(1,2)=-b(34);
B10(1,3)=-b(35);

y10=[ cauchy(x2(1,:),x5(1,:)); cauchy(x2(3,:),x5(1,:)); cauchy(x1(1,:),x8(1,:))];
x10=A10\(B10*sparse(y10));

%%
A11=zeros(1);
A11(1,1)=-b(27)-f(28);

B11=zeros(1,2);
B11(1,1)=-f(27);
B11(1,2)=-b(28);

y11=[ cauchy(x2(2,:),x7(2,:)); cauchy(x3(9,:),x6(1,:))];
x11=A11\(B11*sparse(y11));

%%
mol.CO2=x1(1,:);
mol.GLC=x9(1,:);
mol.G6P=x9(2,:);
mol.F6P=x9(3,:);
mol.FBP=x9(4,:);
mol.GAP=x3(1,:);
mol.DHAP=x3(2,:);
mol.BPG=x3(3,:);
mol.PGA=x3(4,:);
mol.PEP=x3(5,:);
mol.PYR=x3(6,:);
mol.LAC=x3(7,:);
mol.OAA=x5(1,:);
mol.MAL=x5(2,:);
mol.m6PG=x9(5,:);
mol.Ru5P=x7(1,:);
mol.R5P=x7(2,:);
mol.X5P=x7(3,:);
mol.S7P=x11(1,:);
mol.E4P=x6(1,:);
mol.AcCoA=x2(1,:);
mol.CitICit=x10(1,:);
mol.OGA=x8(1,:);
mol.SuccCoA=x5(3,:);
mol.Succ=x5(4,:);
mol.Fum=x5(5,:);
mol.Glu=x8(2,:);
mol.Ala=x3(8,:);
mol.Gln=x8(3,:);