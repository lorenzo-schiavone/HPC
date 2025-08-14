clc; clear; close all
% A = [8.219900  5.204082  8.249664  5.618524;
% 4.630122  3.205535  4.689079  3.162474;
% 4.060645  2.878401  4.205213  2.737441;
% 3.634686  2.751820  3.705338  2.398958;];
% 
% %A = A';
% plot(A, '-s')
% legend({'NO1', 'NO2', 'RCM1', 'RCM2'}, 'FontSize',12)
% ylabel("time (sec)")
% xlabel("np")

% A = [1         1      1   1   ;
%      1.77   1.62   1.76   1.78;
%      2.03   1.81   1.96   2.05;
%      2.26   1.89   2.23   2.34;
% ];
% 
% plot(A, '-s')
% legend({'NO1', 'NO2', 'RCM1', 'RCM2'}, 'FontSize',12)
% ylabel("SpeedUp")
% xlabel("np")

A = [180.900720  141.801716;
    148.415983  105.693213 ;
    151.490037  103.072335 ;
    153.800005  97.531304  ;
];
figure
plot(A, '-s')
legend({'fullQR', 'Givens'}, 'FontSize',12)
ylabel("time (sec)")
xlabel("np")

A =[27.755223   28.141819; 
    20.456823   19.665255; 
    19.955573   19.751338; 
    18.681080   18.420930; ];
figure
plot(A, '-s')
legend({'fullQR(20)', 'Givens(20)'}, 'FontSize',12)
ylabel("time (sec)")
xlabel("np")

A = [1     1      1      1   ;
    1.22   1.34   1.35   1.43 ;
    1.19   1.37   1.39   1.42 ;
    1.18   1.45   1.48   1.52 ;   
];
figure
plot(A, '-s')
legend({'fullQR', 'Givens', 'fullQR(20)', 'Givens(20)'}, 'FontSize',12)
ylabel("SpeedUp")
xlabel("np")