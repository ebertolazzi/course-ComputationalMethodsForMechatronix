addpath('../lib');
close all;

PB1 = Problem1();
S   = Direct_Transcription_1D( PB1, 40 ); % PROBLEM, NUMBER OF DISCRETIZATION POINTS

S.solve( 1000, 0, 3 ); % max iter, lower, upper
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB1.exact(t), 'Linewidth', 3 );

