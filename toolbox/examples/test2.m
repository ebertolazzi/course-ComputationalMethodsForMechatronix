addpath('../lib');
close all;

PB2 = Problem2();
S   = Direct_Transcription_1D( PB2, 30 );

S.solve( 1000, 0, 3 );
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB2.exact(t), 'Linewidth', 3 );

