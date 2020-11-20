addpath('../lib');
close all;

PB3 = Problem3( 1 );
S   = Direct_Transcription_1D( PB3, 20 );

S.solve( 400, -10, 0 );
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB3.exact(t), 'Linewidth', 3 );

