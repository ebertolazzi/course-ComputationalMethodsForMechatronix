addpath('../lib');
close all;


PB4 = Problem4( 0, 1, pi ) %0.9*pi/2 );
S   = Direct_Transcription_1D( PB4, 50 );

S.solve( 400, 0, 1 );
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB4.exact(t), 'Linewidth', 3 );
axis equal;