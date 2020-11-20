addpath('../lib');
close all;


PB6 = Problem6( 0, 1, 1, 1, 10 );
S   = Direct_Transcription_1D( PB6, 20 );

S.solve( 400, -5, 1 );
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB6.exact(t), 'Linewidth', 3 );
axis equal;