addpath('../lib');
close all;


PB6 = Problem6( 0, 1, 1, 1, 2 );
S   = Direct_Transcription_1D( PB6, 40 );

S.solve( 400, -1, 1 );
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB6.exact(t), 'Linewidth', 3 );
axis equal;