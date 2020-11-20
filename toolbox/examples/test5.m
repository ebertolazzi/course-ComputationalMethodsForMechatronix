addpath('../lib');
close all;


PB5 = Problem5();
S   = Direct_Transcription_1D( PB5, 40 );

S.solve( 400, -Inf, Inf );
S.plot();
hold on;
t = 0:0.001:1;
plot( t, PB5.exact(t), 'Linewidth', 3 );
axis equal;