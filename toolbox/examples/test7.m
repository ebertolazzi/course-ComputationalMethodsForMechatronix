addpath('../lib');
close all;
hold off;

max_iter = 400;

PB7 = Problem7( 4 );
npt = 10;
S   = Direct_Transcription_XY( PB7, npt );

ok = S.solve( max_iter, -10, -10, 10, 10 );
subplot(2,2,1);
S.plot();
axis equal;

[t,x,y] = S.get_solution();
xx = zeros( 2*npt-1, 1 );
xx(1:2:end)   = x;
xx(2:2:end-1) = (x(1:end-1)+x(2:end))/2;
yy = zeros( 2*npt-1, 1 );
yy(1:2:end)   = y;
yy(2:2:end-1) = (y(1:end-1)+y(2:end))/2;

npt = 2*npt-1;
ok = S.solve2( npt, xx, yy, max_iter, -10, -10, 10, 10 );
subplot(2,2,2);
S.plot();
axis equal;


[t,x,y] = S.get_solution();
xx = zeros( 2*npt-1, 1 );
xx(1:2:end)   = x;
xx(2:2:end-1) = (x(1:end-1)+x(2:end))/2;
yy = zeros( 2*npt-1, 1 );
yy(1:2:end)   = y;
yy(2:2:end-1) = (y(1:end-1)+y(2:end))/2;

npt = 2*npt-1;
ok = S.solve2( npt, xx, yy, max_iter, -10, -10, 10, 10 );
subplot(2,2,3);
S.plot();
axis equal;

[t,x,y] = S.get_solution();
xx = zeros( 2*npt-1, 1 );
xx(1:2:end)   = x;
xx(2:2:end-1) = (x(1:end-1)+x(2:end))/2;
yy = zeros( 2*npt-1, 1 );
yy(1:2:end)   = y;
yy(2:2:end-1) = (y(1:end-1)+y(2:end))/2;

npt = 2*npt-1;
ok = S.solve2( npt, xx, yy, max_iter, -10, -10, 10, 10 );
subplot(2,2,4);
S.plot();
axis equal;
