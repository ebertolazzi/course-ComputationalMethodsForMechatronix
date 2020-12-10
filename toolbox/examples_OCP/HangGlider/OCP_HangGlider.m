%--------------------------------------------------------------------------%
%                                                                          %
%  Copyright (C) 2020                                                      %
%                                                                          %
%         , __                 , __                                        %
%        /|/  \               /|/  \                                       %
%         | __/ _   ,_         | __/ _   ,_                                %
%         |   \|/  /  |  |   | |   \|/  /  |  |   |                        %
%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       %
%                           /|                   /|                        %
%                           \|                   \|                        %
%                                                                          %
%      Enrico Bertolazzi                                                   %
%      Dipartimento di Ingegneria Industriale                              %
%      Universita` degli Studi di Trento                                   %
%      email: enrico.bertolazzi@unitn.it                                   %
%                                                                          %
%--------------------------------------------------------------------------%

classdef OCP_HangGlider < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    % PARAMETERS OF THE PROBLEM
    rc
    uc
    c0
    c1
    S
    rho
    m
    g
    cmax
    %
    x_i
    y_i
    y_f
    vx_i
    vx_f
    vy_i
    vy_f
    %
    epsi
    epsi1
  end

  methods ( Hidden = true )

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %
    function res = r( self, x )
      rc  = self.rc;
      res = (x/rc-2.5)^2;
    end
    %
    function res = r_D( self, x )
      rc  = self.rc;
      res = 2*(x/rc-2.5)/rc;
    end
    %
    function res = r_DD( self, x )
      rc  = self.rc;
      res = 2/rc^2;
    end
    %
    function res = u( self, x )
      uc  = self.uc;
      rx  = self.r(x);
      res = uc*(1-rx)*exp(-rx);
    end
    %
    function res = u_D( self, x )
      uc  = self.uc;
      rx  = self.r(x);
      t1  = self.r_D(x);
      t3  = self.r(x);
      t4  = exp(-t3);
      res = (t3 - 2) * t4 * t1 * uc;
    end
    %
    function res = u_DD( self, x )
      uc  = self.uc;
      t1  = self.r(x);
      t3  = self.r_DD(x);
      t5  = self.r_D(x);
      t6  = t5 ^ 2;
      t11 = exp(-t1);
      res = -t11 * (t3 * (2-t1) + (t1 - 3) * t6) * uc;
    end
    %
    function res = w( self, x, yd )
      res = yd-self.u(x);
    end
    %
    function res = w_x( self, x, yd )
      res = -self.u_D(x);
    end
    %
    function res = w_yd( self, x, yd )
      res = 1;
    end
    %
    function res = w_xx( self, x, yd )
      res = -self.u_DD(x);
    end
    %
    function res = w_xyd( self, x, yd )
      res = 0;
    end
    %
    function res = w_ydyd( self, x, yd )
      res = 0;
    end
    %
    function res = v( self, x, xd, yd )
      res = hypot( xd, self.w(x,yd) );
    end
    %
    function res = v_x( self, x, xd, yd )
      w   = self.w(x,yd);
      w_x = self.w_x(x,yd);
      res = w*w_x/hypot( xd, w );
    end
    %
    function res = v_xd( self, x, xd, yd )
      w   = self.w(x,yd);
      res = xd/hypot( xd, w );
    end
    %
    function res = v_yd( self, x, xd, yd )
      w    = self.w(x,yd);
      w_yd = self.w_yd(x,yd);
      res  = w*w_yd/hypot( xd, w );
    end
    %
    function res = v_xx( self, x, xd, yd )
      w    = self.w(x,yd);
      w_x  = self.w_x(x,yd);
      w_xx = self.w_xx(x,yd);
      v    = hypot( xd, w );
      res  = (w*w_xx+(w_x*xd/v)^2)/v;
    end
    %
    function res = v_xxd( self, x, xd, yd )
      w    = self.w(x,yd);
      w_x  = self.w_x(x,yd);
      v    = hypot( xd, w );
      res  = -w*xd*w_x/v^3;
    end
    %
    function res = v_xyd( self, x, xd, yd )
      w     = self.w(x,yd);
      w_x   = self.w_x(x,yd);
      w_yd  = self.w_yd(x,yd);
      w_xxd = self.w_xxd(x,yd);
      v     = hypot( xd, w );
      res   = (w*w_xxd+w_x*w_yd*(xd/v)^2)/v;
    end
    %
    function res = v_xdxd( self, x, xd, yd )
      w   = self.w(x,yd);
      v   = hypot( xd, w );
      res = (w/v)^2/v;
    end
    %
    function res = v_xdyd( self, x, xd, yd )
      w    = self.w(x,yd);
      w_yd = self.w_yd(x,yd);
      v    = hypot( xd, w );
      res  = -(xd*w*w_yd)/v^3;
    end
    %
    function res = v_ydyd( self, x, xd, yd )
      w      = self.w(x,yd);
      w_yd   = self.w_yd(x,yd);
      w_ydyd = self.w_ydyd(x,yd);
      v      = hypot( xd, w );
      res    = (w*w_ydyd+(xd*w_yd/v)^2)/v;
    end
    %
    function res = SV2( self, x, xd, yd )
      w   = self.w(x,yd);
      v   = self.v(x,xd,yd);
      res = w*v;
    end
    %
    function res = SV2_x( self, x, xd, yd )
      w   = self.w(x,yd);
      w_x = self.w_x(x,yd);
      v   = self.v(x,xd,yd);
      v_x = self.v_x(x,xd,yd);
      res = w_x*v+w*v_x;
    end
    %
    function res = SV2_xd( self, x, xd, yd )
      w    = self.w(x,yd);
      v_xd = self.v_xd(x,xd,yd);
      res  = w*v_xd;
    end
    %
    %
    function res = SV2_yd( self, x, xd, yd )
      w    = self.w(x,yd);
      w_yd = self.w_yd(x,yd);
      v    = self.v(x,xd,yd);
      v_yd = self.v_yd(x,xd,yd);
      res  = w_yd*v+w*v_yd;
    end
    %
    function res = SV2_xx( self, x, xd, yd )
      t1  = self.w_xx(x, yd);
      t2  = self.v(x, xd, yd);
      t4  = self.w_x(x, yd);
      t5  = self.v_x(x, xd, yd);
      t8  = self.w(x, yd);
      t9  = self.v_xx(x, xd, yd);
      res = t2 * t1 + 2 * t5 * t4 + t9 * t8;
    end
    %
    function res = SV2_xxd( self, x, xd, yd )
      t1  = self.w_x(x, yd);
      t2  = self.v+xd(x, xd, yd);
      t4  = self.w(x, yd);
      t5  = self.v_xxd(x, xd, yd);
      res = t2 * t1 + t5 * t4;
    end
    %
    function res = SV2_xyd( self, x, xd, yd )
      t1  = self.w_xyd(x, yd);
      t2  = self.v(x, xd, yd);
      t4  = self.w_x(x, yd);
      t5  = self.v_yd(x, xd, yd);
      t7  = self.w_yd(x, yd);
      t8  = self.v_x(x, xd, yd);
      t10 = self.w(x, yd);
      t11 = self.v_xyd(x, xd, yd);
      res = t2 * t1 + t11 * t10 + t5 * t4 + t8 * t7;
    end
    %
    function res = SV2_xdxd( self, x, xd, yd )
      t1  = self.w(x, yd);
      t2  = self.v_xdxd(x, xd, yd);
      res = t2 * t1;
    end
    %
    function res = SV2_xdyd( self, x, xd, yd )
      t1  = self.w_yd(x, yd);
      t2  = self.v_xd(x, xd, yd);
      t4  = self.w(x, yd);
      t5  = self.v_xdyd(x, xd, yd);
      res = t2 * t1 + t5 * t4;
    end
    %
    function res = SV2_ydyd( self, x, xd, yd )
      t1  = self.w_ydyd(x, yd);
      t2  = self.v(x, xd, yd);
      t4  = self.w_yd(x, yd);
      t5  = self.v_yd(x, xd, yd);
      t8  = self.w(x, yd);
      t9  = self.v_ydyd(x, xd, yd);
      res = t2 * t1 + 2 * t5 * t4 + t9 * t8;
    end
    %
    function res = CV2( self, x, xd, yd )
      v   = self.v(x,xd,yd);
      res = xd*v;
    end
    %
    function res = CV2_x( self, x, xd, yd )
      t1  = self.v_x(x, xd, yd);
      res = t1 * xd;
    end
    %
    function res = CV2_xd( self, x, xd, yd )
      t1  = self.v(x, xd, yd);
      t2  = self.v_xd(x, xd, yd);
      res = t2 * xd + t1;
    end
    %
    function res = CV2_yd( self, x, xd, yd )
      t1  = self.v_yd(x, xd, yd);
      res = t1 * xd;
    end
    %
    function res = CV2_xx( self, x, xd, yd )
      t1  = self.v_xx(x, xd, yd);
      res = t1 * xd;
    end
    %
    function res = CV2_xxd( self, x, xd, yd )
      t1  = self.v_x(x, xd, yd);
      t2  = self.v_xxd(x, xd, yd);
      res = t2 * xd + t1;
    end
    %
    function res = CV2_xyd( self, x, xd, yd )
      t1  = self.v_xyd(x, xd, yd);
      res = t1 * xd;
    end
    %
    function res = CV2_xdxd( self, x, xd, yd )
      t1 = self.v_xd(x, xd, yd);
      t3 = self.v_xdxd(x, xd, yd);
      res = t3 * xd + 2 * t1;
    end
    %
    function res = CV2_xdyd( self, x, xd, yd )
      t1  = self.v_yd(x, xd, yd);
      t2  = self.v_xdyd(x, xd, yd);
      res = t2 * xd + t1;
    end
    %
    function res = CV2_ydyd( self, x, xd, yd )
      t1  = self.v_ydyd(x, xd, yd);
      res = t1 * xd;
    end
    %
    function [resx,resy] = RHS_XY( self, x, xd, yd, cL )
      rho  = self.rho;
      S    = self.S;
      c0   = self.c0;
      c1   = self.c1;
      tmp  = 0.5*rho*S;

      D    = (c0+c1*cL^2)*tmp;
      L    = cL*tmp;

      SS   = self.SV2(x, xd, yd);
      CC   = self.CV2(x, xd, yd);

      resx = -L*SS-D*CC;
      resy = -D*SS+L*CC;
    end
    %
    function [gradx,grady] = RHS_XY_gradient( self, x, xd, yd, cL )
      rho  = self.rho;
      S    = self.S;
      c0   = self.c0;
      c1   = self.c1;
      tmp  = 0.5*rho*S;

      D    = (c0+c1*cL^2)*tmp;
      L    = cL*tmp;
      D_cL = 2*c1*cL*tmp;
      L_cL = tmp;

      SS   = self.SV2(x, xd, yd);
      S_x  = self.SV2_x(x, xd, yd);
      S_xd = self.SV2_xd(x, xd, yd);
      S_yd = self.SV2_yd(x, xd, yd);

      CC   = self.CV2(x, xd, yd);
      C_x  = self.CV2_x(x, xd, yd);
      C_xd = self.CV2_xd(x, xd, yd);
      C_yd = self.CV2_yd(x, xd, yd);

      gradx = -L*[S_x,S_xd,S_yd]-D*[C_x,C_xd,C_yd];
      grady = -D*[S_x,S_xd,S_yd]+L*[C_x,C_xd,C_yd];

      gradx = [ gradx, -L_cL*SS-D_cL*CC ];
      grady = [ grady, -D_cL*SS+L_cL*CC ];
    end

    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, P )
      x  = X(1);
      y  = X(2);
      vx = X(3);
      vy = X(4);
      cL = U(1);
      T  = P(1);
      m  = self.m;
      g  = self.g;
      [resx,resy] = self.RHS_XY( x, vx, vy, cL );
      RES    = zeros(4,1);
      RES(1) = T*vx;
      RES(2) = T*vy;
      RES(3) = (T/m)*resx;
      RES(4) = T*(resy/m-g);
    end
    %
    function RES = JAC( self, nseg, t, X, U, P )
      RES = sparse( self.nx, self.nx+self.nu+self.np );

      x  = X(1);
      y  = X(2);
      vx = X(3);
      vy = X(4);
      cL = U(1);
      T  = P(1);
      m  = self.m;
      g  = self.g;

      [resx,resy]   = self.RHS_XY( x, vx, vy, cL );
      [gradx,grady] = self.RHS_XY_gradient( x, vx, vy, cL );

      RES(1,3) = T;
      RES(1,6) = vx;

      RES(2,4) = T;
      RES(2,6) = vy;

      Tm = T/m;

      RES(3,1) = Tm*gradx(1);
      RES(3,3) = Tm*gradx(2);
      RES(3,4) = Tm*gradx(3);
      RES(3,5) = Tm*gradx(4);
      RES(3,6) = resx/m;

      RES(4,1) = Tm*grady(1);
      RES(4,3) = Tm*grady(2);
      RES(4,4) = Tm*grady(3);
      RES(4,5) = Tm*grady(4);
      RES(4,6) = resy/m-g;

    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( self.nx, self.nx+self.nu+self.np );

      RES(1,3) = 1;
      RES(1,6) = 1;

      RES(2,4) = 1;
      RES(2,6) = 1;

      RES(3,1) = 1;
      RES(3,3) = 1;
      RES(3,4) = 1;
      RES(3,5) = 1;
      RES(3,6) = 1;

      RES(4,1) = 1;
      RES(4,3) = 1;
      RES(4,4) = 1;
      RES(4,5) = 1;
      RES(4,6) = 1;
    end
  end

  methods

    function self = OCP_HangGlider( )
      nx  = 4; % number of states
      nu  = 1; % number of controls
      np  = 1; % number of free parameters
      nbc = 7; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      % setup the parameters of the model
      self.uc    = 2.5;
      self.rc    = 100;
      self.c0    = 0.034;
      self.c1    = 0.069662;
      self.S     = 14;
      self.rho   = 1.13;
      self.x_i   = 0;
      self.y_i   = 1000;
      self.y_f   = 900;
      self.vx_i  = 13.23;
      self.vx_f  = 13.23;
      self.vy_i  = -1.288;
      self.vy_f  = -1.288;
      self.m     = 100;
      self.g     = 9.81;
      self.cmax  = 1.4;
      self.epsi  = 0; % no penalization of total variations of the control
      self.epsi1 = 0;
    end

    function info = solve( self )

      % load discretization parmeters
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      N    = self.N;
      nseg = self.nseg;
      njmp = self.njmp;
      nbc  = self.nbc;

      utot = (N-nseg)*nu;
      xtot = N*nx;

      % set the lower bound of the states (nx x N)
      x_lb = [          ...
        zeros(1,N);     ...
        zeros(1,N); ...
        0.01*ones(1,N);     ...
        -Inf*ones(1,N); ...
      ];

      % set the upper bound of the states (nx x N)
      x_ub = [         ...
        Inf*ones(1,N); ...
        Inf*ones(1,N); ...
        Inf*ones(1,N); ...
        Inf*ones(1,N); ...
      ];

      % set the lower bound of the controls (nu x (N-nseg))
      u_lb = zeros(1, N-nseg);

      % set the upper bound of the controls (nu x (N-nseg))
      u_ub = self.cmax*ones(1, N-nseg);

      % set the lower bound of the optimization parameters (1 x np)
      p_lb = [0];

      % set the upper bound of the optimization parameters (1 x np)
      p_ub = [200];

      options.lb = self.pack( x_lb, u_lb, p_lb ); % Lower bound on the variables.
      options.ub = self.pack( x_ub, u_ub, p_ub ); % Upper bound on the variables.

      % The constraint functions are bounded to zero
      dim = (N-nseg)*nx + (nseg-1)*njmp + nbc;
      options.cl = zeros(1,dim); %  constraints
      options.cu = zeros(1,dim);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant   = 'no';
      options.ipopt.hessian_constant = 'no';
      options.ipopt.mu_strategy      = 'adaptive';
      options.ipopt.max_iter         = 800;
      options.ipopt.tol              = 1e-10;%
      options.ipopt.linear_solver    = 'mumps';

      % The callback functions.
      funcs.objective         = @(Z) self.NLP_target(Z);
      funcs.gradient          = @(Z) self.NLP_target_gradient(Z);

      funcs.constraints       = @(Z) self.NLP_constraints(Z);
      funcs.jacobian          = @(Z) self.NLP_constraints_jacobian(Z);
      funcs.jacobianstructure = @() self.NLP_constraints_jacobian_pattern();

      if false
        %options.ipopt.derivative_test = 'first-order';
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian           = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure  = @() self.NLP_hessian_pattern();
      else
        options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % GUESS SOLUTION

      Tguess   = 10;
      Vy_guess = (self.y_f-self.y_i)/Tguess;
      Xguess = [                                     ...
        self.x_i + self.vx_i * self.nodes * Tguess;  ...
        self.y_i + (self.y_f-self.y_i) * self.nodes; ...
        self.vx_i * ones(1,N);                       ...
        Vy_guess * ones(1,N);      ...
      ];
      % set the lower bound of the controls (nu x (N-nseg))
      Uguess = 0.5*self.cmax*ones(1, N-nseg);
      Pguess = Tguess;

      x0 = self.pack( Xguess, Uguess, Pguess ); % Lower bound on the variables.

      [self.sol, info] = ipopt(x0,funcs,options);
      self.sol(end)
    end

    function plot( self )
      % example of pliot solution
      nodes = self.nodes;
      X     = self.states();
      Tf    = self.sol(end);

      subplot( 2, 2, 1 );
      plot( Tf*nodes, X(2,:), '-o', 'Linewidth', 2 );
      title('y');

      subplot( 2, 2, 2 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Tf*Unodes, UU(1,:), '-r', 'Linewidth', 2 );
      title('cL');

      subplot( 2, 2, 3 );
      plot( Tf*nodes, X(3,:), '-o', 'Linewidth', 2 );
      title('vx');

      subplot( 2, 2, 4 );
      plot( Tf*nodes, X(4,:), '-o', 'Linewidth', 2 );
      title('vy');

    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, P )
      %L = self.lagrange_zero(nseg, tL, tR, XL, XR, UC, P);
      L = self.epsi1*UC(1)^2;
    end
    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, P )
      %gradL = self.lagrange_zero_gradient(nseg, tL, tR, XL, XR, UC, P);
      gradL = [0,0,0,0,0,0,0,0,self.epsi1*2*UC(1),0];
    end
    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, P )
      %hessL = self.lagrange_zero_hessian(nseg, tL, tR, XL, XR, UC, P);
      dim   = 2*self.nx+self.nu+self.np;
      hessL = sparse(dim,dim);
      hessL(9,9) = self.epsi1*2;
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      %hessL = self.lagrange_zero_hessian_pattern();
      dim   = 2*self.nx+self.nu+self.np;
      hessL = sparse(dim,dim);
      hessL(9,9) = 1;
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( self, tL, tR, XL, XR, TS )
      M = -XR(1);
    end
    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, TS )
      gradM = [ 0,0,0,0, -1, 0,0,0, 0 ];
    end
    %
    function hessM = mayer_hessian( self, tL, tR, XL, XR, TS )
      hessM = sparse( 2*self.nx + self.np, 2*self.nx + self.np );
    end
    %
    function hessM = mayer_hessian_pattern( self )
      hessM = sparse( 2*self.nx + self.np, 2*self.nx + self.np );
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, nseg, tL, tR, XL, XR, UC, PARS )
      C = self.midpoint_ds( nseg, tL, tR, XL, XR, UC, PARS, @self.RHS );
    end
    %
    function CJ = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, PARS )
      CJ = self.midpoint_ds_jacobian( nseg, tL, tR, XL, XR, UC, PARS, @self.JAC );
    end
    %
    function CJ = ds_jacobian_pattern( self )
      CJ = self.midpoint_ds_jacobian_pattern( @self.JAC_pattern );
    end
    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, PARS, L )
      H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, PARS, L );
    end
    %
    function H = ds_hessian_pattern( self )
      dim = 2*self.nx+self.nu+self.np;
      H = sparse(ones(dim,dim));
    end
    %              _   _                           _             _           _
    %  _ __   __ _| |_| |__     ___ ___  _ __  ___| |_ _ __ __ _(_)_ __  ___| |_
    % | '_ \ / _` | __| '_ \   / __/ _ \| '_ \/ __| __| '__/ _` | | '_ \/ __| __|
    % | |_) | (_| | |_| | | | | (_| (_) | | | \__ \ |_| | | (_| | | | | \__ \ |_
    % | .__/ \__,_|\__|_| |_|  \___\___/|_| |_|___/\__|_|  \__,_|_|_| |_|___/\__|
    % |_|
    %
    % Path constraints
    function C = pc( self, t, X, PARS )
      C = zeros(0,1);
    end
    %
    function CJ = pc_jacobian( self, t, X, PARS )
      CJ = zeros(0,self.nx);
    end
    %
    function CJ = pc_jacobian_pattern( self )
      CJ = sparse(0,self.nx);
    end
    %
    function CH = pc_hessian( self, t, X, PARS, L )
      CH = sparse(self.nx,self.nx);
    end
    %
    function CH = pc_hessian_pattern( self )
      CH = sparse(self.nx,self.nx);
    end

    %     _
    %  _ | |_  _ _ __  _ __
    % | || | || | '  \| '_ \
    %  \__/ \_,_|_|_|_| .__/
    %                 |_|
    %
    function ODE = jump( self, nsegL, t, XL, XR, P )
      ODE = self.jump_standard( nsegL, t, XL, XR, P );
    end
    %
    function JAC = jump_jacobian( self, nsegL, t, XL, XR, P )
      JAC = self.jump_standard_jacobian( nsegL, t, XL, XR, P );
    end
    %
    function JAC = jump_jacobian_pattern( self )
      JAC = self.jump_standard_jacobian_pattern();
    end
    %
    function H = jump_hessian( self, nsegL, t, XL, XR, P, L )
      H = self.jump_standard_hessian( nsegL, t, XL, XR, P, L );
    end
    %
    function H = jump_hessian_pattern( self )
      H = self.jump_standard_hessian_pattern();
    end

    %  ___  ___
    % | _ )/ __|
    % | _ \ (__
    % |___/\___|
    %
    function bc = bc( self, tL, tR, XL, XR, P )
      bc = [             ...
        XL(1);           ...
        XL(2)-self.y_i;  ...
        XR(2)-self.y_f;  ...
        XL(3)-self.vx_i; ...
        XR(3)-self.vx_f; ...
        XL(4)-self.vy_i; ...
        XR(4)-self.vy_f  ...
      ];
    end
    %
    function Jac = bc_jacobian( self, tL, tR, XL, XR, P )
      Jac = sparse( self.nbc, 2*self.nx + self.np );
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,6) = 1;
      Jac(4,3) = 1;
      Jac(5,7) = 1;
      Jac(6,4) = 1;
      Jac(7,8) = 1;
    end
    %
    function Jac = bc_jacobian_pattern( self )
      Jac = sparse( self.nbc, 2*self.nx + self.np );
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,6) = 1;
      Jac(4,3) = 1;
      Jac(5,7) = 1;
      Jac(6,4) = 1;
      Jac(7,8) = 1;
    end
    %
    function Hess = bc_hessian( self, tL, tR, XL, XR, P, L )
      Hess = sparse( 2*self.nx + self.np, 2*self.nx + self.np );
    end
    %
    function Hess = bc_hessian_pattern( self )
      Hess = sparse( 2*self.nx + self.np, 2*self.nx + self.np );
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    % tvU
    function tvU = TVU( self, dt, UCL, UCR )
      tvU = self.TVU_standard( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_standard_gradient( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_standard_hessian( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvH = TVU_hessian_pattern( self )
      tvH = self.TVU_standard_hessian_pattern();
    end

  end
end
