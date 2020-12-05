%--------------------------------------------------------------------------%
%                                                                          %
%  Copyright (C) 2018                                                      %
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

classdef OCP_GoddardRocket < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    h_i
    v_i
    v_f
    m_i
    m_f
    %
    Dc
    Tmax
    c
    g0
    %
    hc
    mc
    vc
    %
    epsi
  end

  methods ( Hidden = true )

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %
    function res = D( self, h, v )
      Dc  = self.Dc;
      hc  = self.hc;
      h_i = self.h_i;
      rr  = hc/h_i;
      ee  = exp( rr * (h_i - h) );
      res = Dc * (v^2) * ee;
    end
    %
    function res = D_h( self, h, v )
      Dc   = self.Dc;
      hc   = self.hc;
      h_i  = self.h_i;
      rr   = hc/h_i;
      ee_1 = -rr * exp( rr * (h_i - h) );
      res  = Dc * (v^2) * ee_1;
    end
    %
    function res = D_v( self, h, v )
      Dc  = self.Dc;
      hc  = self.hc;
      h_i = self.h_i;
      rr  = hc/h_i;
      ee  = exp( rr * (h_i - h) );
      res = 2 * Dc * v * ee;
    end
    %
    function res = D_hh( self, h, v )
      Dc   = self.Dc;
      hc   = self.hc;
      h_i  = self.h_i;
      rr   = hc/h_i;
      ee_2 = rr^2 * exp( rr * (h_i - h) );
      res  = Dc * (v^2) * ee_2;
    end
    %
    function res = D_hv( self, h, v )
      Dc   = self.Dc;
      hc   = self.hc;
      h_i  = self.h_i;
      rr   = hc/h_i;
      ee_1 = -rr * exp( rr * (h_i - h) );
      res  = 2 * Dc * v * ee_1;
    end
    %
    function res = D_vv( self, h, v )
      Dc  = self.Dc;
      hc  = self.hc;
      h_i = self.h_i;
      rr  = hc/h_i;
      ee  = exp( rr * (h_i - h) );
      res = 2 * Dc * ee;
    end
    %
    function res = g( self, h )
      g0  = self.g0;
      h_i = self.h_i;
      res = g0 * (h_i/h)^2;
    end
    %
    function res = g_D( self, h )
      g0  = self.g0;
      h_i = self.h_i;
      res =  -2 * g0 * h_i ^ 2 / h ^ 3;
    end
    %
    function res = g_DD( self, h )
      g0  = self.g0;
      h_i = self.h_i;
      res = 6 * g0 * h_i ^ 2 / h ^ 4;
    end
    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, P )
      %
      tf  = P(1);
      c   = self.c;
      h   = X(1);
      v   = X(2);
      m   = X(3);
      T   = U(1);
      D   = self.D(h,v);
      gh  = self.g(h);
      RES = tf*[ v; (T-D)/m-gh; -T/c ];
    end
    %
    function RES = JAC( self, nseg, t, X, U, P )
      RES = sparse( 3, 5 );
      %
      tf   = P(1);
      c    = self.c;
      h    = X(1);
      v    = X(2);
      m    = X(3);
      T    = U(1);
      D    = self.D(h,v);
      Dh   = self.D_h(h,v);
      Dv   = self.D_v(h,v);
      gh   = self.g(h);
      gh_D = self.g_D(h);

      RES(1,2) = tf;
      RES(1,5) = v;

      RES(2,1) = -tf*(Dh/m+gh_D);
      RES(2,2) = -tf*Dv/m;
      RES(2,3) = tf*(D-T)/m^2;
      RES(2,4) = tf/m;
      RES(2,5) = (T-D)/m-gh;

      RES(3,4) = -tf/c;
      RES(3,5) = -T/c;
    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( 3, 5 );
      RES(1,2) = 1;
      RES(1,5) = 1;
      RES(2,1) = 1;
      RES(2,2) = 1;
      RES(2,3) = 1;
      RES(2,4) = 1;
      RES(2,5) = 1;
      RES(3,4) = 1;
      RES(3,5) = 1;
    end
    %
    function H = HESS( self, nseg, t, X, U, P, L )
      H1    = sparse(5,5);
      H2    = sparse(5,5);
      H3    = sparse(5,5);

      tf    = P(1);
      c     = self.c;
      h     = X(1);
      v     = X(2);
      m     = X(3);
      T     = U(1);
      D     = self.D(h,v);
      Dh    = self.D_h(h,v);
      Dv    = self.D_v(h,v);
      Dhh   = self.D_hh(h,v);
      Dhv   = self.D_hv(h,v);
      Dvv   = self.D_vv(h,v);
      gh    = self.g(h);
      gh_D  = self.g_D(h);
      gh_DD = self.g_DD(h);

      % tf*v
      H1(2,5) = 1;
      H1(5,2) = 1;

      % tf*T/m - tf*D/m - tf*gh;

      % -Dh*tf/m - gh_D
      H2(1,1) = -tf*(Dhh/m+gh_DD);
      H2(1,2) = -tf*Dhv/m;
      H2(1,3) =  tf*Dh/m^2;
      H2(1,5) = -(Dh/m+gh_D);
      H2(2,1) = H2(1,2);
      H2(3,1) = H2(1,3);
      H2(5,1) = H2(1,5);

      % -tf*Dv/m
      H2(2,2) = -tf*Dvv/m;
      H2(2,3) =  tf*Dv/m^2;
      H2(2,5) = -Dv/m;
      H2(3,2) = H2(2,3);
      H2(5,2) = H2(2,5);

      % tf*(D-T)/m^2
      H2(3,3) = 2*tf*(T-D)/m^3;
      H2(3,4) = -tf/m^2;
      H2(3,5) = (D-T)/m^2;
      H2(4,3) = H2(3,4);
      H2(5,3) = H2(3,5);

      % tf/m;
      H2(4,5) = 1/m;
      H2(5,4) = H2(4,5);

      % -tf*T/c;
      H3(4,5) = -1/c;
      H3(5,4) = H3(4,5);

      H = L(1)*H1+L(2)*H2+L(3)*H3;
    end
    %
    function H = HESS_pattern( self )
      H = sparse(5,5);

      H(1,1) = 1;
      H(1,2) = 1;
      H(1,3) = 1;
      H(1,5) = 1;

      H(2,1) = 1;
      H(2,2) = 1;
      H(2,3) = 1;
      H(2,5) = 1;

      H(3,1) = 1;
      H(3,2) = 1;
      H(3,3) = 1;
      H(3,4) = 1;
      H(3,5) = 1;

      H(4,3) = 1;
      H(4,4) = 1;
      H(4,5) = 1;

      H(5,1) = 1;
      H(5,2) = 1;
      H(5,3) = 1;
      H(5,4) = 1;
    end
  end

  methods

    function self = OCP_GoddardRocket( )
      nx  = 3; % number of states
      nu  = 1; % number of controls
      np  = 1; % number of free parameters
      nbc = 5; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );

      self.hc   = 500;
      self.mc   = 0.6;
      self.vc   = 620;
      self.g0   = 1;

      self.h_i  = 1;
      self.v_i  = 0;
      self.v_f  = 0;
      self.m_i  = 1;
      self.m_f  = self.m_i * self.mc;

      self.Dc   = 0.5*self.vc*self.m_i/self.g0;
      self.Tmax = 3.5*self.g0*self.m_i;
      self.c    = 0.5*sqrt(self.g0*self.h_i);
      
      self.epsi = 1e-8;

    end

    function info = solve( self )

      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      N    = self.N;
      nseg = self.nseg;
      njmp = self.njmp;
      nbc  = self.nbc;
      Tmax = 10;

      utot = (N-nseg)*nu;
      xtot = N*nx;

      x_lb = [ zeros(1,N); zeros(1,N); zeros(1,N)];
      x_ub = [ Inf*ones(1,N); Inf*ones(1,N); self.m_i*ones(1,N)];
      u_lb = zeros(1, N-nseg);
      u_ub = self.Tmax*ones(1, N-nseg);

      options.lb = self.pack( x_lb, u_lb, 0   ); % Lower bound on the variables.
      options.ub = self.pack( x_ub, u_ub, 100 ); % Upper bound on the variables.

      % The constraint functions are bounded to zero
      dim = (N-nseg)*nx + (nseg-1)*njmp + nbc;
      options.cl = zeros(1,dim); %  constraints
      options.cu = zeros(1,dim);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant   = 'no';
      options.ipopt.hessian_constant = 'no';
      options.ipopt.mu_strategy      = 'adaptive';
      options.ipopt.max_iter         = 400;
      options.ipopt.tol              = 1e-10;%
      options.ipopt.linear_solver    = 'mumps';

      % The callback functions.
      funcs.objective         = @(Z) self.NLP_target(Z);
      funcs.gradient          = @(Z) self.NLP_target_gradient(Z);

      funcs.constraints       = @(Z) self.NLP_constraints(Z);
      funcs.jacobian          = @(Z) self.NLP_constraints_jacobian(Z);
      funcs.jacobianstructure = @() self.NLP_constraints_jacobian_pattern();

      if true
        %options.ipopt.derivative_test = 'first-order';
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian           = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure  = @() self.NLP_hessian_pattern();
      else
        options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      t_f_guess = 1;
      hguess    = ones(1,N) * self.h_i;
      vguess    = (self.nodes.*(1-self.nodes));
      mguess    = (self.m_f*self.nodes+self.m_i*(1-self.nodes));
      Tguess    = ones(1,utot)*self.Tmax/2;

      x0 = self.pack( [hguess;vguess;mguess], Tguess, t_f_guess ); % Lower bound on the variables.

      tic
      [self.sol, info] = ipopt(x0,funcs,options);

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();

      subplot( 2, 2, 1 );
      plot( nodes, X(1,:), '-o', 'Linewidth', 2 );
      title('h');

      subplot( 2, 2, 2 );
      plot( nodes, X(2,:), '-o', 'Linewidth', 2 );
      title('v');

      subplot( 2, 2, 3 );
      plot( nodes, X(3,:), '-o', 'Linewidth', 2 );
      title('m');

      subplot( 2, 2, 4 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Unodes, UU(1,:), '-r', 'Linewidth', 2 );
      title('thrust');

    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, P )
      L = self.lagrange_zero(nseg, tL, tR, XL, XR, UC, P);
    end
    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, P )
      gradL = self.lagrange_zero_gradient(nseg, tL, tR, XL, XR, UC, P);
    end
    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, P )
      hessL = self.lagrange_zero_hessian(nseg, tL, tR, XL, XR, UC, P);
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      hessL = self.lagrange_zero_hessian_pattern();
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( ~, tL, tR, XL, XR, TS )
      hL = XL(1); vL = XL(2); mL = XL(3);
      hR = XR(1); vR = XR(2); mR = XR(3);
      M  = -hR;
    end
    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, TS )
      gradM = [ 0, 0, 0, -1, 0, 0, 0];
    end
    %
    function hessM = mayer_hessian( self, tL, tR, XL, XR, TS )
      dim   = 2*self.nx + self.np;
      hessM = sparse(dim,dim);
    end
    %
    function hessM = mayer_hessian_pattern( self )
      dim   = 2*self.nx + self.np;
      hessM = sparse(dim,dim);
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
      H = self.midpoint_ds_hessian( nseg, tL, tR, XL, XR, UC, PARS, L, @self.HESS );
    end
    %
    function H = ds_hessian_pattern( self )
      H = self.midpoint_ds_hessian_pattern( @self.HESS_pattern );
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
      hL = XL(1); vL = XL(2); mL = XL(3);
      hR = XR(1); vR = XR(2); mR = XR(3);
      bc = [ hL - self.h_i; ...
             vL - self.v_i; ...
             mL - self.m_i; ...
             vR - self.v_f; ...
             mR - self.m_f ];
    end
    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, P )
      Jac = sparse( 5, 7 );
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,5) = 1;
      Jac(5,6) = 1;
    end
    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse( 5, 7 );
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,5) = 1;
      Jac(5,6) = 1;
    end
    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, P, L )
      Hess = zeros(7,7);
    end
    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = zeros(7,7);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    % tvU
    function tvU = TVU( self, dt, UCL, UCR )
      %tvU = self.TVU_zero( dt, UCL, UCR );
      tvU = self.TVU_reg( self.epsi, dt, UCL, UCR );
    end

    function tvG = TVU_gradient( self, dt, UCL, UCR )
      %tvG = self.TVU_zero_gradient( dt, UCL, UCR );
      tvG = self.TVU_reg_gradient( self.epsi, dt, UCL, UCR );
    end

    function tvH = TVU_hessian( self, dt, UCL, UCR )
      %tvH = self.TVU_zero_hessian( dt, UCL, UCR );
      tvH = self.TVU_reg_hessian( self.epsi, dt, UCL, UCR );
    end

    function tvH = TVU_hessian_pattern( self )
      %tvH = self.TVU_zero_hessian_pattern();
      tvH = self.TVU_reg_hessian_pattern();
    end

  end
end
