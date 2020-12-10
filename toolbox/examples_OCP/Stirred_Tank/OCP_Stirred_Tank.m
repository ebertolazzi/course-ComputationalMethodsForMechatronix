%-------------------------------------------------------------------------%
%                                                                         %
%  Copyright (C) 2020                                                     %
%                                                                         %
%         , __                 , __                                       %
%        /|/  \               /|/  \                                      %
%         | __/ _   ,_         | __/ _   ,_                               %
%         |   \|/  /  |  |   | |   \|/  /  |  |   |                       %
%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                      %
%                           /|                   /|                       %
%                           \|                   \|                       %
%                                                                         %
%      Enrico Bertolazzi                                                  %
%      Dipartimento di Ingegneria Industriale                             %
%      Universita` degli Studi di Trento                                  %
%      email: enrico.bertolazzi@unitn.it                                  %
%                                                                         %
%-------------------------------------------------------------------------%

classdef OCP_Stirred_Tank < OCP_NLP
  % http://tomdyn.com/examples/stirredTank.html

  properties (SetAccess = private, Hidden = true)
    x1_i
    x1_f
    x2_i
    x2_f
    prb
  end

  methods ( Hidden = true )

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %

    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, P )
      %
      x1  = X(1);
      x2  = X(2);
      u   = U(1);
      a1  = x1+0.25;
      a2  = x2+0.5;
      a3  = x1+2;
      a4  = a2*exp(21*x1/a3);
      RES = [-2*a1+a4-a1*u;0.5-x2-a4];
    end
    %
    function RES = JAC( self, nseg, t, X, U, P )
      RES = sparse( 2, 3 );
      x1  = X(1);
      x2  = X(2);
      u   = U(1);
      ee  = exp(21*x1/(x1+2));
      RES(1,1) = -u-2+(42.*x2+21.0)*ee/(x1+2)^2;
      RES(1,2) = ee;
      RES(1,3) = -x1-0.25;
      RES(2,1) = -(42*(x2+0.5))*ee/(x1+2)^2;
      RES(2,2) = -1-ee;
      %
    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( 2, 3 );
      RES(1,1) = 1;
      RES(1,2) = 1;
      RES(1,3) = 1;
      RES(2,1) = 1;
      RES(2,2) = 1;
    end
    %
    function H = HESS( self, nseg, t, X, U, P, L )
      H1 = sparse(3,3);
      H2 = sparse(3,3);

      x1 = X(1);
      x2 = X(2);
      u  = U(1);
      ee = exp(21*x1/(x1+2));

      H1(1,1) = -(84*(x2+.5))*ee*(x1-19)/(x1+2)^4;
      H1(1,2) = 42*ee/(x1+2)^2;
      H1(1,3) = -1;
      
      H1(2,1) = H1(1,2);
      H1(3,1) = H1(1,3);
      
      H2(1,1) = (84*(x2+.5))*ee*(x1-19)/(x1+2)^4;
      H2(1,2) = -42*ee/(x1+2)^2;
      H2(2,1) = H2(1,2);

      H = L(1)*H1+L(2)*H2;
    end
    %
    function H = HESS_pattern( self )
      H = sparse(3,3);
      H(1,1) = 1;
      H(1,2) = 1;
      H(1,3) = 1;
      H(2,1) = 1;
      H(3,1) = 1;
    end
    
    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %

    function L = lag( self, nseg, t, X, U, P )
      x1 = X(1);
      x2 = X(2);
      u  = U(1);
      prb = self.prb;
      switch prb
      case {'a','b'}
        L  = (x1^2+x2^2+0.1*u^2)/2;
      case 'c'
        L  = (x1^2+x2^2)/2;
      end
    end
    %
    function gradL = lag_gradient( self, nseg, t, X, U, P )
      x1 = X(1);
      x2 = X(2);
      u  = U(1);
      prb = self.prb;
      switch prb
      case {'a','b'}
        gradL = [ x1, x2, 0.1*u ];
      case 'c'
        gradL = [ x1, x2, 0 ];
      end
    end
    %
    function hessL = lag_hessian( self, nseg, t, X, U, P )
      hessL = sparse(3,3);
      x1 = X(1);
      x2 = X(2);
      u  = U(1);
      hessL(1,1) = 1;
      hessL(2,2) = 1;
      prb = self.prb;
      switch prb
      case {'a','b'}
        hessL(3,3) = 0.1;
      case 'c'
      end
    end
    %
    function hessL = lag_hessian_pattern( self )
      hessL = sparse(3,3);
      hessL(1,1) = 1;
      hessL(2,2) = 1;
      prb = self.prb;
      switch prb
      case {'a','b'}
        hessL(3,3) = 1;
      case 'c'
      end
    end

  end

  methods

    function self = OCP_Stirred_Tank( )
      nx  = 2; % number of states
      nu  = 1; % number of controls
      np  = 0; % number of free parameters
      prb = 'a';
      switch prb
      case 'a'
        nbc = 2; % number of boundary conditions
      case {'b','c'}
        nbc = 4; % number of boundary conditions
      end    
      self@OCP_NLP( nx, nu, np, nbc );
      self.prb  = prb;
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.x1_i = 0.05;
      self.x2_i = 0;
      self.x1_f = 0;
      self.x2_f = 0;
    end

    function info = solve( self )

      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      N    = self.N;
      nseg = self.nseg;
      njmp = self.njmp;
      nbc  = self.nbc;

      utot = (N-nseg)*nu;
      xtot = N*nx;

      x_lb = [ -ones(1,N); -Inf*ones(1,N)];
      x_ub = [ Inf*ones(1,N); Inf*ones(1,N)];
      u_lb = -Inf*ones(1, N-nseg);

      switch self.prb
      case {'a','b'}
        u_ub = Inf*ones(1, N-nseg);
      case 'c'
        u_ub = ones(1, N-nseg);
      end

      options.lb = self.pack( x_lb, u_lb, [] ); % Lower bound on the variables.
      options.ub = self.pack( x_ub, u_ub, [] ); % Upper bound on the variables.

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
      funcs.jacobianstructure = @()  self.NLP_constraints_jacobian_pattern();

      if false
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian                 = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure        = @() self.NLP_hessian_pattern();
      else
        options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      x1guess = zeros(1,N);
      x2guess = zeros(1,N);
      uguess  = zeros(1,utot);

      x0 = self.pack( [x1guess;x2guess], uguess, [] ); % Lower bound on the variables.

      tic
      [self.sol, info] = ipopt(x0,funcs,options);

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      %T     = self.parameters();

      subplot( 2, 1, 1 );
      plot( nodes, X(1,:), nodes, X(2,:), 'Linewidth', 2 );
      title('x1,x2');

      subplot( 2, 1, 2 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Unodes, UU(1,:), 'Linewidth', 2 );
      title('u');

    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, PARS )
      L = self.midpoint_lagrange( nseg, tL, tR, XL, XR, UC, PARS, @self.lag );
    end
    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, PARS )
      gradL = self.midpoint_lagrange_gradient( nseg, tL, tR, XL, XR, UC, PARS, @self.lag_gradient );
    end
    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, PARS )
      hessL = self.midpoint_lagrange_hessian( nseg, tL, tR, XL, XR, UC, PARS, @self.lag_hessian );
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      hessL = self.midpoint_lagrange_hessian_pattern();
    end


    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( self, tL, tR, XL, XR, P )
      switch self.prb
      case 'a'
        x1 = XR(1);
        x2 = XR(2);
        M  = x1*x1 + x2*x2;
      case {'b','c'}
        M = 0;
      end
    end

    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, P )
      switch self.prb
      case 'a'
        x1 = XR(1);
        x2 = XR(2);
        gradM = [ 0, 0, 2*x1, 2*x2 ];
      case {'b','c'}
        gradM = [ 0, 0, 0, 0 ];
      end
     
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( self, tL, tR, XL, XR, ~ )
      hessM = sparse(4,4);
      hessM(3,3) = 1;
      hessM(4,4) = 1;
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( self )
      hessM = sparse(4,4);
      hessM(3,3) = 1;
      hessM(4,4) = 1;
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

    function CJ = pc_jacobian( self, t, X, PARS )
      CJ = sparse(0,self.nx);
    end

    function CJ = pc_jacobian_pattern( self )
      CJ = sparse(0,self.nx);
    end

    function CH = pc_hessian( self, t, X, PARS, L )
      CH = sparse(self.nx,self.nx);
    end

    function CH = pc_hessian_pattern( self )
      CH = sparse(self.nx,self.nx);
    end

    %     _
    %  _ | |_  _ _ __  _ __
    % | || | || | '  \| '_ \
    %  \__/ \_,_|_|_|_| .__/
    %                 |_|
    %
    function JMP = jump( self, nsegL, t, XL, XR, P )
      JMP = self.jump_standard( nsegL, t, XL, XR, P );
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
      H = self.jump_standard_hessian_pattern( );
    end

    %  ___  ___
    % | _ )/ __|
    % | _ \ (__
    % |___/\___|
    %
    function bc = bc( self, tL, tR, XL, XR, P )
      x1L = XL(1); x2L = XL(2);
      x1R = XR(1); x2R = XR(2);
      switch self.prb
      case 'a'
        bc = [ x1L - self.x1_i; ...
               x2L - self.x2_i ];
      case {'b','c'}
        bc = [ x1L - self.x1_i; ...
               x2L - self.x2_i; ...
               x1R - self.x1_f; ...
               x2R - self.x2_f ];
      end    
    end

    %
    function Jac = bc_jacobian( self, tL, tR, XL, XR, ~ )
      switch self.prb
      case 'a'
        Jac = sparse(2,4);
        Jac(1,1) = 1;
        Jac(2,2) = 1;
      case {'b','c'}
        Jac = sparse(4,4);
        Jac(1,1) = 1;
        Jac(2,2) = 1;
        Jac(3,3) = 1;
        Jac(4,4) = 1;
      end    
    end

    %
    function Jac = bc_jacobian_pattern( self )
      switch self.prb
      case 'a'
        Jac = sparse(2,4);
        Jac(1,1) = 1;
        Jac(2,2) = 1;
      case {'b','c'}
        Jac = sparse(4,4);
        Jac(1,1) = 1;
        Jac(2,2) = 1;
        Jac(3,3) = 1;
        Jac(4,4) = 1;
      end    
    end

    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, ~ )
      Hess = sparse(4,4);
    end

    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = sparse(4,4);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    % tvU
    function tvU = TVU( self, dt, UCL, UCR )
      tvU = self.TVU_zero( dt, UCL, UCR );
    end
    %
    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_zero_gradient( dt, UCL, UCR );
    end
    %
    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_zero_hessian( dt, UCL, UCR );
    end
    %
    function tvH = TVU_hessian_pattern( self )
      tvH = self.TVU_zero_hessian_pattern();
    end

  end
end
