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

classdef OCP_NAME_OF_THE_CLASS < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    % PARAMETERS OF THE PROBLEM
    %
    epsi
  end

  methods ( Hidden = true )

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %


    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lag( self, nseg, t, X, U, P )
      ...
    end
    %
    function gradL = lag_gradient( self, nseg, t, X, U, P )
      ...
    end
    %
    function hessL = lag_hessian( self, nseg, t, X, U, P )
      ...
    end
    %
    function hessL = lag_hessian_pattern( self )
      ...
    end

    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, P )
      % rhs of teh ODE x' = f(x,u,p,t)
      ...
      ...
    end
    %
    function RES = JAC( self, nseg, t, X, U, P )
      RES = sparse( self.nx, self.nx+self.nu+self.np );
      ...
      ...
    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( self.nx, self.nx+self.nu+self.np );
      ...
      ...
    end
    %
    function H = HESS( self, nseg, t, X, U, P, L )
      H = sparse( self.nx+self.nu+self.np, self.nx+self.nu+self.np );
     ...
     ...
    end
    %
    function H = HESS_pattern( self )
      H = sparse( self.nx+self.nu+self.np, self.nx+self.nu+self.np );
      ...
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
      % setup the parameters of the model

      self.epsi  = 0; % no penalization of total variations of the control
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
      Tmax = self.Tmax;

      utot = (N-nseg)*nu;
      xtot = N*nx;

      % set the lower bound of the states (nx x N)
      x_lb = [ zeros(1,N); ...
               zeros(1,N); ...
               zeros(1,N)];

      % set the upper bound of the states (nx x N)
      x_ub = [ Inf*ones(1,N); ...
               Inf*ones(1,N); ...
               Inf*ones(1,N)];

      % set the lower bound of the controls (nu x (N-nseg))
      u_lb = zeros(1, N-nseg);

      % set the upper bound of the controls (nu x (N-nseg))
      u_ub = ones(1, N-nseg);

      % set the lower bound of the optimization parameters (1 x np)
      p_lb = [];

      % set the upper bound of the optimization parameters (1 x np)
      p_ub = [];

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
      options.ipopt.max_iter         = 1000;
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
        %options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % GUESS SOLUTION

      x0 = self.pack( Xguess, Uguess, Pguess ); % Lower bound on the variables.

      tic
      [self.sol, info] = ipopt(x0,funcs,options);

      self.sol(end)

    end

    function plot( self )
      % example of pliot solution
      nodes = self.nodes;
      X     = self.states();
      Tf    = self.sol(end);

      subplot( 2, 2, 1 );
      plot( Tf*nodes, X(1,:), '-o', 'Linewidth', 2 );
      title('h');

      subplot( 2, 2, 2 );
      plot( Tf*nodes, X(2,:), '-o', 'Linewidth', 2 );
      title('v');

      subplot( 2, 2, 3 );
      plot( Tf*nodes, X(3,:), '-o', 'Linewidth', 2 );
      title('m');

      subplot( 2, 2, 4 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Tf*Unodes, UU(1,:), '-r', 'Linewidth', 2 );
      title('thrust');

    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, P )
      L = self.midpoint_lagrange( nseg, tL, tR, XL, XR, UC, PARS, @self.lag );
    end
    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, P )
      gradL = self.midpoint_lagrange_gradient( nseg, tL, tR, XL, XR, UC, PARS, @self.lag_gradient );
    end
    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, P )
      hessL = self.midpoint_lagrange_hessian( nseg, tL, tR, XL, XR, UC, PARS, @self.lag_hessian );
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      hessL = self.midpoint_lagrange_hessian_pattern( nseg, tL, tR, XL, XR, UC, PARS, @self.lag_hessian_pattern );
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( self, tL, tR, XL, XR, TS )
    end
    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, TS )
    end
    %
    function hessM = mayer_hessian( self, tL, tR, XL, XR, TS )
    end
    %
    function hessM = mayer_hessian_pattern( self )
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
      ....
      ....
    end
    %
    function Jac = bc_jacobian( self, tL, tR, XL, XR, P )
      ...
      ...
    end
    %
    function Jac = bc_jacobian_pattern( self )
      ...
      ...
    end
    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, P, L )
      ...
      ...
    end
    %
    function Hess = bc_hessian_pattern( ~ )
      ...
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    % tvU
    function tvU = TVU( self, dt, UCL, UCR )
      tvU = self.TVU_reg( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_reg_gradient( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_reg_hessian( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvH = TVU_hessian_pattern( self )
      tvH = self.TVU_reg_hessian_pattern();
    end

  end
end
