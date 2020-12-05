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

classdef OCP_BangBang < OCP_NLP

  properties ( SetAccess = private, Hidden = true )
    x_i
    v_i
    x_f
    v_f
  end

  methods ( Hidden = true )
    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, T )
      RES = [ T*X(2); T*U ];
    end
    %
    function RES = JAC( self, nseg, t, X, U, T )
      RES = sparse( 2, 4 );
      RES(1,2) = T;
      RES(1,4) = X(2);
      RES(2,3) = T;
      RES(2,4) = U;
    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( 2, 4 );
      RES(1,2) = 1;
      RES(1,4) = 1;
      RES(2,3) = 1;
      RES(2,4) = 1;
    end
    %
    function H = HESS( self, nseg, t, X, U, T, L )
      H = sparse(4,4);
      H(2,4) = L(1); H(4,2) = L(1);
      H(3,4) = L(2); H(4,3) = L(2);
    end
    %
    function H = HESS_pattern( self )
      H = sparse(4,4);
      H(2,4) = 1; H(4,2) = 1;
      H(3,4) = 1; H(4,3) = 1;
    end
  end

  methods

    function self = OCP_BangBang( )
      nx  = 2; % number of states
      nu  = 1; % number of controls
      np  = 1; % number of free parameters
      nbc = 4; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.x_i = 0;
      self.x_f = 1;
      self.v_i = 0;
      self.v_f = 0;
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

      options.lb = self.pack( -Inf*ones( nx, N ), -ones( nu, N-nseg ), 0    ); % Lower bound on the variables.
      options.ub = self.pack(  Inf*ones( nx, N ),  ones( nu, N-nseg ), Tmax ); % Upper bound on the variables.

      % The constraint functions are bounded to zero
      nc = (N-nseg)*nx+(nseg-1)*njmp+nbc;
      options.cl = zeros(1,nc); % constraints
      options.cu = zeros(1,nc);

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

      if true
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian           = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure  = @() self.NLP_hessian_pattern();
      else
        %options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      xguess = (self.x_i+(self.x_f-self.x_i)*self.nodes(:)).';
      vguess = zeros(1,N);
      uguess = zeros(1,N-nseg);

      x0 = self.pack( [ xguess; vguess], uguess, 1 );

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      T     = self.parameters();

      subplot( 3, 1, 1 );
      plot( nodes, X(1,:), 'Linewidth', 2 );
      title('x');

      subplot( 3, 1, 2 );
      plot( nodes, X(2,:), 'Linewidth', 2 );
      title('v');

      subplot( 3, 1, 3 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Unodes, UU(1,:), 'Linewidth', 2 );
      title('u');

    end

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
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, P )
      L = self.lagrange_zero( nseg, tL, tR, XL, XR, UC, P );
    end

    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, P )
      gradL = self.lagrange_zero_gradient( nseg, tL, tR, XL, XR, UC, P );
    end

    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, P )
      hessL = self.lagrange_zero_hessian( nseg, tL, tR, XL, XR, UC, P );
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
    function M = mayer( ~, tL, tR, XL, XR, T )
      M = T; % only one parameter
    end

    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, T )
      dim   = 2*self.nx + self.np;
      gradM = zeros(1,dim);
      gradM(2*self.nx+1) = 1;
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( self, tL, tR, XL, XR, T )
      dim   = 2*self.nx + self.np;
      hessM = sparse(dim,dim);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( self )
      dim   = 2*self.nx + self.np;
      hessM = sparse(dim,dim);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, nseg, tL, tR, XL, XR, UC, T )
      C = self.midpoint_ds( nseg, tL, tR, XL, XR, UC, T, @self.RHS );
    end

    %
    function CJ = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, T )
      CJ = self.midpoint_ds_jacobian( nseg, tL, tR, XL, XR, UC, T, @self.JAC );
    end

    %
    function CJ = ds_jacobian_pattern( self )
      CJ = self.midpoint_ds_jacobian_pattern( @self.JAC_pattern );
    end

    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, T, L )
      H = self.midpoint_ds_hessian( nseg, tL, tR, XL, XR, UC, T, L, @self.HESS );
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
      xL = XL(1); vL = XL(2);
      xR = XR(1); vR = XR(2);
      bc = [ xL - self.x_i; ...
             xR - self.x_f; ...
             vL - self.v_i; ...
             vR - self.v_f ];
    end

    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, P )
      Jac      = sparse(4,5);
      Jac(1,1) = 1;
      Jac(2,3) = 1;
      Jac(3,2) = 1;
      Jac(4,4) = 1;
    end

    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac      = sparse(4,5);
      Jac(1,1) = 1;
      Jac(2,3) = 1;
      Jac(3,2) = 1;
      Jac(4,4) = 1;
    end

    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, P, L )
      Hess = sparse(5,5);
    end

    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = sparse(5,5);
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

    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_zero_gradient( dt, UCL, UCR );
    end

    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_zero_hessian( dt, UCL, UCR );
    end

    function tvH = TVU_hessian_pattern( self, dt, UCL, UCR )
      tvH = self.TVU_zero_hessian_pattern();
    end

  end
end
