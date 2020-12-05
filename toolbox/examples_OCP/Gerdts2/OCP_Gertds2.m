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

classdef OCP_Gertds2 < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x1_i
    x2_i
    x3_i
    x1_f
    x2_f
  end

  methods

    function self = OCP_Gertds2( )
      nx  = 3; % number of states
      nu  = 1; % number of controls
      np  = 0; % number of free parameters
      nbc = 5; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.x1_i = 0;
      self.x2_i = 1;
      self.x3_i = 0;
      self.x1_f = 0;
      self.x2_f = -1;
    end

    function info = solve( self )

      xones = ones(1,self.N*self.nx);
      uones = ones(1,(self.N-self.nseg)*self.nu);

      options.lb = [ -xones*Inf, -uones*Inf ];  % Lower bound on the variables.
      options.ub = [ reshape( [1/9;Inf;Inf] * ones(1,self.N), 1, self.N*self.nx ), ...
                     ones(1,self.N-self.nseg)*Inf ];  % Upper bound on the variables.

      % The constraint functions are bounded to zero
      nc = (self.N-self.nseg)*self.nx+(self.nseg-1)*self.njmp+self.nbc;
      options.cl = zeros(1,nc); %  constraints
      options.cu = zeros(1,nc);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant   = 'no';
      options.ipopt.hessian_constant = 'no';
      options.ipopt.mu_strategy      = 'adaptive';
      options.ipopt.max_iter         = 400;
      options.ipopt.tol              = 1e-10;
      options.ipopt.linear_solver    = 'mumps';

      % The callback functions.
      funcs.objective         = @(Z) self.NLP_target(Z);
      funcs.gradient          = @(Z) self.NLP_target_gradient(Z);

      funcs.constraints       = @(Z) self.NLP_constraints(Z);
      funcs.jacobian          = @(Z) self.NLP_constraints_jacobian(Z);
      funcs.jacobianstructure = @() self.NLP_constraints_jacobian_pattern();

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
      x0 = zeros(self.nx*self.N+self.nu*(self.N-self.nseg),1);

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      T     = self.parameters();

      subplot( 2, 1, 1 );
      plot( nodes, X(1,:), nodes, X(2,:), nodes, X(3,:), 'Linewidth', 2 );
      title('x');

      subplot( 2, 1, 2 );
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
    function L = lagrange( ~, ~, tL, tR, XL, XR, UC, ~ )
      L = 0;
    end
    %
    function gradL = lagrange_gradient( self, ~, tL, tR, XL, XR, UC, ~ )
      gradL = zeros(1,2*self.nx+self.nu);
    end
    %
    function hessL = lagrange_hessian( self, ~, tL, tR, XL, XR, UC, ~ )
      dim = 2*self.nx+self.nu;
      hessL = sparse( dim, dim );
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      dim = 2*self.nx+self.nu;
      hessL = sparse( dim, dim );
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( ~, tL, tR, XL, XR, ~ )
      M = XR(3);
    end

    %
    function gradM = mayer_gradient( ~, tL, tR, XL, XR, ~ )
      gradM = [0,0,0,0,0,1];
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( ~, tL, tR, XL, XR, ~ )
      hessM = zeros(6,6);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( ~ )
      hessM = zeros(6,6);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, ~, tL, tR, XL, XR, UC, ~ )
      % ----------
      DT  = tR-tL;
      % ----------
      C    = zeros(3,1);
      C(1) = (XR(1) - XL(1))/DT - (XL(2)+XR(2))/2;
      C(2) = (XR(2) - XL(2))/DT - UC(1);
      C(3) = (XR(3) - XL(3))/DT - 0.5*UC(1)^2;
    end

    %
    function JAC = ds_jacobian( self, ~, tL, tR, XL, XR, UC, ~ )
      % ----------
      DT  = tR-tL;
      % ----------
      JAC = sparse( [ -1/DT,  -0.5,     0, 1/DT, -0.5,    0,  0; ...
                          0, -1/DT,     0,    0, 1/DT,    0, -1; ...
                          0,     0, -1/DT,    0,    0, 1/DT, -UC(1) ] );
    end

    %
    function JAC = ds_jacobian_pattern( self )
      JAC = sparse( [ 1, 1, 0, 1, 1, 0, 0; ...
                      0, 1, 0, 0, 1, 0, 1; ...
                      0, 0, 1, 0, 0, 1, 1 ] );
    end

    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, P, L )
      if false
        H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, P, L );
      else
        H      = sparse(7,7);
        H(7,7) = -L(3);
      end
    end

    %
    function H = ds_hessian_pattern( self )
      H      = sparse(7,7);
      H(7,7) = 1;
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
      CJ = zeros(0,self.nx);
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
    function bc = bc( self, tL, tR, XL, XR, ~ )
      bc = [ XL(1) - self.x1_i; ...
             XL(2) - self.x2_i; ...
             XL(3) - self.x3_i; ...
             XR(1) - self.x1_f; ...
             XR(2) - self.x2_f ];
    end
    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, ~ )
      Jac = sparse(5,6);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,4) = 1;
      Jac(5,5) = 1;
    end
    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse(5,6);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,4) = 1;
      Jac(5,5) = 1;
    end
    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, L )
      Hess = zeros(6,6);
    end
    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = zeros(6,6);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    %
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
