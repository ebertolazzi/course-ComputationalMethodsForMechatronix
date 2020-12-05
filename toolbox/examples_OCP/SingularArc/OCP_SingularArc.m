%-------------------------------------------------------------------------%
%                                                                         %
%  Copyright (C) 2018                                                     %
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

classdef OCP_SingularArc < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x1_i
    x2_i
    x2_f
    x3_i
    x3_f
    Tmin
    Tmax
  end

  methods

    function self = OCP_SingularArc( )
      nx  = 3; % number of states
      nu  = 1; % number of controls
      np  = 1; % number of free parameters
      nbc = 5; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.x1_i = pi/2;
      self.x2_i = 4;
      self.x3_i = 0;
      self.x2_f = 0;
      self.x3_f = 0;
      self.Tmin = 1;
      self.Tmax = 100;
    end

    function info = solve( self )

      utot = (self.N-self.nseg)*self.nu;
      xtot = self.N*self.nx;

      xones      = ones(1,xtot);
      uones      = ones(1,utot);
      options.lb = [ -10*xones, -2*uones, self.Tmin ]; % Lower bound on the variables.
      options.ub = [  10*xones,  2*uones, self.Tmax ]; % Upper bound on the variables.

      % The constraint functions are bounded to zero
      dim = (self.N-self.nseg)*self.nx + (self.nseg-1)*self.njmp + self.nbc;
      options.cl = zeros(1,dim); %  constraints
      options.cu = zeros(1,dim);

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
      funcs.jacobianstructure = @()  self.NLP_constraints_jacobian_pattern();

      if true
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian                 = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure        = @() self.NLP_hessian_pattern();
      else
        %options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      ss      = (self.nodes-self.nodes(1))./(self.nodes(end)-self.nodes(1));
      x1guess = self.x1_i*ones(size(ss));
      x2guess = self.x2_i+(self.x2_f-self.x2_i)*ss;
      x3guess = self.x3_i*ones(size(ss));%+(self.x2_f-self.x3_i)*ss;
      uguess  = 0*ones(1,self.N-self.nseg);

      x0 = [ reshape( [ x1guess, x2guess, x3guess], self.N*self.nx ,1 ); ...
             reshape( uguess, (self.N-self.nseg)*self.nu ,1 ); ...
             50 ];

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      %T     = self.parameters();

      subplot( 2, 1, 1 );
      plot( nodes, X(1,:), nodes, X(2,:), nodes, X(3,:),'Linewidth', 2 );
      title('x1,x2,x3');

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
    function L = lagrange( self, ~, tL, tR, XL, XR, UC, TF )
      L = 0;
    end

    %
    function gradL = lagrange_gradient( self, ~, tL, tR, XL, XR, UC, TF )
      gradL = [ 0, 0, 0, 0, 0, 0, 0, 0 ];
    end

    %
    function hessL = lagrange_hessian( ~, ~, tL, tR, XL, XR, UC, ~ )
      hessL = sparse(8,8);
    end

    %
    function hessL = lagrange_hessian_pattern( ~ )
      hessL = sparse(8,8);
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( ~, tL, tR, XL, XR, TF )
      M = TF;
    end

    %
    function gradM = mayer_gradient( ~, tL, tR, XL, XR, TF )
      gradM = [ 0, 0, 0, 0, 0, 0, 1 ];
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( ~, tL, tR, XL, XR, ~ )
      hessM = sparse(7,7);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( ~ )
      hessM = sparse(7,7);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, ~, tL, tR, XL, XR, UC, TF )
      x1L = XL(1); x2L = XL(2); x3L = XL(3);
      x1R = XR(1); x2R = XR(2); x3R = XR(3);
      x1m = (x1L+x1R)/2;
      x2m = (x2L+x2R)/2;
      x3m = (x3L+x3R)/2;
      u   = UC(1);
      % ----------
      DT  = tR - tL;
      % ----------
      C     = zeros(3,1);
      C(1)  = (x1R - x1L)/DT - TF*u;
      C(2)  = (x2R - x2L)/DT - TF*cos(x1m);
      C(3)  = (x3R - x3L)/DT - TF*sin(x1m);
    end

    %
    function JAC = ds_jacobian( self, ~, tL, tR, XL, XR, UC, TF )
      x1L = XL(1); x2L = XL(2); x3L = XL(3);
      x1R = XR(1); x2R = XR(2); x3R = XR(3);
      x1m = (x1L+x1R)/2;
      x2m = (x2L+x2R)/2;
      x3m = (x3L+x3R)/2;
      u   = UC(1);
      % ----------
      DT = tR - tL;
      JAC = sparse(3,8);

      JAC(1,1) = -1/DT;
      JAC(1,4) = 1/DT;
      JAC(1,7) = -TF;
      JAC(1,8) = -u;

      JAC(2,2) = -1/DT;
      JAC(2,5) = 1/DT;
      JAC(2,1) = 0.5*TF*sin(x1m);
      JAC(2,4) = 0.5*TF*sin(x1m);
      JAC(2,8) = -cos(x1m);

      JAC(3,3) = -1/DT;
      JAC(3,6) = 1/DT;
      JAC(3,1) = -0.5*TF*cos(x1m);
      JAC(3,4) = -0.5*TF*cos(x1m);
      JAC(3,8) = -sin(x1m);
    end

    %
    function JAC = ds_jacobian_pattern( self )
      JAC = sparse(3,8);

      JAC(1,1) = 1;
      JAC(1,4) = 1;
      JAC(1,7) = 1;
      JAC(1,8) = 1;

      JAC(2,2) = 1;
      JAC(2,5) = 1;
      JAC(2,1) = 1;
      JAC(2,4) = 1;
      JAC(2,8) = 1;

      JAC(3,3) = 1;
      JAC(3,6) = 1;
      JAC(3,1) = 1;
      JAC(3,4) = 1;
      JAC(3,8) = 1;
    end

    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, P, L )
      H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, P, L );
    end

    %
    function H = ds_hessian_pattern( self )
      H = sparse(ones(8,8));
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
    function bc = bc( self, tL, tR, XL, XR, ~ )
      x1L = XL(1); x2L = XL(2); x3L = XL(3);
      x1R = XR(1); x2R = XR(2); x3R = XR(3);
      bc = [ x1L - self.x1_i; ...
             x2L - self.x2_i; ...
             x3L - self.x3_i; ...
             x2R - self.x2_f; ...
             x3R - self.x3_f ];
    end

    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, ~ )
      Jac = sparse(5,7);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,5) = 1;
      Jac(5,6) = 1;
    end

    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse(5,7);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,5) = 1;
      Jac(5,6) = 1;
    end

    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, L )
      Hess = sparse(7,7);
    end

    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = sparse(7,7);
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
