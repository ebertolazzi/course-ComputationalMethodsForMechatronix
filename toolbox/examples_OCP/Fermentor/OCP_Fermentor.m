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


%
% test "Fermentor" from:
%
% Imperial College London Optimal Control Software User Guide (ICLOCS)
% Paola Falugi Eric Kerrigan Eugene van Wyk
% Department of Electrical and Electronic Engineering,
% Imperial College London London England, UK iclocs@imperial.ac.uk
% 6 May 2010
%
%
classdef OCP_Fermentor < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x1_i
    x2_i
    x3_i
    x4_i

    x1_max
    x2_max
    x3_max
    x4_max

    x1_min
    x2_min
    x3_min
    x4_min

    u_min
    u_max
  end

  methods

    function self = OCP_Fermentor( )
      nx  = 4; % number of states
      nu  = 1; % number of controls
      np  = 0; % number of free parameters
      nbc = 4; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.x1_i = 1.5;
      self.x2_i = 0;
      self.x3_i = 0;
      self.x4_i = 7;

      self.x1_max = 40;
      self.x2_max = 50;
      self.x3_max = 25;
      self.x4_max = 10;

      self.x1_min = 0;
      self.x2_min = 0;
      self.x3_min = 0;
      self.x4_min = 0;

      self.u_min = 0;
      self.u_max = 50;
    end

    function info = solve( self )

      N   = self.N;
      nx  = self.nx;
      nu  = self.nu;
      nbc = self.nbc;

      XLB = reshape( [ self.x1_min; self.x2_min; self.x3_min; self.x4_min ] * ones(1,N), 1, N*nx );
      XUB = reshape( [ self.x1_max; self.x2_max; self.x3_max; self.x4_max ] * ones(1,N), 1, N*nx );

      ULB = self.u_min * ones(1,N-self.nseg);
      UUB = self.u_max * ones(1,N-self.nseg);

      options.lb = [ XLB, ULB ];
      options.ub = [ XUB, UUB ];

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
      s  = (0:N-1)/(N-1);
      X1 = 1.5 + (30-1.5)*s;
      X2 = 0   + (8.5-0)*s;
      X3 = zeros(1,N);
      X4 = 7   + (10-7)*s;
      UU = 2   + (10-2)*s;

      x0 = [ reshape( [ X1; X2; X3; X4], 1, 4*N ), UU(1:end-self.nseg) ].'; 

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )

      nodes = self.nodes;
      X     = self.states();

      subplot( 2, 1, 1 );
      plot( nodes, X(1,:), nodes, X(2,:), nodes, X(3,:), nodes, X(4,:), 'Linewidth', 2 );
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
    function res = h1( self, x1, x3 )
      res = 0.11*(x3/(0.006*x1+x3));
    end

    function res = h1_1( self, x1, x3 )
      res = -0.66e-3*x3/(0.6e-2*x1+x3)^2;
    end

    function res = h1_3( self, x1, x3 )
      res = 0.66e-3*x1/(0.6e-2*x1+x3)^2;
    end

    function res = h1_11( self, x1, x3 )
      res = 0.792e-5*x3/(0.6e-2*x1+x3)^3;
    end

    function res = h1_13( self, x1, x3 )
      res = (-0.396e-5*x1+0.66e-3*x3)/(0.6e-2*x1+x3)^3;
    end

    function res = h1_33( self, x1, x3 )
      res = -0.132e-2*x1/(0.6e-2*x1+x3)^3;
    end

    function res = h2( self, x3 )
      res = 0.0055*(x3/(0.0001+x3*(1+10*x3)));
    end

    function res = h2_D( self, x3 )
      res = (5.5*10^(-7)-0.55e-1*x3^2)/(0.1e-3+10.*x3^2+x3)^2;
    end

    function res = h2_DD( self, x3 )
      res = (-0.11e-5-0.33e-4*x3+1.1*x3^3)/(0.1e-3+10.*x3^2+x3)^3;
    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( ~, ~, tL, tR, XL, XR, UC, ~ )
      L = 0.00001*UC(1)^2;
    end
    %
    function gradL = lagrange_gradient( self, ~, tL, tR, XL, XR, UC, ~ )
      gradL = [0,0,0,0,0,0,0,0,0.00002*UC(1)];
    end
    %
    function hessL = lagrange_hessian( self, ~, tL, tR, XL, XR, UC, ~ )
      dim = 2*self.nx+self.nu;
      hessL = sparse( dim, dim );
      hessL(9,9) = 0.00002;
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      dim = 2*self.nx+self.nu;
      hessL = sparse( dim, dim );
      hessL(9,9) = 1;
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( ~, tL, tR, XL, XR, ~ )
      M = -XR(2)*XR(4);
    end
    %
    function gradM = mayer_gradient( ~, tL, tR, XL, XR, ~ )
      gradM = -[0,0,0,0,0,XR(4),0,XR(2)];
    end
    %
    function hessM = mayer_hessian( ~, tL, tR, XL, XR, ~ )
      hessM      = sparse(8,8);
      hessM(6,8) = -1;
      hessM(8,6) = -1;
    end
    %
    function hessM = mayer_hessian_pattern( ~ )
      hessM      = sparse(8,8);
      hessM(6,8) = 1;
      hessM(8,6) = 1;
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function RES = ODE( self, ~, t, X, U, P )
      RES = [ self.h1(X(1),X(3))*X(1) - U*(X(1)/(500*X(4))); ...
              self.h2(X(3))*X(1) - 0.01*X(2)- U*(X(2)/(500*X(4))); ...
              - self.h1(X(1),X(3))*X(1)/0.47 - self.h2(X(3))*X(1)/1.2 ...
              - X(1)*(0.029*X(3)/(0.0001+X(3)))+ U/X(4)*(1-X(3)/500.); ...
              U/500 ];
    end

    function JAC = JODE( self, ~, t, X, U, P )
      JAC = sparse(4,5);
      JAC(1,1) = self.h1_1(X(1),X(3))*X(1) + self.h1(X(1),X(3)) - U/(500*X(4));
      JAC(1,3) = self.h1_3(X(1), X(3))*X(1);
      JAC(1,4) = U*X(1)/(500*X(4)^2);
      JAC(1,5) = -X(1)/(500*X(4));

      JAC(2,1) = self.h2(X(3));
      JAC(2,2) = -0.01-U/(500*X(4));
      JAC(2,3) = self.h2_D(X(3))*X(1);
      JAC(2,4) = U*X(2)/(500*X(4)^2);
      JAC(2,5) = -X(2)/(500*X(4));

      JAC(3,1) = -2.127659574*self.h1_1(X(1), X(3))*X(1)-2.127659574*self.h1(X(1), X(3))-.8333333333*self.h2(X(3))-0.29e-1*X(3)/(0.1e-3+X(3));
      JAC(3,3) = -2.127659574*self.h1_3(X(1), X(3))*X(1)-.8333333333*self.h2_D(X(3))*X(1)-0.29e-1*X(1)/(0.1e-3+X(3))+0.29e-1*X(1)*X(3)/(0.1e-3+X(3))^2-0.2000000000e-2*U/X(4);
      JAC(3,4) = U*(-1.+0.2e-2*X(3))/X(4)^2;
      JAC(3,5) = (1.-0.2e-2*X(3))/X(4);

      JAC(4,5) = 1/500;
    end

    function JAC = JODE_pattern( self )
      JAC = sparse(4,5);
      JAC(1,1) = 1;
      JAC(1,3) = 1;
      JAC(1,4) = 1;
      JAC(1,5) = 1;

      JAC(2,1) = 1;
      JAC(2,2) = 1;
      JAC(2,3) = 1;
      JAC(2,4) = 1;
      JAC(2,5) = 1;

      JAC(3,1) = 1;
      JAC(3,3) = 1;
      JAC(3,4) = 1;
      JAC(3,5) = 1;

      JAC(4,5) = 1;
    end

    function C = ds( self, nseg, tL, tR, XL, XR, UC, P )
      C = self.midpoint_ds( nseg, tL, tR, XL, XR, UC, P, @self.ODE );
    end
    %
    function JAC = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, P )
      JAC = self.midpoint_ds_jacobian( nseg, tL, tR, XL, XR, UC, P, @self.JODE );
    end
    %
    function JAC = ds_jacobian_pattern( self )
      JAC = self.midpoint_ds_jacobian_pattern( @self.JODE_pattern );
    end
    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, P, L )
      H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, P, L );
    end
    %
    function H = ds_hessian_pattern( self )
      dim = 2*self.nx + self.nu + self.np;
      H = sparse( ones( dim, dim ) );
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
             XL(4) - self.x4_i ];
    end
    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, ~ )
      Jac = sparse(4,8);
      Jac(1:4,1:4) = speye(4,4);
    end
    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse(4,8);
      Jac(1:4,1:4) = speye(4,4);
    end
    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, L )
      Hess = sparse(8,8);
    end
    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = sparse(8,8);
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
