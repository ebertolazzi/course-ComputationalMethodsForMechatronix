%
% Matlab code for the Course:
%
%     Computational Methods for Mechatronix
%
% by
% Enrico Bertolazzi
% Dipartimento di Ingegneria Industriale
% Universita` degli Studi di Trento
% email: enrico.bertolazzi@unitn.it
%
classdef Direct_Transcription_1D < handle

  properties (SetAccess = protected, Hidden = true)
    problem;   % class storing the problem to be solved
    n;         % number of discretisation points
    a;         % left border
    b;         % right border
    ipopt_ck;
    t_sol;     % sampled t
    x_sol;     % sampled approximate solution
  end

  methods
    function self = Direct_Transcription_1D( prob, n )
      self.n        = n;
      self.problem  = prob;
      self.a        = prob.a;
      self.b        = prob.b;
      self.ipopt_ck = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = eval_F( self, x )
      %
      %  evaluate the target: the discrete integral of L(x,x',t)
      %  using midpoint rule (scaled by 1/dt)
      %
      a  = self.a;
      b  = self.b;
      n  = self.n;
      % ----------
      F  = 0;
      dt = (b-a)/(n-1);
      for k=1:n-1
        tk = a + dt * (k-1/2);
        xm = (x(k+1)+x(k))/2;
        xd = (x(k+1)-x(k))/dt;
        F  = F + self.problem.eval_L( xm, xd, tk );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function G = eval_GradF( self, x )
      %
      %  evaluate the gradient of the target
      %
      a  = self.a;
      b  = self.b;
      n  = self.n;
      dt = (b-a)/(n-1);
      % ----------
      G  = zeros( n, 1 );
      for k=1:n-1
        tk  = a + dt * (k-1/2);
        xm  = (x(k+1)+x(k))/2;
        xd  = (x(k+1)-x(k))/dt;
        L_1 = self.problem.eval_L_D_1( xm, xd, tk );
        L_2 = self.problem.eval_L_D_2( xm, xd, tk );
        G(k)   = G(k)   + L_1/2 - L_2/dt;
        G(k+1) = G(k+1) + L_1/2 + L_2/dt;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function IC = eval_int_IC( self, x )
      %
      %  evaluate the integral constraints: the discrete integral of IC(x,x',t)
      %  using midpoint rule (scaled by 1/dt)
      %
      ni = self.problem.get_number_of_ic();
      if ni > 0
        a  = self.a;
        b  = self.b;
        n  = self.n;
        dt = (b-a)/(n-1);
        % ----------
        IC = zeros(ni,1);
        for k=1:n-1
          tk = a + dt * (k-1/2);
          xm = (x(k+1)+x(k))/2;
          xd = (x(k+1)-x(k))/dt;
          IC = IC + self.problem.eval_IC( xm, xd, tk );
        end
      else
        IC = [];
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function GradIC = eval_Grad_int_IC( self, x )
      %
      %  evaluate the gradient(s) of the integral constraints
      %
      ni = self.problem.get_number_of_ic();
      n  = self.n;
      if ni > 0
        a  = self.a;
        b  = self.b;
        dt = (b-a)/(n-1);
        % ----------
        GradIC = zeros( ni, n );
        for k=1:n-1
          tk  = a + dt * (k-1/2);
          xm  = (x(k+1)+x(k))/2;
          xd  = (x(k+1)-x(k))/dt;
          L_1 = self.problem.eval_IC_D_1( xm, xd, tk );
          L_2 = self.problem.eval_IC_D_2( xm, xd, tk );
          GradIC(:,k)   = GradIC(:,k)   + L_1/2 - L_2/dt;
          GradIC(:,k+1) = GradIC(:,k+1) + L_1/2 + L_2/dt;
        end
      else
        GradIC = zeros(0,n);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function C = eval_C( self, x )
      %
      %  evaluate the equality constraints as union
      %  of the boundary conditions and integral constraints
      %
      a = self.a;
      b = self.b;
      n = self.n;
      % ----------
      dt  = (b-a)/(n-1);
      dxa = (x(2)-x(1))/dt;
      dxb = (x(n)-x(n-1))/dt;
      C   = [ ...
        self.problem.eval_BC( x(1), dxa, x(n), dxb );...
        self.eval_int_IC( x ) ...
      ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function JC = eval_JC( self, x )
      %
      %  evaluate the equality constraints as union
      %  of the boundary conditions and integral constraints
      %
      a  = self.a;
      b  = self.b;
      n  = self.n;
      dt = (b-a)/(n-1);
      nb = self.problem.get_number_of_bc();
      ni = self.problem.get_number_of_ic();
      % ----------
      dxa  = (x(2)-x(1))/dt;
      dxb  = (x(n)-x(n-1))/dt;

      D1   = self.problem.eval_BC_D_1( x(1), dxa, x(n), dxb );
      D2   = self.problem.eval_BC_D_2( x(1), dxa, x(n), dxb );
      D3   = self.problem.eval_BC_D_3( x(1), dxa, x(n), dxb );
      D4   = self.problem.eval_BC_D_4( x(1), dxa, x(n), dxb );

      JC           = sparse( nb+ni, n );
      JC(1:nb,1)   = D1 - D2/dt;
      JC(1:nb,2)   = D2/dt;
      JC(1:nb,n-1) = - D4/dt;
      JC(1:nb,n)   = D3 + D4/dt;
      if ni > 0
        JC(nb+1:nb+ni,:) = self.eval_Grad_int_IC( x );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function JC = eval_JC_pattern( self )
      %
      %  evaluate jacobian of the equality constraints
      %
      n  = self.n;
      nb = self.problem.get_number_of_bc();
      ni = self.problem.get_number_of_ic();
      % ----------
      JC = sparse( nb+ni, n );
      JC(1:nb,1)   = ones(nb,1);
      JC(1:nb,2)   = ones(nb,1);
      JC(1:nb,n-1) = ones(nb,1);
      JC(1:nb,n)   = ones(nb,1);
      if ni > 0
        JC(nb+1:nb+ni,:) = ones(ni,n);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = solve( self, max_iter, xmin, xmax )
      a  = self.a;
      b  = self.b;
      n  = self.n;
      dt = (b-a)/(n-1);
      % ----------
      x_guess = zeros( n, 1 );
      t_sol   = a+(0:n-1)*dt;
      for k=1:n
        x_guess(k) = self.problem.eval_Guess( t_sol(k) );
      end

      options = {};

      n  = self.n;
      nb = self.problem.get_number_of_bc();
      ni = self.problem.get_number_of_ic();

      options.ub = xmax*ones(n,1);
      options.lb = xmin*ones(n,1);

      % The constraint functions are bounded to zero
      options.cl = zeros(nb+ni,1); %  constraints
      options.cu = zeros(nb+ni,1);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant      = 'no';
      options.ipopt.hessian_constant    = 'no';
      options.ipopt.mu_strategy         = 'adaptive';
      options.ipopt.max_iter            = max_iter;
      options.ipopt.tol                 = 1e-8;
      options.ipopt.derivative_test_tol = 1e-4;
      if self.ipopt_ck
        options.ipopt.derivative_test = 'first-order';
      else
        options.ipopt.derivative_test = 'none';
      end
      %options.ipopt.derivative_test_perturbation = 1e-8;

      % The callback functions.
      funcs.objective         = @( x ) self.eval_F( x );
      funcs.gradient          = @( x ) self.eval_GradF( x );
      funcs.constraints       = @( x ) self.eval_C( x );
      funcs.jacobian          = @( x ) self.eval_JC( x );
      funcs.jacobianstructure = @(   ) self.eval_JC_pattern();

      %options.ipopt.jacobian_approximation = 'finite-difference-values';
      options.ipopt.hessian_approximation  = 'limited-memory';
      %options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      %options.ipopt.limited_memory_update_type = 'sr1';
      options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      options.ipopt.limited_memory_max_history = 40;

      [x_sol, info] = ipopt( x_guess, funcs, options );
      info;
      self.x_sol = x_sol;
      self.t_sol = t_sol;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self )
      plot( self.t_sol, self.x_sol, '-o', 'Linewidth', 2 );
    end
  end
end
