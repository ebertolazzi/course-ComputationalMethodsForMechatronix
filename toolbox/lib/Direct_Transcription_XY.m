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
classdef Direct_Transcription_XY < handle

  properties (SetAccess = protected, Hidden = true)
    problem;   % class storing the problem to be solved
    n;         % number of discretisation points
    a;         % left border
    b;         % right border
    ipopt_ck;
    t_sol;     % sampled t
    x_sol;     % sampled approximate solution
    y_sol;     % sampled approximate solution
  end

  methods
    function self = Direct_Transcription_XY( prob, n )
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
    function [x,y,a,b,n] = get( self, xy )
      a  = self.a;
      b  = self.b;
      n  = self.n;
      % ----------
      x  = xy(1:2:end);
      y  = xy(2:2:end);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = eval_F( self, xy )
      %
      %  evaluate the target: the discrete integral of L( x, y, x', y', t)
      %  using midpoint rule (scaled by 1/dt)
      %
      [x,y,a,b,n] = self.get( xy );
      F  = 0;
      dt = (b-a)/(n-1);
      for k=1:n-1
        tk = a + dt * (k-1/2);
        xm = (x(k+1)+x(k))/2;
        ym = (y(k+1)+y(k))/2;
        xd = (x(k+1)-x(k))/dt;
        yd = (y(k+1)-y(k))/dt;
        F  = F + self.problem.eval_L( xm, ym, xd, yd, tk );
      end
      %
      ni = self.problem.get_number_of_ic();
      if ni > 0
        IC = self.eval_int_IC( xy );
        F  = F + dot(IC,IC);
      end
      %F = F + self.eval_Penalty( xy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function G = eval_GradF( self, xy )
      %
      %  evaluate the gradient of the target
      %
      [x,y,a,b,n] = self.get( xy );
      dt = (b-a)/(n-1);
      G  = zeros( 1, 2*n );
      kk = 1;
      for k=1:n-1
        tk = a + dt * (k-1/2);
        xm = (x(k+1)+x(k))/2;
        ym = (y(k+1)+y(k))/2;
        xd = (x(k+1)-x(k))/dt;
        yd = (y(k+1)-y(k))/dt;
        L_1 = self.problem.eval_L_D_1( xm, ym, xd, yd, tk );
        L_2 = self.problem.eval_L_D_2( xm, ym, xd, yd, tk );
        L_3 = self.problem.eval_L_D_3( xm, ym, xd, yd, tk );
        L_4 = self.problem.eval_L_D_4( xm, ym, xd, yd, tk );
        G(kk)   = G(kk)   + L_1/2 - L_3/dt;
        G(kk+2) = G(kk+2) + L_1/2 + L_3/dt;
        G(kk+1) = G(kk+1) + L_2/2 - L_4/dt;
        G(kk+3) = G(kk+3) + L_2/2 + L_4/dt;
        kk      = kk+2;
      end
      ni = self.problem.get_number_of_ic();
      if ni > 0
        IC      = self.eval_int_IC( xy );
        Grad_IC = self.eval_Grad_int_IC( xy );
        for k=1:ni
          G = G + IC(k)*Grad_IC(k,:);
        end
      end
      %G = G + self.eval_GradPenalty( xy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function P = eval_Penalty( self, xy )
      %
      %  evaluate the target: the discrete integral of L( x, y, x', y', t)
      %  using midpoint rule (scaled by 1/dt)
      %
      [x,y,a,b,n] = self.get( xy );
      P  = 0;
      dt = (b-a)/(n-1);
      for k=2:n-1
        dxL = x(k)-x(k-1); dxR = x(k+1)-x(k);
        dyL = y(k)-y(k-1); dyR = y(k+1)-y(k);
        dL  = dxL^2+dyL^2; dR  = dxR^2+dyR^2;
        P   = P + (dL-dR)^2;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function GradP = eval_GradPenalty( self, xy )
      %
      %  evaluate the target: the discrete integral of L( x, y, x', y', t)
      %  using midpoint rule (scaled by 1/dt)
      %
      [x,y,a,b,n] = self.get( xy );
      GradP = zeros(1,2*n);
      dt = (b-a)/(n-1);
      for k=2:n-1
        dxL = x(k)-x(k-1); dxR = x(k+1)-x(k);
        dyL = y(k)-y(k-1); dyR = y(k+1)-y(k);
        dL  = dxL^2+dyL^2; dR  = dxR^2+dyR^2;
        tmp = 4*(dL-dR);
        gLx = -tmp*dxL;
        gCx = tmp*(dxL+dxR);
        gRx = -tmp*dxR;
        gLy = -tmp*dyL;
        gCy = tmp*(dyL+dyR);
        gRy = -tmp*dyR;
        GradP(2*(k-1)-1) = GradP(2*(k-1)-1)+gLx;
        GradP(2*(k)-1)   = GradP(2*(k)-1)  +gCx;
        GradP(2*(k+1)-1) = GradP(2*(k+1)-1)+gRx;
        GradP(2*(k-1))   = GradP(2*(k-1))  +gLy;
        GradP(2*(k))     = GradP(2*(k))    +gCy;
        GradP(2*(k+1))   = GradP(2*(k+1))  +gRy;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function IC = eval_int_IC( self, xy )
      %
      %  evaluate the integral constraints: the discrete integral of IC(x,x',t)
      %  using midpoint rule (scaled by 1/dt)
      %
      ni = self.problem.get_number_of_ic();
      if ni > 0
        [x,y,a,b,n] = self.get( xy );
        dt = (b-a)/(n-1);
        % ----------
        IC = zeros(ni,1);
        for k=1:n-1
          tk = a + dt * (k-1/2);
          xm = (x(k+1)+x(k))/2;
          ym = (y(k+1)+y(k))/2;
          xd = (x(k+1)-x(k))/dt;
          yd = (y(k+1)-y(k))/dt;
          IC = IC + self.problem.eval_IC( xm, ym, xd, yd, tk );
        end
      else
        IC = [];
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function GradIC = eval_Grad_int_IC( self, xy )
      %
      %  evaluate the gradient(s) of the integral constraints
      %
      ni = self.problem.get_number_of_ic();
      if ni > 0
        [x,y,a,b,n] = self.get( xy );
        dt = (b-a)/(n-1);
        GradIC = zeros( ni, 2*n );
        kk = 1;
        for k=1:n-1
          tk = a + dt * (k-1/2);
          xm = (x(k+1)+x(k))/2;
          ym = (y(k+1)+y(k))/2;
          xd = (x(k+1)-x(k))/dt;
          yd = (y(k+1)-y(k))/dt;
          L_1 = self.problem.eval_IC_D_1( xm, ym, xd, yd, tk );
          L_2 = self.problem.eval_IC_D_2( xm, ym, xd, yd, tk );
          L_3 = self.problem.eval_IC_D_3( xm, ym, xd, yd, tk );
          L_4 = self.problem.eval_IC_D_4( xm, ym, xd, yd, tk );
          GradIC(:,kk)   = GradIC(:,kk)   + L_1/2 - L_3/dt;
          GradIC(:,kk+2) = GradIC(:,kk+2) + L_1/2 + L_3/dt;
          GradIC(:,kk+1) = GradIC(:,kk+1) + L_2/2 - L_4/dt;
          GradIC(:,kk+3) = GradIC(:,kk+3) + L_2/2 + L_4/dt;
          kk             = kk+2;
        end
      else
        GradIC = zeros(0,2*n);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function C = eval_C( self, xy )
      %
      %  evaluate the equality constraints as union
      %  of the boundary conditions and integral constraints
      %
      [x,y,a,b,n] = self.get( xy );
      % ----------
      dt  = (b-a)/(n-1);
      dxa = (x(2)-x(1))/dt;
      dya = (y(2)-y(1))/dt;
      dxb = (x(n)-x(n-1))/dt;
      dyb = (y(n)-y(n-1))/dt;
      C   = [ ...
        self.problem.eval_BC( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );...
        self.eval_int_IC( xy ) ...
      ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function JC = eval_JC( self, xy )
      %
      %  evaluate the equality constraints as union
      %  of the boundary conditions and integral constraints
      %
      [x,y,a,b,n] = self.get( xy );
      dt = (b-a)/(n-1);
      nb = self.problem.get_number_of_bc();
      ni = self.problem.get_number_of_ic();
      % ----------
      dxa = (x(2)-x(1))/dt;
      dya = (y(2)-y(1))/dt;
      dxb = (x(n)-x(n-1))/dt;
      dyb = (y(n)-y(n-1))/dt;

      D1   = self.problem.eval_BC_D_1( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D2   = self.problem.eval_BC_D_2( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D3   = self.problem.eval_BC_D_3( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D4   = self.problem.eval_BC_D_4( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D5   = self.problem.eval_BC_D_5( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D6   = self.problem.eval_BC_D_6( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D7   = self.problem.eval_BC_D_7( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );
      D8   = self.problem.eval_BC_D_8( x(1), y(1), dxa, dya, x(n), y(n), dxb, dyb );

      JC           = sparse( nb+ni, 2*n );
      JC(1:nb,1)   = D1 - D3/dt;
      JC(1:nb,3)   = D3/dt;
      JC(1:nb,2)   = D2 - D4/dt;
      JC(1:nb,4)   = D4/dt;

      JC(1:nb,2*n-3) = - D7/dt;
      JC(1:nb,2*n-1) = D5 + D7/dt;
      JC(1:nb,2*n-2) = - D8/dt;
      JC(1:nb,2*n)   = D6 + D8/dt;

      if ni > 0
        JC(nb+1:nb+ni,:) = self.eval_Grad_int_IC( xy );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function JC = eval_JC_pattern( self )
      %
      %  evaluate jacobian of the equality constraints
      %
      %
      %  evaluate the equality constraints as union
      %  of the boundary conditions and integral constraints
      %
      a  = self.a;
      b  = self.b;
      n  = self.n;
      nb = self.problem.get_number_of_bc();
      ni = self.problem.get_number_of_ic();

      JC           = sparse( nb+ni, 2*n );
      JC(1:nb,1)   = ones(nb,1);
      JC(1:nb,3)   = ones(nb,1);
      JC(1:nb,2)   = ones(nb,1);
      JC(1:nb,4)   = ones(nb,1);

      JC(1:nb,2*n-3) = ones(nb,1);
      JC(1:nb,2*n-1) = ones(nb,1);
      JC(1:nb,2*n-2) = ones(nb,1);
      JC(1:nb,2*n)   = ones(nb,1);

      if ni > 0
        JC(nb+1:nb+ni,:) = ones(ni,2*n);
      end

    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = solve( self, max_iter, xmin, ymin, xmax, ymax )
      a  = self.a;
      b  = self.b;
      n  = self.n;
      dt = (b-a)/(n-1);
      % ----------
      x_guess = zeros( n, 1 );
      y_guess = zeros( n, 1 );
      t_sol   = a+(0:n-1)*dt;
      for k=1:n
        [x,y] = self.problem.eval_Guess( t_sol(k) );
        x_guess(k) = x;
        y_guess(k) = y;
      end
      ok = self.solve2( n, x_guess, y_guess, max_iter, xmin, ymin, xmax, ymax );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = solve2( self, n, x_guess, y_guess, max_iter, xmin, ymin, xmax, ymax )
      self.n = n;
      a     = self.a;
      b     = self.b;
      dt    = (b-a)/(n-1);
      t_sol = a+(0:n-1)*dt;
      xy_guess = zeros( 2*n, 1 );
      for k=1:n
        xy_guess(2*k-1) = x_guess(k);
        xy_guess(2*k)   = y_guess(k);
      end
      %xy_guess = reshape( [x_guess(:).';y_guess(:).'], 2*n, 1 );

      options = {};

      nb = self.problem.get_number_of_bc();
      ni = self.problem.get_number_of_ic();

      options.lb = zeros(1,2*n);
      options.lb(1:2:2*n) = xmin*ones(n,1);
      options.lb(2:2:2*n) = ymin*ones(n,1);

      options.ub = zeros(1,2*n);
      options.ub(1:2:2*n) = xmax*ones(n,1);
      options.ub(2:2:2*n) = ymax*ones(n,1);

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
      funcs.objective         = @( xy ) self.eval_F( xy );
      funcs.gradient          = @( xy ) self.eval_GradF( xy );
      funcs.constraints       = @( xy ) self.eval_C( xy );
      funcs.jacobian          = @( xy ) self.eval_JC( xy );
      funcs.jacobianstructure = @(    ) self.eval_JC_pattern();

      %options.ipopt.jacobian_approximation = 'finite-difference-values';
      options.ipopt.hessian_approximation  = 'limited-memory';
      %options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      %options.ipopt.limited_memory_update_type = 'sr1';
      options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      options.ipopt.limited_memory_max_history = 40;

      [xy_sol, ok] = ipopt( xy_guess, funcs, options );
      self.x_sol = xy_sol(1:2:end);
      self.y_sol = xy_sol(2:2:end);
      self.t_sol = t_sol;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [t,x,y] = get_solution( self )
      t = self.t_sol;
      x = self.x_sol;
      y = self.y_sol;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self )
      plot( self.x_sol, self.y_sol, '-o', 'Linewidth', 2 );
    end
  end
end
