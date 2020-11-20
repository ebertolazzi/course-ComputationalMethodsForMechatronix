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
classdef Problem1 < Problem_Base_1D

  properties (SetAccess = protected)
    a;
    b;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %     / a
    %    |     x'(t)^2 + 10 * t * x(t) dt
    %   / b
    %
    %   x(a) = 1, x(b) = 2
    %
    %
    % evaluate the function L(x,x',t)
    function x = eval_Guess( self, t )
      x = 1+t;
    end

    % evaluate the function L(x,x',t)
    function L = eval_L( self, x, z, t )
      L = z*z + 10*t*x;
    end

    % evaluate the function DL(x,x',t) / Dx
    function L_D_1 = eval_L_D_1( self, x, z, t )
      L_D_1 = 10*t;
    end

    % evaluate the function DL(x,x',t) / Dx'
    function L_D_2 = eval_L_D_2( self, x, z, t )
      L_D_2 = 2*z;
    end

    % evaluate the function BC(x(a),x'(a),x(b),x'(b))
    function BC = eval_BC( self, xa, xa_dot, xb, xb_dot )
       BC = [ xa - 1; xb - 2];
    end

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx(a)
    function BC_D_1 = eval_BC_D_1( self, xa, xa_dot, xb, xb_dot )
      BC_D_1 = [ 1; 0];
    end

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx'(a)
    function BC_D_2 = eval_BC_D_2( self, xa, xa_dot, xb, xb_dot )
      BC_D_2 = [ 0; 0];
    end

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx(b)
    function BC_D_3 = eval_BC_D_3( self, xa, xa_dot, xb, xb_dot )
      BC_D_3 = [ 0; 1];
    end

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx'(b)
    function BC_D_4 = eval_BC_D_4( self, xa, xa_dot, xb, xb_dot )
      BC_D_4 = [ 0; 0];
    end

    % evaluate the function IC(x,x',t)
    function IC = eval_IC( self, x, x_dot, t )
      IC = [];
    end

    % evaluate the function DIC(x,x',t)/Dx
    function IC_D_1 = eval_IC_D_1( self, x, x_dot, t )
      IC_D_1 = [];
    end

    % evaluate the function DIC(x,x',t)/Dx'
    function IC_D_2 = eval_IC_D_2( self, x, x_dot, t )
      IC_D_2 = [];
    end

  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Problem1()
      self@Problem_Base_1D( 'Example taked for Calogero exercize', 2, 0 );
      self.a = 0;
      self.b = 1;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t )
      res = (5/6)*t.^3 + t/6 + 1;
    end
  end
end
