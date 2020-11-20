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
classdef Problem2 < Problem_Base_1D

  properties (SetAccess = protected)
    a;
    b;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %     / a
    %    |     (1+t)*x'(t)^2) dt
    %   / b
    %
    %  x(a) = 0
    %  x(b) = 1
    %
    % evaluate the function L(x,x',t)
    function x = eval_Guess( self, t )
      x = 0;
    end

    % evaluate the function L(x,x',t)
    function L = eval_L( self, x, z, t )
      L = (1+t)*z*z;
    end

    % evaluate the function DL(x,x',t) / Dx
    function L_D_1 = eval_L_D_1( self, x, z, t )
      L_D_1 = 0;
    end

    % evaluate the function DL(x,x',t) / Dx'
    function L_D_2 = eval_L_D_2( self, x, z, t )
      L_D_2 = 2*(1+t)*z;
    end

    % evaluate the function BC(x(a),x'(a),x(b),x'(b))
    function BC = eval_BC( self, xa, xa_dot, xb, xb_dot )
       BC = [ xa; xb-1];
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
    function self = Problem2( a, b )
      self@Problem_Base_1D( '1D-test2', 2, 0 );
      self.a = a;
      self.b = b;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t )
      res = log(1+t)/log(2);
    end
  end
end
