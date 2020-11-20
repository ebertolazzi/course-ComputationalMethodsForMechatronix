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
classdef Problem5 < Problem_Base_1D

  properties (SetAccess = protected)
    a;
    b;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %     / a
    %    |     L( x, x', t ) dt
    %   / b
    %
    % evaluate the function L(x,x',t)
    function y = eval_Guess( self, x )
      y = 1;
    end

    % evaluate the function L(y,y',x)
    function L = eval_L( self, y, z, x )
      L = z*z;
    end

    % evaluate the function DL(y,y',x) / Dy
    function L_D_1 = eval_L_D_1( self, y, z, t )
      L_D_1 = 0;
    end

    % evaluate the function DL(y,y',t) / Dy'
    function L_D_2 = eval_L_D_2( self, y, z, x )
      L_D_2 = 2*z;
    end

    % evaluate the function BC(y(a),y'(a),y(b),y'(b))
    function BC = eval_BC( self, ya, za, yb, zb )
       BC = [ ya-1; yb-1];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy(a)
    function BC_D_1 = eval_BC_D_1( self, ya, za, yb, zb )
      BC_D_1 = [ 1; 0];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy'(a)
    function BC_D_2 = eval_BC_D_2( self, ya, za, yb, zb)
      BC_D_2 = [ 0; 0];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy(b)
    function BC_D_3 = eval_BC_D_3( self, ya, za, yb, zb )
      BC_D_3 = [ 0; 1];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy'(b)
    function BC_D_4 = eval_BC_D_4( self, ya, za, yb, zb )
      BC_D_4 = [ 0; 0];
    end

    % evaluate the function IC(x,x',t)
    function IC = eval_IC( self, y, z, t )
      IC = y-2;
    end

    % evaluate the function DIC(x,x',t)/Dx
    function IC_D_1 = eval_IC_D_1( self, y, z, t )
      IC_D_1 = 1;
    end

    % evaluate the function DIC(x,x',t)/Dx'
    function IC_D_2 = eval_IC_D_2( self, y, z, t )
      IC_D_2 = 0;
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Problem5()
      self@Problem_Base_1D( 'isoperimetric example', 2, 1 );
      self.a   = 0;
      self.b   = 1;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, x )
      res = 5/2-6*(x-1/2).^2;
    end
  end
end
