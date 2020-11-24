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
classdef Problem7 < Problem_Base_XY

  properties (SetAccess = protected)
    a;
    b;
    len;
  end

  methods
    % evaluate the function L(x,x',t)
    function [x,y] = eval_Guess( self, t )
      x = t;
      y = 0;
      % set y in such a way maintain length
      if (t>0) && (t<1)
        y = (1-self.len)/2; 
      end
    end

    % evaluate the function L( x, y, x', y', t )
    function L = eval_L( self, x, y, xp, yp, t )
      L = x*yp-y*xp;
    end
    function res = eval_L_D_1( self, x, y, xp, yp, t )
      res = yp;
    end
    function res = eval_L_D_2( self, x, y, xp, yp, t )
      res = -xp;
    end
    function res = eval_L_D_3( self, x, y, xp, yp, t )
      res = -y;
    end
    function res = eval_L_D_4( self, x, y, xp, yp, t )
      res = x;
    end
    % evaluate the function BC(y(a),y'(a),y(b),y'(b))
    function res = eval_BC( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ xa; ya; xb-1; yb ];
    end
    function res = eval_BC_D_1( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 1; 0; 0; 0];
    end
    function res = eval_BC_D_2( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 1; 0; 0];
    end
    function res = eval_BC_D_3( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 0; 0; 0];
    end
    function res = eval_BC_D_4( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 0; 0; 0];
    end
    function res = eval_BC_D_5( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 0; 1; 0];
    end
    function res = eval_BC_D_6( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 0; 0; 1];
    end
    function res = eval_BC_D_7( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 0; 0; 0];
    end
    function res = eval_BC_D_8( self, xa, ya, xpa, ypa, xb, yb, xpb, ypb )
      res = [ 0; 0; 0; 0];
    end
    % evaluate the function IC( x, y, x', y', t )
    function IC = eval_IC( self, x, y, xp, yp, t )
      len = self.len;
      IC  = sqrt(xp*xp+yp*yp)-len;
    end
    function res = eval_IC_D_1( self, x, y, xp, yp, t )
      res = 0;
    end
    function res = eval_IC_D_2( self, x, y, xp, yp, t )
      res = 0;
    end
    function res = eval_IC_D_3( self, x, y, xp, yp, t )
      res = xp/sqrt(xp*xp+yp*yp);
    end
    function res = eval_IC_D_4( self, x, y, xp, yp, t )
      res = yp/sqrt(xp*xp+yp*yp);
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Problem7( len )
      self@Problem_Base_XY( 'DIDO problem', 4, 1 );
      self.len = len;
      self.a   = 0;
      self.b   = 1;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [x,y] = exact( self, t )
      x = t;
      y = zeros(size(t));
    end
  end
end
