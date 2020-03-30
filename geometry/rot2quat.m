% rot2quat - converts a rotation matrix (3x3) to a unit quaternion (4x1)
%
%    q = rot2quat(R)
% 
%    R - 3x3 rotation matrix
%    q - 4x1 unit quaternion
%        q = [cos(theta/2) sin(theta/2) * v]
%        theta: rotation angle
%        v    :  unit rotation axis, |v| = 1
%
%    
% See also: quat2rot
%           A8 in Horn's paper, http://www.j3d.org/matrix_faq/matrfaq_latest.html

% Manolis Lourakis 2011

function q = rot2quat(R)

	% find the maximum of the 4 quantities
  [m, i]=max([
	    1.0 + R(1,1) + R(2,2) + R(3,3);
	    1.0 + R(1,1) - R(2,2) - R(3,3);
	    1.0 - R(1,1) + R(2,2) - R(3,3);
	    1.0 - R(1,1) - R(2,2) + R(3,3);
      ]);

  S=sqrt(m)*2;
	switch i
	case 1
    q=[
		  S*0.25;
		  (R(3,2) - R(2,3))/S;
		  (R(1,3) - R(3,1))/S;
		  (R(2,1) - R(1,2))/S;
      ];
	case 2
    q=[
		  (R(3,2) - R(2,3))/S;
		  S*0.25;
		  (R(2,1) + R(1,2))/S;
		  (R(1,3) + R(3,1))/S;
      ];
	case 3
    q=[
		  (R(1,3) - R(3,1))/S;
		  (R(2,1) + R(1,2))/S;
		  S*0.25;
		  (R(3,2) + R(2,3))/S;
      ];
	case 4
    q=[
		  (R(2,1) - R(1,2))/S;
		  (R(1,3) + R(3,1))/S;
		  (R(3,2) + R(2,3))/S;
		  S*0.25;
      ];
	otherwise % should not happen 
    disp('Warning: rot2quat: internal error');
	end

	% enforce unit length
  q=q./sqrt(q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);

return
