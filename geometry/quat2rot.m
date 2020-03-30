% quat2rot - convert a quaternion(4x1) to a rotation matrix (3x3) 
%
%    R = quat2rot(q)
% 
%    q - 4x1 quaternion
%    R - 3x3 rotation matrix
%        q = [cos(theta/2), sin(theta/2) * v]
%        theta: rotation angle
%        v    : unit rotation axis, |v| = 1
%
% See also: rot2quat
%           Horn's paper

% Manolis Lourakis 2011

function R = quat2rot(q)

  p=q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2;
  if(p<0.99999 || p>1.00001), 
    disp('Warning: quat2rot: quaternion is not of unit norm');
    disp(1-p);
  end
  q=q./sqrt(p); % normalize
   
  % q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2 = 1
  R=zeros(3);
  R(1,1)=1-2*(q(3)*q(3)+q(4)*q(4));
  R(1,2)=2*(q(2)*q(3)-q(1)*q(4));
  R(1,3)=2*(q(2)*q(4)+q(1)*q(3));

  R(2,1)=2*(q(2)*q(3)+q(1)*q(4));
  R(2,2)=1-2*(q(2)*q(2)+q(4)*q(4));
  R(2,3)=2*(q(3)*q(4)-q(1)*q(2));

  R(3,1)=2*(q(2)*q(4)-q(1)*q(3));
  R(3,2)=2*(q(3)*q(4)+q(1)*q(2));
  R(3,3)=1-2*(q(2)*q(2)+q(3)*q(3));
   
return
