% Model selection with Torr's GRIC
% Computes GRIC the scores for the supplied fundamental matrix and homography;
% smallest is better. Returns 1 or 2 in model depending on whether the points
% are best decribed by a fundamental matrix or a homography.
% See "An assessment of information criteria for motion model selection", CVPR 1997
%
% m1, m2 are matched point pairs, including outliers (2xN)
% F, H are fundamental matrix and homography estimates and
% sigma is the assumed variance of the error, e.g. 0.8

% Manolis Lourakis 2018
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece

function [model gricF gricH] = gric_modsel(m1, m2, F, H, sigma)

  Nt = size(m1, 2);
  R = 4; % data dimension (point pairs)

  % make point sets homogeneous
  if(size(m1, 1)==2), m1=[m1; ones(1, size(m1, 2))]; end
  if(size(m2, 1)==2), m2=[m2; ones(1, size(m2, 2))]; end

  % fundamental
  K = 7; % number of parameters
  D = 3; % manifold dimension

  resF = F_sampson_distance_sqr(F, m1, m2);
  resF = min(resF/sigma^2, ones(size(resF))*2*(R-D));
  gF = sum(resF) + Nt*D*log(R) + K*log(R*Nt);


  % homography
  K = 8; % number of parameters
  D = 2; % manifold dimension

  resH = vgg_H_sampson_distance_sqr(H, m1, m2);
  resH = min(resH/sigma^2, ones(size(resH))*2*(R-D));
  gH = sum(resH) + Nt*D*log(R) + K*log(R*Nt);

  if(gF<gH), model=1; else model=2; end

  if(nargout>=2), gricF=gF; end
  if(nargout>=3), gricH=gH; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  d = F_sampson_distance_sqr(F, m1, m2);

% sampson(F, m1, m2)  valuate the first order approximation of the geometric error
% of the fit of F  with respect to a set of matched points m1, m2.
%
% Returns an approximation of the **squared** geometric distance from the
% 4d joint-space point [m1;m2] to the F manifold. 
%
% Reference: Luong, Faugeras. IJCV 1996 (called 'Gradient criterion")  

% Author: A. Fusiello. Extracted and adapted from a piece of code by Peter Kovesi 
% The interface is the same as vgg_H_sampson_distance_sqr (see)

% Minimal changes for vectorization by M. Lourakis

%lng = size(m1,2);

%m2tFm1 = zeros(1,lng);
%for n = 1:lng
%    m2tFm1(n) = m2(:,n)'*F*m1(:,n);
%end

Fm1 = F*m1;
Ftm2 = F'*m2;

m2tFm1 = sum(m2.*Fm1);

% Evaluate squared distances
d =  m2tFm1.^2 ./ (Fm1(1,:).^2 + Fm1(2,:).^2 + Ftm2(1,:).^2 + Ftm2(2,:).^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = vgg_H_sampson_distance_sqr(H,X1,X2)

% d = vgg_H_sampson_distance_sqr(H,X1,X2)
%
% Returns an approximation of the squared geometric distance from the
% 4d joint-space point [X1;X2] to the H manifold. See Torr CVIU 97.

if (size(X1) ~= size(X2))
  error('Point sets not same size!');
end

% ensure that the third coord is 1 avoiding unnecessary divisions
if(sum(X1(3,:)~=1)>0), p1=X1./repmat(X1(3,:), 3,1); else, p1=X1; end
if(sum(X2(3,:)~=1)>0), p2=X2./repmat(X2(3,:), 3,1); else, p2=X2; end

%p1 = X1 ./ repmat(X1(3,:),3,1);
%p2 = X2 ./ repmat(X2(3,:),3,1);

%alg = vgg_H_algebraic_distance(H,p1,p2);
N = size(p1,2);
Dx = [ p1' .* repmat(p2(3,:)',1,3) , zeros(N,3) , -p1' .* repmat(p2(1,:)',1,3) ];
Dy = [ zeros(N,3) , p1' .* repmat(p2(3,:)',1,3) , -p1' .* repmat(p2(2,:)',1,3) ];
h = reshape(H',9,1);
alg = [Dx * h , Dy * h]';

%N = size(X1,2);

%h = reshape(H',9,1);

G1 = [ H(1,1) - p2(1,:) * H(3,1) ; ...
  H(1,2) - p2(1,:) * H(3,2) ; ...      
  -p1(1,:) * H(3,1) - p1(2,:) * H(3,2) - H(3,3) ; ...
  zeros(1,N) ];

G2 = [ H(2,1) - p2(2,:) * H(3,1) ; ...
  H(2,2) - p2(2,:) * H(3,2) ; ...
  zeros(1,N) ; ...
  -p1(1,:) * H(3,1) - p1(2,:) * H(3,2) - H(3,3) ];

magG1 = sqrt(sum(G1 .* G1));
magG2 = sqrt(sum(G2 .* G2));
magG1G2 = sum(G1 .*  G2);

alpha = acos( magG1G2 ./ (magG1 .* magG2) );

D1 = alg(1,:) ./ magG1;
D2 = alg(2,:) ./ magG2;

d = (D1.*D1 + D2.*D2 - 2 * D1 .* D2 .* cos(alpha)) ./ sin(alpha);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
