% nviewTriangulatePts - single point triangulation from image matches and camera matrices using N views
%
%    [pt3D, inliers] = nviewTriangulatePts(Ps, pts, imgsizes, howto, refine)
% 
%    Ps         - N finite camera matrices, specified either as a 3x4xN matrix or as a cell array of N 3x4 matrices
%    pts        - 2xN array of corresponding image projections
%    imgsizes   - (optional) 2xN matrix whose columns are image sizes, used for preconditioning in DLT;
%                 e.g. [640 640 640; 480 480 480] for three 640x480 images.
%                 If equal to [], a simple scaling to the columns of the DLT matrix is performed.
%    howto      - optional argument that specifies how triangulation is to be performed:
%                    * if 'dlt', the DLT method is employed to combine all image projections (default);
%                    * if 'inh', the N view inhomogeneous method is employed.
%                    * if 'midp', the N view midpoint method is employed.
%                    * if 'rob', a method tolerant to mismatches is employed; rob-xx sets the repr. error threshold to xx pix.
%    refine     - optional argument that specifies how triangulation is to be refined:
%                    * if 'vgg', the point is refined by minimizing the reprojection error in N views
%                      using VGG's vgg_X_from_xP_nonlin (which relies on Gauss-Newton);
%                    * if 'lm', the point is refined by minimizing the reprojection error in N views
%                      using Levenberg-Marquardt (this requires marquardt.m);
%                    * if 'ang', the point is refined with the N view angular minimization criterion (L-M/Gauss-Newton);
%                    * if 'none', no refinement is performed.
%                 The refinement operates on a point initially triangulated as specified by 'howto'.
%
%                 The recommended combinations for 'howto' and 'refine' are 'dlt'-'vgg' and 'midp'-'vgg'
%
% Returns
%    pt3D     - 3x1 triangulated 3D point (Euclidean)
%    inliers  - (optional) array of inlying point indices. Relevant only when howto='rob', equal to 1:N otherwise

% see HZ2, p.312, 314, 318 & "Triangulation" by Hartley-Sturm
%
% Manolis Lourakis 2007-20
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece

% July  2018  - Original version. (v. 1.0)
% Sep   2018  - Added p-norm & robust method. (v. 1.1)
% Sep   2018  - Added extra checks to X_from_xP_rob [ray angle, #inliers]. (v. 1.2)
% Oct   2018  - Added X_from_xP_linLS, X84 rule to X_from_xP_rob, cheirality test. (v. 1.3)
% Feb   2020  - Added generalized weighted midpoint calculation to X_from_xP_rob. (v. 1.4)

function [pt3D, inliers] = nviewTriangulatePts(Ps, pts, imgsizes, howto, refine)

  if(iscell(Ps))
    Ps=cat(3, Ps{:});
  end
  N=size(Ps,3);
  if(N<2)
    error('nviewTriangulatePts: cannot reconstruct 3D from one image');
  end

  if(N~=size(pts, 2))
    error('nviewTriangulatePts: Numbers of camera matrices and image projections must match');
  end

  if((size(pts, 1)~=2 && size(pts, 1)~=3))
    error('nviewTriangulatePts: can only handle 2D correspondences (got %dx%d)', size(pts, 1), N);
  end

  if(nargin<=2)
    imgsizes=[];
  else
    if(~isnumeric(imgsizes))
      error('nviewTriangulatePts: "imgsizes" must be a numeric array!');
    end

    [r, c]=size(imgsizes);
    if(r~=2 || c~=N)
      if(r*c==2)
        %warning('nviewTriangulatePts: assuming identical sizes for all images');
        imgsizes=repmat(imgsizes(:), 1, N);
      elseif(r*c~=0)
        error('nviewTriangulatePts: "imgsizes" must be a 2xN array');
      end
    end
  end

  if(nargin<=3)
    howto='dlt';
  else
    if(~ischar(howto))
      error('nviewTriangulatePts: "howto" must be a string!');
    end
    howto=lower(howto);
  end

  if(size(pts, 1)==3)
    % make point sets inhomogeneous: divide by third coord if non-unit
    if(sum(pts(3,:)~=1)>0), pts=pts(1:2,:)./repmat(pts(3,:), 2,1); else, pts=pts(1:2,:); end
  end

  if(nargout>1)
    inliers=1:size(Ps,3);
  end

  % initial triangulation
  if(strncmp(howto, 'dlt', 3))
    pt3D=X_from_xP_lin(Ps, pts, imgsizes);
  elseif(strncmp(howto, 'inh', 3))
    pt3D=X_from_xP_Lp(Ps, pts); % L2
  elseif(strncmp(howto, 'midp', 3))
    pt3D=X_from_xP_mid(Ps, pts);
  elseif(strncmp(howto, 'rob', 3))
    thr=sscanf(howto, 'rob-%f');
    if isempty(thr), thr=1.0; end
    [pt3D, inl]=X_from_xP_rob(Ps, pts, imgsizes, 0.05, thr); % common perpendicular segment < 0.05*baseline, reprojection error <= thr
    if(sum(isinf(pt3D))), return; end

    % remove outliers
    Ps=Ps(:,:,inl);
    pts=pts(:,inl);
    if ~isempty(imgsizes), imgsizes=imgsizes(:,inl); end

    if(nargout>1), inliers=inl; end
  else
    error('nviewTriangulatePts: unknown triangulation method!');
  end

  if(nargin<=4 || strncmp(refine, 'none', 3))
    return;
  end

  % refinement
  if(strncmp(refine, 'vgg', 3))
    if(isempty(imgsizes)), imgsizes=repmat([2;2], 1, N); end % workaround: vgg_X_from_xP_nonlin cannot handle empty "imgsizes"
    pt3D=vgg_X_from_xP_nonlin(pts, Ps, imgsizes, [pt3D;1]);
    pt3D=pt3D(1:3)./pt3D(4);
  elseif(strncmp(refine, 'lm', 2))
    pt3D=X_from_xP_nonlin_levmar(pts, Ps, imgsizes, [pt3D;1]);
    pt3D=pt3D(1:3)./pt3D(4);
  elseif(strncmp(refine, 'ang', 3))
    pt3D=X_from_xP_ang(Ps, pts, pt3D);
  else
    error('nviewTriangulatePts: unknown triangulation refinement method!');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_from_xP_lin  Estimation of 3D point from image matches and camera matrices with
%   the DLT (Linear-Eigen) method
%
%   X =  X_from_xP_lin(P,x,imsize) computes 3D point X (column 3-vector)
%   from its projections in K images x (2-by-K matrix) and camera matrices P (K-cell
%   of 3-by-4 matrices). Image sizes imsize (2-by-K matrix) are needed for preconditioning.
%
%   See also vgg_X_from_xP_lin, X_from_xP_linLS.

% lourakis at ics.forth.gr, 2018

function X = X_from_xP_lin(P,u,imsize)

  K=size(P,3);

  if nargin>2 && ~isempty(imsize)
    for k=1:K
      H = [2/imsize(1,k) 0 -1
           0 2/imsize(2,k) -1
           0 0              1];
      P(:,:,k)=H*P(:,:,k);
      u(:,k)=H(1:2,1:2)*u(:,k) + H(1:2,3);
    end
  end

  A=zeros(2*K, 4);
  for k=1:K
    Pk=P(:,:,k);
    A(2*(k-1)+1:2*k, :)=[u(:,k)*Pk(3,:)-Pk(1:2,:)];
  end

  %{
  A=zeros(3*K, 4+K);
  for k=1:K
    %A(3*(k-1)+1:3*k, 1:4)=P(:,:,k);
    %A(3*(k-1)+1:3*k, 4+k)=-[u(:,k); 1];
    A(3*(k-1)+1:3*k, [1:4, 4+k])=[P(:,:,k), -[u(:,k); 1]];
  end
  %}

  if(nargin<=2 || isempty(imsize)) % no scaling specified, use diagonal scaling matrix
    D=diag(1./max(abs(A))); % A = A*D*inv(D)
    A=A*D;
    [dummy dummy V]=svd(A, 0);
    V(:,end)=D*V(:,end); % undo normalization
    X=V(1:3,end)./V(4,end);
  else
    [dummy dummy V]=svd(A, 0);
    X=V(1:3,end)./V(4,end);
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_from_xP_linLS  Estimation of 3D point from image matches and camera matrices with
%   the Linear-LS method (i.e., 4th element assumed 1) which is faster than DLT for non-ideal points.
%
%   X =  X_from_xP_linLS(P,x,imsize) computes 3D point X (column 3-vector)
%   from its projections in K images x (2-by-K matrix) and camera matrices P (K-cell
%   of 3-by-4 matrices). Image sizes imsize (2-by-K matrix) are needed for preconditioning.
%
%   See also X_from_xP_lin.

% lourakis at ics.forth.gr, 2018

function X = X_from_xP_linLS(P,u,imsize)

  K=size(P,3);

  if nargin>2 && ~isempty(imsize)
    for k=1:K
      H = [2/imsize(1,k) 0 -1
           0 2/imsize(2,k) -1
           0 0              1];
      P(:,:,k)=H*P(:,:,k);
      u(:,k)=H(1:2,1:2)*u(:,k) + H(1:2,3);
    end
  end

  A=zeros(2*K, 4);
  for k=1:K
    Pk=P(:,:,k);
    A(2*(k-1)+1:2*k, :)=[u(:,k)*Pk(3,:)-Pk(1:2,:)];
  end

  if(nargin<=2 || isempty(imsize)) % no scaling specified, use diagonal scaling matrix
    D=diag(1./max(abs(A))); % A = A*D*inv(D)
    A=A*D;

    X=-A(:,1:3)\A(:, 4);
    X=(D(1:3,1:3)*X)./D(4,4); % undo normalization
  else
    X=-A(:,1:3)\A(:, 4);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_from_xP_Lp  Estimation of 3D point from image matches and camera matrices
%   assuming that the sought point is not at infinity and solving a p-norm minimization problem. 
%   If p==2 (L2 norm), the 3D point is computed with least squares and corresponds to the
%                      inhomogeneous method of HZ (i.e., Linear-LS of HS).
%   For all other values of p, the solution is obtained by solving min ||A*X-b||_p with IRLS
%
%   X =  X_from_xP_Lp(P,x,p,maxit,X0) computes 3D point X (column 3-vector)
%   from its projections in K images x (2-by-K matrix) and camera matrices P (K-cell
%   of 3-by-4 matrices).
%
%   See also X_from_xP_lin.

% lourakis at ics.forth.gr, 2018

function X = X_from_xP_Lp(P,u,p,maxit,X0)
% p, maxit & X0 are used only when p~=2

  K=size(P,3);
  Ab=zeros(2*K, 4);
  for k=1:K
    Pk=P(:,:,k);
    Ab(2*(k-1)+1:2*k, :)=[u(:,k)*Pk(3,:)-Pk(1:2,:)];
  end

  A= Ab(:, 1:3);
  b=-Ab(:, 4);

  % solve A*X=b

  if(nargin<3 || p==2) % ordinary least squares
    % SVD
    %[U,S,V]=svd(A,0);
    %c=U'*b;
    %X=V*[c(1)/S(1,1); c(2)/S(2,2); c(3)/S(3,3)];

    % QR
    %[Q,R]=qr(A,0);
    %X=R\(Q'*b);

    % LU
    [L,U,perm]=lu(A);
    X=U\(L\(perm*b));
  else
    % minimize the Lp norm ||A*X-b||_p with IRLS
    % see https://cnx.org/exports/92b90377-2b34-49e4-b26f-7fe572db78a1@12.pdf/iterative-reweighted-least-squares-12.pdf
    if nargin<4, maxit=10; end

    if nargin<5 % no initial solution, solve least squares A*X=b with QR
      [Q,R]=qr(A,0);
      X=R\(Q'*b);
    else
      X=X0;
    end
    m=size(b,1);
    wprev=Inf(m,1);

    for i=1:maxit
      e=abs(A*X - b);
      e(find(e<1E-05))=1E-05; % ignore small values
      w=e.^((p-2)/2);
      W=diag(w/sum(w)); % normalized weight matrix
      WA=W*A;
      X=(WA'*WA) \ (WA'*(W*b)); % weighted L2 solution
      %[Q, R]=qr(WA); X=R\(Q'*W*b); % as above with QR

      if(norm(w-wprev)<1E-10*m), break; end % converged

      wprev=w;
    end
    %fprintf('\nX_from_xP_Lp: IRLS terminated after %d iterations\n', i);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_from_xP_mid  Estimation of 3D point from image matches and camera matrices with
%   the midpoint method (see Ramalingam, Lodha & Sturm: A generic structure-from-motion framework)
%
%   X =  X_from_xP_mid(P,x) computes 3D point X (column 3-vector)
%   from its projections in K images x (2-by-K matrix) and camera matrices P (K-cell
%   of 3-by-4 matrices).
%
%   See also Triangulation (Ch. 7) in Szeliski, p.17 in http://www.cse.psu.edu/~rtc12/CSE486/lecture21.pdf
%   and Slabaugh et al.: Optimal Ray Intersection For Computing 3D Points From N-View Correspondences
%   This might not perform well if some of the cameras are closer to the 3D point than others.

% lourakis at ics.forth.gr, 2018

function X = X_from_xP_mid(P,u)

  K=size(P,3);

  A=zeros(3, 3);
  b=zeros(3, 1);

  for k=1:K
    Pk=P(:,:,k);
    M1=inv(Pk(1:3, 1:3));
    c=-M1*Pk(:, 4);

    uh=M1*[u(:,k); 1];
    tmp=eye(3,3) - (uh*uh')./(uh'*uh);
    A=A+tmp;
    b=b+tmp*c;
  end

  % solve A*X=b

  % SVD
  %[U,S,V]=svd(A,0);
  %c=U'*b;
  %X=V*[c(1)/S(1,1); c(2)/S(2,2); c(3)/S(3,3)];

  % QR
  %[Q,R]=qr(A,0);
  %X=R\(Q'*b);

  % LU
  [L,U,perm]=lu(A);
  X=U\(L\(perm*b));

  % Cholesky
  %R=chol(A);
  %X=R\(R'\b);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_from_xP_rob  Estimation of 3D point from image matches and camera matrices with
%   a robust method: a preliminary reconstruction is obtained by considering camera pairs
%   and selecting the midpoint of the smallest segment with a) sufficiently small length
%   with respect to the baseline and b) sufficiently large angle between the rays.
%   This point is projected on all cameras and outliers are identified from the
%   reprojection error. Finally, the point is re-estimated from all inliers.
%
%   X =  X_from_xP_rob(P,x,imsize) computes 3D point X (column 3-vector)
%   from its projections in K images x (2-by-K matrix) and camera matrices P (K-cell
%   of 3-by-4 matrices). Image sizes imsize (2-by-K matrix) are needed for preconditioning.
%
%   See also X_from_xP_mid, X_from_xP_Lp

% lourakis at ics.forth.gr, 2018

function [X inl] = X_from_xP_rob(P,u,imsize, baselnfrac, maxreprerr)

  K=size(P,3);

  % find a pair of camera matrices and image projections such as the segment perpendicular
  % to the two corresponding rays is short enough compared to their baseline

  % previous cutoff for all combinations was 12 and npairs=5*K
  if(K<=16) % generate all index pairs from the set 1..K, see http://matlabtricks.com/post-43/an-alternative-of-nchoosek-for-selecting-all-pairs
    [y, x]=find(tril(logical(ones(K)), -1));
    pairs=[x, y];
    npairs=size(pairs, 1);
  else % randomly select npairs
    npairs=7*K;
    % this guarantees uniqueness, see https://stackoverflow.com/questions/15793172/efficiently-generating-unique-pairs-of-integers
    kk=randperm(K/2*(K-1),npairs);
    y=floor(sqrt(8*(kk-1) + 1)/2 + 3/2);
    x=kk - (y-1).*(y-2)/2;
    pairs=[x;y]';
    %sortrows(pairs,[1 2])

    % see https://stackoverflow.com/questions/15684928/can-randperm-generate-several-random-permutations
    %pairs=cell2mat(arrayfun(@(x)randperm(K,2),(1:npairs)','UniformOutput',0));
  end
  
  mind=Inf;
  baselnfrac=baselnfrac*baselnfrac;
  for k=1:npairs
    i=pairs(k, 1); j=pairs(k, 2);
    
    Pi=P(:,:,i);
    Pj=P(:,:,j);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%{
    % two-view midpoint
    M0=Pi(1:3, 1:3);
    M0=inv(M0);
    c0=-M0*Pi(:, 4);

    M1=Pj(1:3, 1:3);
    M1=inv(M1);
    c1=-M1*Pj(:, 4);

    A=[M0*[u(:,i);1], -M1*[u(:,j);1]];
    b=c1-c0;
    blen=b'*b;

    [Q, R]=qr(A); a=R\(Q'*b); % LS with QR
    p0=c0+a(1)*A(:,1);
    p1=c1-a(2)*A(:,2);
    %%}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    % alternative: Lee & Civera generalized weighted midpoint ("Triangulation: Why Optimize?")
    [K0 R0]=myrq(Pi(1:3,1:3)); invK0=inv(K0); t0=invK0*Pi(:,4);
    [K1 R1]=myrq(Pj(1:3,1:3)); invK1=inv(K1); t1=invK1*Pj(:,4);
    R=R1*R0'; t=t1-R*t0; % rel. motion from Pi to Pj
    blen=t'*t;
    c0=-R0'*t0; c1=-R1'*t1;

    f0=invK0*[u(:,i);1]; nf0=f0/norm(f0);
    f1=invK1*[u(:,j);1]; nf1=f1/norm(f1);

    p=cross(R*nf0, nf1); q=cross(R*nf0, t); r=cross(nf1, t);
    magsqr=p'*p;
    % L&C formulation of midpoint
    %l0=p'*r/magsqr;
    %l1=p'*q/magsqr;

    % generalized weighted midpoint
    l0=sqrt((r'*r)/magsqr); % norm(r)/norm(p)
    l1=sqrt((q'*q)/magsqr); % norm(q)/norm(p)


    p0=t+l0*R*nf0;
    p1=l1*nf1;
    % p0, p1 are in Pj's coord. system
    p0=R1'*(p0-t1);
    p1=R1'*(p1-t1);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    seg=p1-p0; % segment perpendicular to both rays
    d=seg'*seg;

    %if(d<baselnfrac*blen), X=(p0+p1)*0.5; break; end % find first segment shorter than a fraction of the baseline

    % search for the segment of minimal length among those of acceptable length & ray angle defined by the pair indices
    if(d<mind && d<baselnfrac*blen)
      M=(p0+p1)*0.5; % midpoint

      % divide perp. segment proportionally to camera center distances from its endpoints
      %lam=norm(p1-c1)/norm(p0-c0);
      %M=p0 + (p1-p0)/(lam+1); % == p0*lam/(lam+1) + p1*1/(lam+1)

      % check cheirality
      sn=[Pi(3,:); Pj(3,:)]*[M;1];
      if(any(sn<0)), continue; end

      % check angle defined by midpoint and camera centers
      v0=M-c0; v0=v0/norm(v0);
      v1=M-c1; v1=v1/norm(v1);
      if(abs(v0'*v1)>0.984), continue; end % min acceptable ray angle ~ 10.2 deg; 0.997 would correspond to 4.4 deg, 0.998 to 3.6 deg, 0.999 to 2.5 deg, etc

      Xmin=M;
      mind=d;
    end
  end

  if(~isinf(mind))
    X=Xmin;
  else
    warning('X_from_xP_rob: no image pair suitable for triangulation found!');
    X=Inf(3, 1);
    inl=[];
    return;
  end

  % project X on each camera and compute reprojection error
  sqrepr=zeros(K,1);
  for k=1:K
    Pk=P(:,:,k);
    m=Pk*[X;1]; m=m(1:2)/m(3);
    d=m-u(:,k);
    sqrepr(k)=d'*d;
  end

  sqthresh=maxreprerr*maxreprerr;
  if maxreprerr>0
    inl=find(sqrepr<=sqthresh);
    if(numel(inl)<=2 && K>2) % thresholding did not find sufficient inliers
      %{
      % adjust to retain 0.75*K smallest repr. errors
      maxreprerr=-75;
      warning('X_from_xP_rob: too strict threshold, adjusting to keep %.0f%% of camera-point pairs', maxreprerr);
      %}

      % alternative: find inliers with the X84 rejection rule
      repr=sqrt(sqrepr);
      loc=median(repr);
      repr2=abs(repr-loc);
      scl=1.4826*2.5*median(repr2); % accept upto 2.5 standard deviations (3.7 MADs)
      inl=find(repr2<=scl);

      if(length(inl)<2) % return the two smallest
        [dummy idx]=sort(repr2);
        inl=idx(1:2);
      end
      warning('X_from_xP_rob: too strict threshold, invoking the X84 rule');
    end 
  end 

  if maxreprerr<0 % interpret as percent, e.g. -70 implies retaining 0.7*K smallest repr. errors. Note: might be reached with -75 from above!
    maxreprerr=min(-maxreprerr, 100);

    srtd=sort(sqrepr);
    sqthresh=srtd(uint32((maxreprerr*K)/100));

    if sqthresh>10*10 % sanity check
      warning('X_from_xP_rob: unexpectedly large reprojection threshold %f', sqrt(sqthresh));
    end

    inl=find(sqrepr<=sqthresh);
    %if isempty(inl), error('X_from_xP_rob: no inliers found! (threshold %f)', sqrt(sqthresh)); end % should not happen
  end

  P=P(:,:,inl);
  u=u(:,inl);
  if nargin>2 && ~isempty(imsize)
    imsize=imsize(:,inl);
  end

  X=X_from_xP_Lp(P, u, 1.5, 15, X); % triangulate using all inliers
  %X=X_from_xP_linLS(P, u, imsize);

  %if (X(3)<0), X=Inf(3, 1); end
end

% RQ decomposition of 3x3 matrix
% see https://math.stackexchange.com/questions/1640695/rq-decomposition
% Kovesi's rq3 should be faster
function [R,Q]=myrq(A)
  P=[0 0 1; 0 1 0; 1 0 0];
  [Q,R]=qr((P*A)');
  R=P*R'*P;
  Q=P*Q';

  % make sure that the diagonal elements of R are positive
  for n=1:3
    if R(n,n)<0
      R(:,n)=-R(:,n);
      Q(n,:)=-Q(n,:);
    end
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Estimation of 3D point from image matches and camera matrices, nonlinear.
%  X = X_from_xP_nonlin_levmar(x,P,imsize,X0) computes max. likelihood estimate of projective
%  3D point X (column 4-vector) from its projections in K images x (2-by-K matrix)
%  and camera matrices P (K-cell of 3-by-4 matrices). Image sizes imsize (2-by-K matrix)
%  are needed for preconditioning. X0 it the initial estimate of X.
%  By minimizing reprojection error, LM iterations.
%
%  See also vgg_X_from_xP_nonlin.

% lourakis at ics.forth.gr, 2018

function X = X_from_xP_nonlin_levmar(u,P,imsize,X)

  if iscell(P)
    P = cat(3,P{:});
  end
  K = size(P,3);

  % precondition
  if ~isempty(imsize)
    for k = 1:K
      H = [2/imsize(1,k) 0 -1
           0 2/imsize(2,k) -1
           0 0              1];
      P(:,:,k) = H*P(:,:,k);
      u(:,k) = H(1:2,1:2)*u(:,k) + H(1:2,3);
    end
  end

  % Parametrize X such that X = T*[Y;1]; thus x = P*T*[Y;1] = Q*[Y;1]
  [dummy,dummy,T] = svd(X',0);
  T = T(:,[2:end 1]);
  for k = 1:K
    Q(:,:,k) = P(:,:,k)*T;
  end

  % Levenberg-Marquardt
  Y = [0;0;0];
  [Y info]=marquardt(@objfunc_nonlin, {u, Q}, Y, [1E-04 1E-14 1E-14 30]);
  %fprintf('\nX_from_xP_nonlin_levmar: LM terminated after %d iterations, reason %d\n', info(5), info(6));

  %{
  % Newton
  Y = [0;0;0];
  eprev = inf;
  for n = 1:10
    [e,J] = objfunc_nonlin(Y,{u,Q});
    if 1-norm(e)/norm(eprev) < 1000*eps
      break;
    end
    eprev = e;  
    Y = Y - (J'*J)\(J'*e);
  end
  %}

  X = T*[Y;1];

  return
end

function [e, J] = objfunc_nonlin(Y, parms)
  u=parms{1};
  Q=parms{2};

  K = size(Q,3);
  e = [];
  J = [];
  %e = zeros(2*K,1);
  %J = zeros(2*K,3);
  for k = 1:K
    q = Q(:,1:3,k);
    x0 = Q(:,4,k);
    x = q*Y + x0;
    e = [e; x(1:2)/x(3)-u(:,k)];
    J = [J; [x(3)*q(1,:)-x(1)*q(3,:)
             x(3)*q(2,:)-x(2)*q(3,:)]/x(3)^2];

    %e(2*k-1:2*k)   = [x(1:2)/x(3)-u(:,k)];
    %J(2*k-1:2*k,:) = [x(3)*q(1,:)-x(1)*q(3,:);   x(3)*q(2,:)-x(2)*q(3,:)]/x(3)^2;
  end
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X_from_xP_ang  Estimation of 3D point from image matches and camera matrices with
%   angular triagulation (see Recker et al.: "Statistical Angular Error-Based Triangulation
%   for Efficient and Accurate Multi-View Scene Reconstruction"
%
%   X =  X_from_xP_ang(P,x) computes 3D point X (column 3-vector)
%   from its projections in K images x (2-by-K matrix) and camera matrices P (K-cell
%   of 3-by-4 matrices).
%

% lourakis at ics.forth.gr, 2018

function X = X_from_xP_ang(P,u,X)

  K=size(P,3);

  % starting point
  if(nargin<=2), X=X_from_xP_mid(P,u); end

  % Levenberg-Marquardt
  [X info]=marquardt(@objfunc_ang, {P, u}, X, [1E-04 1E-14 1E-14 50]);
  %fprintf('\nX_from_xP_ang: LM terminated after %d iterations, reason %d\n', info(5), info(6));
  return;

  %{
  % Gauss-Newton; might fail to converge correctly
  e=zeros(K, 1);
  J=zeros(K, 3);
  prevnrm=inf;
  for j=1:20
    [e, J]=objfunc_ang(X, {P, u});

    nrm=norm(e);
    if(1-nrm/prevnrm<1000*eps)
      break;
    end
    prevnrm=nrm;

    % d=(J'*J)\(J'*e); % normal eqs
    %[Q, R]=qr(J); d=R\(Q'*e); % QR
    [U S V]=svd(J,0);
    if(S(1,1)<S(3,3)*1000) % well-conditioned J, solve with its SVD
      d=V*diag(1./[S(1,1) S(2,2) S(3,3)])*U'*e;
    else
      % solve with truncated SVD
      d=V(:,1:2)*diag(1./[S(1,1) S(2,2)])*U(:,1:2)'*e;

      %{
      % Tikhonov regularization
      svals=diag(S);
      delta=.001; % regularization parameter
      Sdelta=diag(svals./(svals.^2+delta));
      d=V*Sdelta*U'*e;
      %}
    end
    X=X - d;
  end
  %fprintf('\nX_from_xP_ang: GN terminated after %d iterations\n', j);
  %}
end


% compute eq.(3) and its Jacobian in e & J
function [e, J] = objfunc_ang(X, parms)
  P=parms{1};
  u=parms{2};
  K=size(P,3);
  e=zeros(K, 1);
  J=zeros(K, 3);
    for k=1:K
      Pk=P(:,:,k);
      M1=inv(Pk(1:3, 1:3));
      c=-M1*Pk(:, 4);

      w=M1*[u(:,k); 1];
      w=w./norm(w);

      v=X-c;
      vmag=norm(v);
      vdotw=v'*w;

      J(k,:)=-w'./vmag + v'.*(vdotw/(vmag^3));
      e(k)=1-vdotw/vmag;
    end
end


%--------------------------------------------------------------------------
% finite difference Jacobian for verification
function J = fjac_ang(X, parms)
  P=parms{1};
  u=parms{2};

  K=size(P,3);

  delta=1e-6*(max(1,sqrt(norm(X))));
  m=3;
  J=zeros(K, m);
  dX=zeros(m, 1);
  for k=1:K
      Pk=P(:,:,k);
      M1=inv(Pk(1:3, 1:3));
      c=-M1*Pk(:, 4);

      w=M1*[u(:,k); 1];
      w=w./norm(w);

      v=X-c;
      %e=1-(v'*w)/norm(v);
    for j=1:m
      dX(j)=delta;
      v=X+dX-c; ep=1-(v'*w)/norm(v);
      v=X-dX-c; em=1-(v'*w)/norm(v);
      %J(k,j)=(ep - e)/delta; % fwd
      J(k,j)=(ep - em)/(2*delta); % central
      dX(j)=0;
    end
  end
end
