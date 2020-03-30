% multimg_matches - build a set of multi-image feature matches (i.e., matches in N>2 views forming trajectories)
%                   from a set of pairwise matches: if (i, j) and (j, k) are pairwise matches,
%                   then the multi-image match includes i, j, k, etc
%
%    mims = multimg_matches(pmatches)
%
%    pmatches  - cell array of cell arrays such that pmatches{i}{j} is a 2 x M matrix of feature
%                indices between frames i,j. Note that pmatches is upper triangular, i.e. cells
%                with i>=j are empty
%
% Returns
%    mims      - N x nframes array where each row corresponds to one multi-image match and nonnegative
%                elements are feature indices in the frames corresponding to their columns

% Constructs a directed graph whose vertices are fetures and its edges correspond to pairwise matches
% and finds multi-image matches as connected components with BFS


% Manolis Lourakis 2018
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece

% October  2018  - Original version. (v. 1.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mims] = multimg_matches(pmatches)

  if(~iscell(pmatches))
    error('multimg_matches: expected a cell array for argument "pmatches"!');
  end

  nframes=length(pmatches)+1; % #frames

  % count the pairwise matches and find the matched feature ids in every frame
  npmatches=zeros(nframes, nframes);
  feats=cell(nframes, 1);
  for i=1:nframes-1
    for j=i+1:nframes
      pmatchesij=pmatches{i}{j};
      if isempty(pmatchesij), continue; end
      npmatches(i, j)=size(pmatchesij, 2);
      feats{i}=union(feats{i}, pmatchesij(1,:));
      feats{j}=union(feats{j}, pmatchesij(2,:));
    end
  end

  nfeat=zeros(nframes, 1); % #features per frame
  for i=1:nframes, nfeat(i)=length(feats{i}); end

  frstart=[0; cumsum(nfeat)]; % offset of each frame in ftinfr (see below); size nframes+1
  totfeat=frstart(nframes+1);

  % ftinfr maps feature numbers in 1..totfeat to feature ids in individual frames
  % the inverse mapping, i.e. feature ft in frame i is binsearch(ftinfr, frstart(i), frstart(i+1), ft)
  ftinfr=zeros(totfeat, 1);
  for k=1:nframes
    for i=1:nfeat(k)
      ftinfr(frstart(k)+i)=feats{k}(i);
    end
  end

tic;
  totpmatches=sum(npmatches(:)); % #pmatches in all pairs
  G=sparse([],[],[], totfeat, totfeat, totpmatches);
  for i=1:nframes-1
    fprintf('. '); % show progress
    for j=i+1:nframes
      pmatchesij=pmatches{i}{j};
      for k=1:length(pmatchesij)
        m=pmatchesij(:, k);
        ii=binsearch(ftinfr, frstart(i), frstart(i+1), m(1));
        jj=binsearch(ftinfr, frstart(j), frstart(j+1), m(2));

        if(ii==-1 || jj==-1), error('multimg_matches: internal error!'); end % should not happen

        G(ii, jj)=1; % G is directed!
      end
    end
  end
  fprintf('\n');
toc

  fprintf('multimg_matches: found %d matched features, matches graph has %d edges (density %.3f%%)\n', totfeat, totpmatches, nnz(G)*100/totfeat^2);

tic;
  % BFS for connected components
  [rp, ci]=sparse_to_csr(G+G'); % adding G' makes the CRS matrix undirected!
  visited=zeros(totfeat, 1);
  ncc=0;
  mims=zeros(floor(totfeat/10), nframes); % arbitrary size!

  for v=1:totfeat
    if(~visited(v))
      d=bfs_csr(rp, ci, v);
      cc=find(d>-1);
      visited(cc)=1;
      if(length(cc)<2), continue; end;
      fids=ftinfr(cc);
      mim=-1*ones(1, nframes);
      for i=1:length(fids)
        fr=sum(frstart<cc(i));
        mim(fr)=fids(i);
      end
      ncc=ncc+1;
      mims(ncc,:)=mim;
    end
  end

  mims(ncc+1:end, :)=[]; % remove unused rows
toc

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% binary search for val in mat(left:right)
function [idx] = binsearch(mat, left, right, val)

  idx=-1;

  while left<=right
    mid=ceil((left+right)/2);

    d=mat(mid)-val;
    if d>0
      right=mid-1;
    else if d<0
          left=mid+1;
        else % d==0
          idx=mid;
          return;
      end
    end
  end

end

% bfs.m & sparse_to_csr.m are from gaimc at https://www.mathworks.com/matlabcentral/fileexchange/24134-gaimc-graph-algorithms-in-matlab-code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d] = bfs_csr(rp,ci,u)
% BFS Compute breadth first search distances for a graph
%
% [d] = bfs_csr(rp,ci,u) returns the distance (d)
% for each vertex in the graph in a breadth first search
% starting from vertex u.
%   d = -1 if vertex i is not reachable from u
%   rp and ci are the row pointer & column index for the graph in compressed row storage (see sparse_to_csr)
% (COMMENTED OUT --ML) pred is the predecessor array.  pred(i) = 0 if vertex (i)
% is in a component not reachable from u and i != u.
%
%
% Example:
%   load_gaimc_graph('bfs_example.mat') % use the dfs example from Boost
%   [rp ci]=sparse_to_csr(A);
%   d = bfs_csr(rp,ci,1)
%
% See also DFS

% Manolis Lourakis, Oct. 2018
% Based on bfs.m by David F. Gleich
% Copyright, Stanford University, 2008-20098

% History
% 2018-10-11: Initial coding

n=length(rp)-1;
d=-1*ones(n,1);
%pred=zeros(1,n);
sq=zeros(n,1); sqt=0; sqh=0; % search queue and search queue tail/head

% start bfs at u
sqt=sqt+1; sq(sqt)=u;
d(u)=0;
%pred(u)=u;
while sqt-sqh>0
    sqh=sqh+1; v=sq(sqh); % pop v off the head of the queue
    for ri=rp(v):rp(v+1)-1
        w=ci(ri);
        if d(w)<0
            sqt=sqt+1; sq(sqt)=w;
            d(w)=d(v)+1;
            %pred(w)=v;
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rp ci ai ncol]=sparse_to_csr(A,varargin)
% SPARSE_TO_CSR Convert a sparse matrix into compressed row storage arrays
% 
% [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp), column index
% (ci) and value index (ai) arrays of a compressed sparse representation of
% the matrix A.
%
% [rp ci ai] = sparse_to_csr(i,j,v,n) returns a csr representation of the
% index sets i,j,v with n rows.
%
% Example:
%   A=sparse(6,6); A(1,1)=5; A(1,5)=2; A(2,3)=-1; A(4,1)=1; A(5,6)=1; 
%   [rp ci ai]=sparse_to_csr(A)
%
% See also CSR_TO_SPARSE, SPARSE

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-07: Initial version
% 2008-04-24: Added triple array input
% 2009-05-01: Added ncol output
% 2009-05-15: Fixed triplet input

error(nargchk(1, 5, nargin, 'struct'))
retc = nargout>1; reta = nargout>2;

if nargin>1
    if nargin>4, ncol = varargin{4}; end
    nzi = A; nzj = varargin{1};
    if reta && length(varargin) > 2, nzv = varargin{2}; end    
    if nargin<4, n=max(nzi); else n=varargin{3}; end
    nz = length(A);
    if length(nzi) ~= length(nzj), error('gaimc:invalidInput',...
            'length of nzi (%i) not equal to length of nzj (%i)', nz, ...
            length(nzj)); 
    end
    if reta && length(varargin) < 3, error('gaimc:invalidInput',...
            'no value array passed for triplet input, see usage'); 
    end
    if ~isscalar(n), error('gaimc:invalidInput',...
            ['the 4th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
    if nargin < 5, ncol = max(nzj); 
    elseif ~isscalar(ncol), error('gaimc:invalidInput',...
            ['the 5th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
else
    n = size(A,1); nz = nnz(A); ncol = size(A,2);
    retc = nargout>1; reta = nargout>2;
    if reta,     [nzi nzj nzv] = find(A); 
    else         [nzi nzj] = find(A);
    end
end
if retc, ci = zeros(nz,1); end
if reta, ai = zeros(nz,1); end
rp = zeros(n+1,1);
for i=1:nz
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
if ~retc && ~reta, rp=rp+1; return; end
for i=1:nz
    if reta, ai(rp(nzi(i))+1)=nzv(i); end
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
