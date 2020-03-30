function [descriptorsOut] = rootSift(descriptorsIn)
%
x = descriptorsIn;

x=x'; 
% y=y'; % transpose for faster access
for i=1:size(x,2),  x(:,i)=sqrt(x(:,i)/norm(x(:,i),1)); end
% for i=1:size(y,2),  y(:,i)=sqrt(y(:,i)/norm(y(:,i),1)); end
x=x';

% y=y';

descriptorsOut = x;

end

