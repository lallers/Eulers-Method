function du= euler3(t,u)
%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
du = zeros(2,1);
du(1)= u(2);
du(2)= 7.5*cos(t)-0.05*u(2) - u(1)^3;
end