function [z,angle,step] = translation(z0,angle0,multiplier,num)
% The function is to apply tranlation operation on a data matrix which has
% periodic boundary and generate a seqence of matrix after different
% tranlation steps.
%===========================================================
% The input parameters:
% z0: the input data of surface 
% angle0: The angle of the chosen direction in degree
% multiplier: unit step * multiplier = the chosen step size
% num: the number of steps
%===========================================================
% The output results:
% z: the seqence of matrix after different steps of translations
% angel: the angle with reasonable step size which is closest to the chosen one
% step: the step size 

    z = cell(1,num);
    [elex,eley] = size(z0);
    x = round(cos(angle0 / 180 * pi)*100); % turn the input angle into the  
    y = round(sin(angle0 / 180 * pi)*100); % integer coordinates on the circle 
    l = gcd(x,y); % with radius 100
    x = x/l;
    y = y/l; % find the smallest step by being divided by gcd
    if x*multiplier*num > elex || y*multiplier*num > eley 
        disp('The smallest step of translation for the given angle is too large compared to the size of array.')
        return;
    end
    step = norm(x*multiplier,y*multiplier);
    angle = atan(y/x)/pi *180;
    for i = 1:num
        z0 = circshift(z0,x*multiplier*i,2); % use circshift() to tranlate array
        z0 = circshift(z0,-y*multiplier*i,1);
        z{1,i} = z0;
    end
end
