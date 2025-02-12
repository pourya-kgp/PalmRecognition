% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : CircularGaborFilter.m
% Description   : Applies a Circular Gabor Filter (with zero DC component) to an input image
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

function [reGabor, imGabor] = CircularGaborFilter(image, filSz, theta, u, sigma)

% Inputs:
% - image   : Input image
% - filSz   : Size of the filter
% - theta   : Orientation of the function (in radians)
% - u       : Frequency of the sinusoidal wave
% - sigma   : Standard deviation of the Gaussian envelope
%
% Outputs:
% - reGabor : Real part ot filtered image
% - imGabor : Imaginary part ot filtered image
%
% Circular Gabor filter description:
% The circular Gabor filter is an effective tool for texture analysis and has the following
% general form:
%                   1              -(x^2 + y^2)
% Gabor(x,y) = ------------ * exp (------------ + 2*pi*i*(u*x*cos(theta) + u*y*sin(theta)))
%              2*pi*sigma^2          2*sigma^2
% Where
% - i     : sqrt(-1), imaginary unit
% - u     : Frequency of the sinusoidal wave
% - theta : Orientation of the wave (in radians)
% - sigma : Standard deviation of the Gaussian envelope.
%
% To make it more robust against brightness, a discrete Gabor filter is turned to zero DC
% (direct current) with the application of the following formula:
% Gab_DC0 = Gab - sum(Gab)/filSz^2
%
% Reference:
% - Zhang, D., Kong, W., You, J., & Wong, M. (2003). "Online Palmprint Identification". 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence. 
 
% -------------------- Adjusting the parameters
% Ensuring the image is in double format for calculations
if ~isa(image, 'double')
    image = double(image);
end

% Ensuring the filter size is odd
if mod(filSz, 2) == 0
    filSz = filSz + 1;
end

% -------------------- Generating the filter           
% Generating a grid of coordinates for the filter
[x, y] = meshgrid(-fix(filSz/2):fix(filSz/2), fix(filSz/2):-1:fix(-filSz/2));

% Defining the Gabor filter
Gab = 1/(2*pi*sigma^2)*exp(-.5*(x.^2+ y.^2)/sigma^2 + 2*pi*1i*(u*x*cos(theta)+u*y*sin(theta)));

% Adjusting the Gabor filter to zero DC
DC0 = sum(Gab(:))/filSz^2;
Gab_DC0 = Gab - real(DC0);

% -------------------- Filtering the input image
% Apply the Gabor filter to the image
reGabor = conv2(image, double(real(Gab_DC0)), 'valid'); % Real part
imGabor = conv2(image, double(imag(Gab_DC0)), 'valid'); % Imaginary part

end