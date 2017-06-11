#ifndef ANGLES_H
#define ANGLES_H

enum AngleShape {ConstantShape, EllipticalShape, InverseEllipticalShape};

/**
function [theta, kmax, dcf_theta] = calc_angles(theta0, theta_width, FOV, varargin)
% [theta, kmax, dcf_theta] = calc_angles([theta0, theta_width] FOV, F1, F2, ...,
%                                        [KFCN, K1, K2, ...])
%
% Calculates a set of angles and k-space extents for desired
% radial imaging field-of-view (FOV)
% Some FOV/KFCN functions are provided in the shape_fcns directory.
% Use "help shape_fcns" for more information on these functions.
%
% Inputs:
%   theta0 (optional) - initial starting angle, in radians.  Default is 0.
%   theta_width (optional, requires theta0 when used) - size of range of
%         angles to design.  Used in 3D PR design.  Default is pi.
%   FOV - function handle of the desired FOV (@fcn_name).
%         This function must be pi-periodic and convex
%   F1, F2, ... - Input parameters to FOV function
%   KFCN (optional) - function of desired kmax (defaults to constant)
%   K1, K2, ... - Inputs to KFCN function, see comments below on values here
%
% Outputs:
%   theta - resulting angles for desired FOV/KFCN
%           angles are for a full-spoke (0 < theta < pi)
%   kmax - trajectory extents for corresponding theta's
%   dcf_theta - angular density compensation weighting function
%
% For FOV, the inputs to the external functions correspond to the pixel sizes.
% For KFCN, they correspond to the inverse of the resolution pixel size.
%  (KFCN = @const, K1 = 0.5 will give a resolution size of 2 pixels,
%   while K1 = 1 will give a resolution of 1 pixel)
%
% Examples:
%   % Add shape_fcns directory to the path
%   addpath([pwd '/shape_fcns'])
%
%   % Design a circular FOV
%   X = 80;
%   [theta, kmax, dcf_theta] = calc_angles(@const, X);
%   % theta will consist of X * pi / 2 equally spaced angles
%
%   % Design a 2D elliptical FOV
%   X = 80; Y = 120;
%   [theta, kmax, dcf_theta] = calc_angles(@ellipse, X, Y);
%
%   % Rectangular FOV with rectangular kmax pattern
%   [theta, kmax, dcf_theta] = calc_angles(@rect, X, Y, @rect, 1, X/Y);
%   % NOTE: resolution will be 1 pixel in x, and Y/X in y
%
%   % Oval FOV with initial angle of pi/2
%   theta0 = pi/2; R = 40;
%   [theta, kmax, dcf_theta] = calc_angles(theta0, @oval, X, Y, R);
%
%   % Diamond FOV using half projection acquisitions
%   [theta, kmax, dcf_theta] = calc_angles(0, 2*pi, @diamond, X, Y);
%
% Paul Gurney and Peder Larson, 12/12/2005, updated 5/30/2006
% (c) 2006, Board of Trustees, Leland Stanford Junior University
*/
void calculateAngles(float theta0, float theta_width, enum AngleShape FOV, float *F, enum AngleShape KFCN, float *K, float **theta, float **kmax, float **dcf_theta, int *ntheta);

void calculateAngles3D(int halfProjection, enum AngleShape thetaShape, enum AngleShape phiShape, int *N, float **thetaOut, float **phiOut, float **kmaxOut, float **dcfOut, int *Nproj);

float getExtent(enum AngleShape shape, float angle, float *axisLengths);

#endif // ANGLES_H
