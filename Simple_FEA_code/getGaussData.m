function [ weight, gauss_points ] = getGaussData( n_gp )
% a function returns wieght and gauss points

if n_gp == 2
    weight = [1, 1];
    gauss_points = [1/sqrt(3), -1/sqrt(3)];
elseif n_gp == 3
    weight = [5/9, 8/9, 5/9];
    gauss_points=[-sqrt(3/5), 0, sqrt(3/5)];
elseif n_gp == 4
    weight = [ (18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36 ];
    gauss_points = [-sqrt( (15+2*sqrt(30))/35 ), -sqrt( (15-2*sqrt(30))/35 ), sqrt( (15-2*sqrt(30))/35 ), sqrt( (15+2*sqrt(30))/35 )];
elseif n_gp == 5
    weight = [ 128/225, (322+13*sqrt(70))/900, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900,(322-13*sqrt(70))/900, ];
    gauss_points = [0, 1/3 * sqrt(5-2*sqrt(10/7)) , -1/3 * sqrt(5-2*sqrt(10/7)), 1/3 * sqrt(5+2*sqrt(10/7)), -1/3 * sqrt(5+2*sqrt(10/7))];
end
end

