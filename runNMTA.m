% Copyright (C) 2017 Andreas Bertsatos <andreas.bertsatos@gmail.com>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.
%
%
% This is an example script for running a non-metric trait analysis based on
% a number of population samples. The functions nmtanalysis.m, checkNMTA.m,
% smithsMMD.m, cellrange.m and fishextest.m are required as well as the 
% io-package needs to be installed in GNU Octave.
%
%
% load external library .io for manipulating odf files
pkg load io;
% define the name of each population to be included in the analysis
% the .csv file extension is added automatically
population = {'population_A','population_B','population_C',...
              'population_D','population_E','population_F'};
% calculate frequency, laterality, sexual dimorphism
% and inter trait correlation for each population
% and save resultfs in a single structure container
for p=1:length(population)
  a(p) = nmtanalysis(char(population(p)),5);
endfor

%define threshold of p-value for laterality test
p_lat = 0.05;
%define threshold of p-value for sexual dimorphism test
p_sex = 0.05;
%define threshold of p-value for inter trait correlation test
p_int = 0.01;


freq_table = checkNMTA(a,population,p_lat,p_sex,p_int);
% calculate Smith's MMD from the resulted frequency table
MMD = smithsMMD(freq_table,population);

disp('');disp('Resulted Smith''s MMD variance matrix for populations.');
disp('');disp(MMD.variance);
