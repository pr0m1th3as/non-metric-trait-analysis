% Copyright (C) 2016 Andreas Bertsatos <andreas.bertsatos@gmail.com>
% Copyright (C) 2016 Maria-Eleni Chovalopoulou <chovalopoulou.eleni@gmail.com>
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
function MMD = smithsMMD(varargin)
% function MMD = smithsMMD(freq_table, population, 'angular transformation')
%
% This function takes the frequency table for a number of traits that
% correspond to the populations contained in the 'population' string cell
% and calculates the Smith's MMD incorporating the Freeman and Tukey angular 
% transformation. Except for the MMD distance matrix, which is stored
% in 'MMD.distance' structure, the function also calculates the variance
% matrix, stored in 'MMD.variance'. Subsequently, the "standardized distances"
% are calculated and stored in 'MMD.significance', which provide a significance
% measure for the calculated distances based on that when the hypothesis of
% equality of proportions holds, that is, when the populations are not divergent,
% the MMD may be regarded as significant at a significance level of approximately
% 2.5% when larger than twice the standard deviation (Sjøvold, 1977). Since the
% distance matrix produced by MMD is likely to be subsequently used as input for
% cluster analysis (hence the phylip format), another matrix is calculated and 
% stored under 'MMD.phylip', where all negative distances are set to zero.
% Furthermore, all MMDs that are less than twice their standard deviations are
% set to zero since these estimates of the underlying population differences
% are statistically nonsignificant (Harris and Sjøvold, 2004).
%  
% The distance matrix, stored in 'MMD.phylip', is printed in the Phylip format
% for direct use with the online T-REX web server (Boc. et al., 2012) at
% http://www.trex.uqam.ca
% The original distance matrix, stored in 'MMD.distance', is also printed in
% the Phylip format for subsequent used if desired. All four matrices are
% returned in a single structure container called 'MMD'.
%
% The MMD calculation follows the formula with the {(1/n_i+0.5)+(1/n_j+0.5)}
% correction term (attributable to Freeman and Tukey, 1950) as it performs
% better at extreme trait frequencies (Green and Suchey, 1976). The Freeman
% and Tukey angular transformation is used by default in this function.
% However, the Anscombe (1948) angular transformation is also available for
% calculating theta values if selected in the input arguments. The results are
% essentially identical, but Anscombe transformation may be preferable under
% certain circumstances (Harris and Sjøvold, 2004).
% 
% The input array must be 3 dimensional with each corresponding
% dimension as follows:
%
% freq_table(population, trait, frequency)
%
% where population and trait are scalars and frequency is a 1x2 matrix
% containing number of trait occurence and number of observations for
% given trait and population.
%
% The string cell must contrain the names of the populations used in the analysis
% in the same order as they appear in the frequency table.
%
% The angular transformation can be defined as a character string. The two
% options are 'Freeman and Tukey' or 'Anscombe' and can be parsed into the 
% fuction as an optional argument. If no argument for angular transformation
% is defined by the user, Freeman and Tukey is used by default.
%
% References
%
% Anscombe FJ. 1948. The transformation of Poisson, binomial and negative-binomial
%   data. Biometrika 35:246-254.
% Freeman MF, Tukey JW. 1950. Transformations related to the angular and square
%   root. Ann Math Stat 21: 607-611.
% Green R.F., Suchey J.M. (1976), The use of inverse sine transformations in
%   the analysis of non-metric cranial data, Am J Phys Anthropol 45:61-68
% Harris E.F., Sjøvold T. (2004), Calculation of Smith’s Mean Measure of Divergence
%   for intergroup comparisons using nonmetric data, Dental Anthropology 17(3):83-93.
% Sjøvold T. (1977), Non-metrical divergence between skeletal populations,
%   Ossa 4:suppl. 1

% Check the number of input arguments and parse them into the appropriate
% variables. If no angular transformation is defined then the Freeman & Tukey
% will be used by default. Additional checks are made on whether the number of
% populations is consistent between frequency table and population list or
% whether the angular transformation has been identified correctly.

if (length(varargin)<2)
  printf('\nInsufficient input arguments.\n');
  printf('Both frequency table and population list are required.\n');
  MMD = [];
  return;
elseif (length(varargin{1}(:,1,1))!=length(varargin{2}))
  printf('\nFrequency table and population list do not match.\n');
  MMD = [];
  return;
elseif (length(varargin)==2) && (length(varargin{1}(:,1,1))==length(varargin{2}))
  freq_table = varargin{1};
  population = varargin{2};
  thetatrans = 'Freeman and Tukey';
  printf('\nFreeman and Tukey angular transformation is used by default.\n');
elseif (length(varargin)==3) && ((strcmp(varargin{3},'Freeman and Tukey')) ...
       || (strcmp(varargin{3},'Anscombe')))
  freq_table = varargin{1};
  population = varargin{2};
  thetatrans = varargin{3};
  printf('\nUser defined %s angular transformation is used.\n',thetatrans);
elseif (length(varargin)==3) && !((strcmp(varargin{3},"Freeman and Tukey")) ...
       || (strcmp(varargin{3},'Anscombe')))
  printf('\nUndefined angular transformation.\n');
  MMD = [];
  return;
endif
% Calculate MMD distance and variance matrices.
for p1=1:length(freq_table(:,1,1))-1
  for p2 = p1+1:length(freq_table(:,1,1))
    distance(p1,p2-1,:) = 0;
    variance(p1,p2-1,:) = 0;
      for i=1:length(freq_table(1,:,1))
      k1 = freq_table(p1,i,1);
      n1 = freq_table(p1,i,2);
      k2 = freq_table(p2,i,1);
      n2 = freq_table(p2,i,2);
      if (strcmp(thetatrans,'Freeman and Tukey'))
        theta1 = 0.5*asin(1-(2*k1)/(n1+1)) + 0.5*asin(1-2*(k1+1)/(n1+1));
        theta2 = 0.5*asin(1-(2*k2)/(n2+1)) + 0.5*asin(1-2*(k2+1)/(n2+1));
      elseif (strcmp(thetatrans,'Anscombe'))
        theta1 = asin(1-(2*((k1+0.375)/(n1+0.75))));
        theta2 = asin(1-(2*((k2+0.375)/(n2+0.75))));
      endif
      traitMMD = ((theta1 - theta2)^2)-((1/(n1+0.5))+(1/(n2+0.5)));
      distance(p1,p2-1,:) = distance(p1,p2-1,:) + traitMMD;
      variance(p1,p2-1,:) = variance(p1,p2-1,:) + ((1/(n1+0.5))+(1/(n2+0.5)))^2;
      endfor
    distance(p1,p2-1,:) = (distance(p1,p2-1,:))/length(freq_table(1,:,1));
    variance(p1,p2-1,:) = (2*(variance(p1,p2-1,:)))/((length(freq_table(1,:,1)))^2);
  endfor
endfor

% Store the calculated distance and variance in MMD structure.
MMD.distance = distance;
MMD.variance = variance;

% Duplicate the upper triangular matrix to its mirror while leaving all
% diagonal elements zero.
MMD.distance = [zeros(length(MMD.distance),1) MMD.distance];
MMD.distance = [MMD.distance;zeros(1,length(MMD.distance))];
MMD.distance = triu(MMD.distance)+triu(MMD.distance,1)';
% Make similar rearrangement to the MMD.variance matrix and compute the 
% significance MMD matrix for examining statistical significance at a
% level of ~2.5% for values greater than 2.
MMD.variance = [zeros(length(MMD.variance),1) MMD.variance];
MMD.variance = [MMD.variance;zeros(1,length(MMD.variance))];
MMD.variance = triu(MMD.variance)+triu(MMD.variance,1)';
% Calculate the significance MMD matrix.
MMD.significance = MMD.distance./(sqrt(MMD.variance));
% Correct the NaN values due to division by zero.
MMD.significance(isnan(MMD.significance))=0;
% Parse distances from MMD.distance to MMD.phylip, while setting to zero all
% negative or nonsignificant distances.
MMD.phylip = MMD.distance;
MMD.phylip(MMD.distance<0)=0;
MMD.phylip(MMD.significance<2)=0;

% Check number of populations and save their names in a string cell.
p_number = length(population);
for i=1:p_number
  p_name{i} = char(population(i));
endfor
% Find the longest name in the list
for i=1:p_number
  namelen(i) = length(p_name{i});
endfor
maxlen = max(namelen);
% and pad space letters as required for all population names to be the 
% same length as max length while preventing octave to issue a warning
% related to cstrcat function that does apply to this case.
for i=1:p_number
  if (length(p_name{i})<maxlen)  
    p_name{i} = cstrcat(p_name{i},blanks(maxlen-length(p_name{i})));
  endif
endfor

% print the MMD distance matrix in the phylip format with only the
% distances that are statistically significant
printf('\nThe statistically significant MMD distance matrix in the Phylip format.\n');
printf('%i\n',p_number);
for i=1:p_number
  format short;
  printf('%s%s',p_name{i},char(disp(MMD.phylip(i,:))));
endfor

% print the original MMD matrix in the phylip format with the initially
% calculated distance values
printf('\nThe original MMD distance matrix in the Phylip format.\n');
printf('%i\n',p_number);
for i=1:p_number
  format short;
  printf('%s%s',p_name{i},char(disp(MMD.distance(i,:))));
endfor

endfunction
