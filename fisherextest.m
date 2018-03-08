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
function P2_tail = fisherextest(a,b,c,d)
% P2_tail=fisherextest(a,b,c,d)
%
% This function performs the Fisher exact probability test for a table of 
% frequencies, which correspond to a (2x2) contigency table. Calling the
% the function with an output argument will return the 2-tail P-value
% of the Fisher's exact test. When called with no output arguments it will
% display the result as in the following example
%
%       >> fisherextest(5,2,9,15);
%
%       Fisher's exact test for frequency cells:
%       -------------------------------------------
%       a = 5, b = 2
%       c = 9, d = 15
%       -------------------------------------------
%       2-tail P-value: 0.1975478927
%
%
%
% Agresti, A. (1992), A Survey of Exact Inference for Contegency Tables. 
%          Statistical Science,7:131-153.
% Feldman, S.E. (1963), Short cut calculation of the Fisher-Yates "exact test".
%          Psycometrika, 28(3):289-291.
% Fisher, R.A. (1934), Statistical Methods for Research Workers. Chapter 12. 
%          5th Ed., Scotland:Oliver & Boyd.
% Howell, I.P.S. (Internet homepage), http://www.fiu.edu/~howellip/Fisher.pdf
% Zar, J.H. (1999), Biostatistical Analysis (2nd ed.). NJ: Prentice-Hall,
%          Englewood Cliffs. p. 543-555. 
%

% examine the input of arguments to ensure that 4 frequencies are provided
% and determine the method of calculation of the p-value based on the 5th
% argument. If "mid" argument is passed into the function the mid-P-value
% approach is used.
if (nargin < 4)
  error('You need to input four arguments.');
  return;
endif
if ((nargin == 4)&&((sum([a b c d] < 0))||(sum(!(abs(mod([a b c d],1))==0)))))
  error('All frequencies should be positive integers.');
  return;
end;

% store the frequencies in the contigency table format and calculate
% the marginal totals R1, R2, C1, C2 as well as total observations n.

contigency_table = [a b;c d];
R1 = a+b;
R2 = c+d;
C1 = a+c;
C2 = b+d;
n = C1+C2;

% calculate the probability for the observed frequencies table
% using the logarithms of binomial coefficients instead of factorials
factor1 = sum([log(1:R1) log(1:R2) log(1:C1) log(1:C2) -log(1:n)]);
factor2 = sum([log(1:a) log(1:b) log(1:c) log(1:d)]);
P_observed = exp(factor1-factor2);


% find the lowest frequency of the contigency table and the related
% minimum marginal total to identify the pair of frequencies that will
% produce all possible combinations of the contigency table from the
% left (negative) to the right (positive) extreme tail. Consequently, all
% frequencies for all table combinations are calculated in a single
% (m+1 x 2 x 2) matrix, where m1 is the minimum marginal total.

lower_freq = min(min(contigency_table));
[row,column] = find(contigency_table == lower_freq);
row = row(1);
column = column(1);
if ((row == 1)&&(column == 1))
    m1 = min([R1 C1]);
    if (m1 == R1)
      ctable = [0 m1 C1 C2-m1];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[j -j -j j];
      endfor
    else
      ctable = [0 R1 m1 R2-m1];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[j -j -j j];
      endfor
    endif
  elseif ((row == 1)&&(column == 2))
    m1 = min([R1 C2]);
    if (m1 == R1)
      ctable = [m1 0 C1-m1 C2];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[-j j j -j];
      endfor
      else
      ctable = [R1 0 R2-m1 m1];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[-j j j -j];
      endfor
    endif
  elseif ((row == 2)&&(column == 1))
    m1 = min([R2 C1]);
    if (m1 == R2)
      ctable = [C1 C2-m1 0 m1];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[-j j j -j];
      endfor
    else
      ctable = [m1 R1-m1 0 R2];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[-j j j -j];
      endfor
    endif
  elseif ((row == 2)&&(column == 2))
    m1 = min([R2 C2]);
    if (m1 == R2)
      ctable = [C1-m1 C2 m1 0];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[j -j -j j];
      endfor
    else
      ctable = [R1-m1 m1 R2 0];
      ctable_array = zeros(m1+1,4);
      for i=1:m1+1
        j = i-1;
        ctable_array(i,:) = ctable.+[j -j -j j];
      endfor
    endif
endif

% Calculate the probabilities for all possible tables
for i=1:m1+1
  a_ct = ctable_array(i,1);
  b_ct = ctable_array(i,2);
  c_ct = ctable_array(i,3);
  d_ct = ctable_array(i,4);
  factor2 = sum([log(1:a_ct) log(1:b_ct) log(1:c_ct) log(1:d_ct)]);
  P_ctable_array(i) = exp(factor1-factor2);
endfor

% Sum all probabilities lower or equal to the probability of the
% observed frequencies
P2_tail = sum(P_ctable_array(find(P_ctable_array <= P_observed)));


% if output argument is given, then the output argument is
% the 2 tailed probability. When 3 output arguments are requested the P-values
% are parsed as follows: negative tail, positive tail and 2-tail probability. 
if (nargout == 1)
  varargout{1} = P2_tail;
elseif (nargout == 0)
  fprintf('\nFisher''s exact test for frequency cells:\n');
  fprintf('-------------------------------------------\n');
  fprintf('a = %i, b = %i\n', [a b]);
  fprintf('c = %i, d = %i\n', [c d]);
  fprintf('-------------------------------------------\n');
  fprintf('2-tail P-value: %10.10f\n',P2_tail);
end;
endfunction

