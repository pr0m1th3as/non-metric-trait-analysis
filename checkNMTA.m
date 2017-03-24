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
function freq_table = checkNMTA(a,population,p_lat,p_sex,p_int)
% freq_table = checkNMTA2(a,population,p_lat,p_sex,p_int)
%
% This function takes the results of the nmtanalysis.m function contained
% in structure (a) and checks for lateraliry, sexual dimorphism and trait
% intercorrelation based on the significance levels given in arguments
% (p_lat), (p_sex), (p_int) respectively. The second argument (population)
% is a cell array of strings containing the names of the populations
% being analyzed in the (a) structure. Concequently, the names should be
% in the respective order.
%
% After checking the results of the former analysis contained in the
% structure variable (a), the funtion will calculate a 3-dimensional matrix
% with the appropriate data format for Smith's MMD calculation utilizing the
% smithsMMD.m function (Bertsatos & Chovalopoulou, 2016). The first dimension
% corresponds to populations, the second corresponds to traits that meet
% specific criteria and the third to frequency of each trait containing
% number of occurences and number of observations accordingly.
%
% Regarding laterality, the side with the highest expression is included in
% the frequency table, whereas the traits that do not exhibit laterality are
% included based on the individual count method (Turner et al., 1991). The
% traits that exhibit sexual dimorphism are excluded from the frequency table
% automatically, whereas the function outputs on the console the cases of 
% intercorrelated traits found across the population samples and the user is
% asked to select from a list box which traits should be excluded irregardless
% of which pairs were found intercorrelated. If the sexually dimorphic traits
% should be included in the frequency table, the value of (p_sex) should be 
% set to 0.
% 
% Before returning the final frequency matrix, the function checks that each
% trait exhibits statistical significant difference in at least one pair of
% populations across the populations included and if a trait fails to exhibit
% such difference in any population pair it is excluded, since it is considered
% redundant to MMD calculations (Harris & Sjøvold, 2004).
%
%
%
% References:
%
% Bertsatos, A., & Chovalopoulou, M.-E. (2016). Technical note: 
%          A GNU octave function for smith’s mean measure of divergence.
%          Bioarchaeology of the Near East, 10, 69-73.
% Harris E.F., Sjøvold T. (2004), Calculation of Smith’s Mean Measure of
%          Divergence for intergroup comparisons using nonmetric data, Dental
%          Anthropology 17(3):83-93.
% Turner, C.G., Nichol, C.R., & Scott, G.R. (1991). Scoring procedures for key
%          morphological traits of the permanent dentition: The Arizona State
%          University dental anthropology system. In: Kelley MA, Larsen CS.,
%          Advances in dental anthropology. New York: Wiley & Sons. p. 13-31.
%

% check for each population which traits are observable
for i=1:length(a)
q=logical(a(i).uni(:,2)!=0);
traits_observable(i,:) = q';
endfor
% build a vector with traits present in all populations
present_traits = logical(sum(traits_observable)==length(a));
% display the traits observable across all populations
printf('\n%i traits out of the %i examined occur in all populations.\n', ...
      sum(present_traits==1),length(a(1).uni(:,2)))
printf('Traits present are:\n%s\n',disp(find(present_traits==1)));


% check for each population which traits do exhibit laterality for
% a given threshold of p-value defined in p_lat
printf('Laterality examined at significance level of %.2f\n',p_lat);
for i=1:length(a)
  q=logical(a(i).lat(:,6)<p_lat);
  traits_with_laterality(i,:) = q';
endfor
% calculate the number of unilateral traits present in the sample
ult = length(present_traits)-length(traits_with_laterality(1,:));
% display which traits exhibit laterality for each population
for i=1:length(a)
  if ((sum(traits_with_laterality(i,:)==1))==0)
    printf('  No traits exhibit laterality in %s', ...
      disp(char(population(i))));
  elseif ((sum(traits_with_laterality(i,:)==1))==1)
    printf('  Trait%s exhibits laterality in %s', ...
      strtrunc(disp(find(traits_with_laterality(i,:)==1)+ult), ...
      length(disp(find(traits_with_laterality(i,:)==1)))-1), ...
      disp(char(population(i))));
  elseif ((sum(traits_with_laterality(i,:)==1))>1)
    printf('  Traits%s exhibit laterality in %s', ...
      strtrunc(disp(find(traits_with_laterality(i,:)==1)+ult), ...
      length(disp(find(traits_with_laterality(i,:)==1)))-1), ...
      disp(char(population(i))));
  endif
endfor


% check for each population which traits do not exhibit sexual
% dimorphism for a given threshold of p-value defined in p_sex
printf('\nSexual dimorphism examined at significance level of %.2f\n',p_sex);
for i=1:length(a)
q=logical(a(i).sex(:,6)>p_sex);
traits_with_sex_dimorphism(i,:) = q';
% display which traits exhibit sexual dimorphism for each population
  if ((sum(traits_with_sex_dimorphism(i,:)==0))==0)
    printf('  No traits exhibit sexual dimorphism in %s', ...
      disp(char(population(i))));
  elseif ((sum(traits_with_sex_dimorphism(i,:)==0))==1)
    printf('  Trait%s exhibits sexual dimorphism in %s', ...
      strtrunc(disp(find(traits_with_sex_dimorphism(i,:)==0)), ...
      length(disp(find(traits_with_sex_dimorphism(i,:)==0)))-1), ...
      disp(char(population(i))));
  elseif ((sum(traits_with_sex_dimorphism(i,:)==0))>1)
    printf('  Traits%s exhibit sexual dimorphism in %s', ...
      strtrunc(disp(find(traits_with_sex_dimorphism(i,:)==0)), ...
      length(disp(find(traits_with_sex_dimorphism(i,:)==0)))-1), ...
      disp(char(population(i))));
  endif
endfor
% build a vector with traits that do not exhibit sex dimorphism in all populations
pooled_sex_traits = logical(sum(traits_with_sex_dimorphism)==length(a));


% check for each population which traits are intercorrelated

int_cor_ind = 1;
traits_intercorrelated(int_cor_ind,:) = [0 0 0];
for i=1:length(a)
q=logical(a(i).int(:,:,6)>p_int);
[row,column]=find(xor(!q,tril(ones(length(q),length(q)),-1)));
  for len=1:length(row)
    traits_intercorrelated(int_cor_ind,:) = [i row(len) column(len)+1]';
    int_cor_ind = int_cor_ind + 1;
  endfor
endfor

% Contruct a 3-dimensional matrix for Smith's MMD calculation.
% The first dimension corresponds to populations, the second corresponds
% to traits that meet specific criteria and the third to frequency of each
% trait containing number of occurences and number of observations accordingly
% initialize loop for every population examined
for i=1:length(a)
  % initialize index for trait to be included in the final frequency table
  trait_index = 1;
  % initialize loop for each trait to be examined
  for j=1:length(a(1).uni(:,1))
    % make decision tree based on the criteria which traits are
    % to be included in the final frequency table
    if (j<=ult) && (present_traits(j)==1) && (pooled_sex_traits(j)==1)
      freq_table(i,trait_index,:) = a(i).uni(j,[1 2]);
      % save number of examined trait and increment index
      index(i,trait_index,:) = j;
      trait_index = trait_index + 1;
    elseif (j>ult) && (present_traits(j)==1) && ...
           (pooled_sex_traits(j)==1) && (traits_with_laterality(i,j-ult)==0)
      freq_table(i,trait_index,:) = a(i).uni(j,[1 2]);
      % save number of examined trait and increment index
      index(i,trait_index,:) = j;
      trait_index = trait_index + 1;
    elseif (j>ult) && (present_traits(j)==1) && ...
           (pooled_sex_traits(j)==1) && (traits_with_laterality(i,j-ult)==1)
        % check if left side trait occurs more than right side and use left
        if (a(i).fre(ult-1+((j-ult)*2),1)>a(i).fre(ult+((j-ult)*2),1))
          % use left side frequency
          freq_table(i,trait_index,:) = a(i).fre(ult-1+((j-ult)*2),[1 2]);
          % save number of examined trait and increment index
          index(i,trait_index,:) = j;
          trait_index = trait_index + 1;
        % check if left side trait occurs less of equal than right side and use right
        elseif (a(i).fre(ult-1+((j-ult)*2),1)<=a(i).fre(ult+((j-ult)*2),1))
          % use right side frequency
          freq_table(i,trait_index,:) = a(i).fre(ult+((j-ult)*2),[1 2]);
          % save number of examined trait and increment index
          index(i,trait_index,:) = j;
          trait_index = trait_index + 1;
        endif
    endif
  endfor
endfor

% since the traits inlcuded in the frequency table are the same
% across all populations, reduce the index matrix to a single row vector
T_index = index(1,:);
% display the traits included in the frequency table along with the number
% of occurences and observations for each trait
printf('\nThe following traits meet the specific criteria for laterality, \n');
printf('sex dimorphism and consistent occurence across all populations.\n');
printf('%s\n',disp(T_index));
printf('The following table contains number of occurences for each\n');
printf('trait across all populations.\n');
printf('%s\n',disp(freq_table(:,:,1)));
printf('The following table contains number of observations for each\n');
printf('trait across all populations.\n');
printf('%s\n',disp(freq_table(:,:,2)));

% display the significance level at which trait intercorrelation was examined
printf('\nInter trait correlation examined at significance level of %.2f\n',p_int);
% display the traits found to be significantly intercorrelated for each
% population and ask the user if any trait should be eliminated
printf('The following trait pairs were found to be significantly\n');
printf('intercorrelated on specified populations.\n');
printf('%s\n',disp(traits_intercorrelated));

% prepare a cell array of strings of the traits included in the
% frequency table so far
for i=1:length(T_index)
  traitlist(i) = {num2str(T_index(i))};
endfor
% prompt the user to select the traits, which are to be excluded
printf('Please, select from the list which traits should be excluded from.\n');
printf('the frequency table due to intercorrelation. You can choose one or\n');
printf('multiple traits at once. If you wish to include all traits deselect\n');
printf('all traits and press ok.\n');
% display a list dialog for the user to select the traits that should be
% excluded from the frequency table
T_exclude = listdlg("ListString",traitlist,"SelectionMode", "Multiple");

% display the traits excluded due to inter correlation after user's selection
printf('\nTraits excluded due to inter trait correlation.\n');
printf('%s\n',disp(T_index(T_exclude)));

% remove user selected traits from the frequency table and the trait list
freq_table(:,T_exclude,:) = [];
T_index(T_exclude) = [];
% and display the result
printf('The following traits are left in the frequency table\n');
printf('%s\n',disp(T_index));
printf('The following table contains number of occurences for each,\n');
printf('trait across all populations.\n');
printf('%s\n',disp(freq_table(:,:,1)));
printf('The following table contains number of observations for each,\n');
printf('trait across all populations.\n');
printf('%s\n',disp(freq_table(:,:,2)));


% Last test before feeding the frequency table into the Smith's MMD
% calculations examines that each trait exhibits statistical significant
% difference in at least one pair of populations.
[p,t,f]=size(freq_table);
for i=1:t
  T_sig_dif(i) = 0;
  for j=1:p-1
    for k=j+1:p
      pr_P1 = freq_table(j,i,1);
      npr_P1 = freq_table(j,i,2)-freq_table(j,i,1);
      pr_P2 = freq_table(k,i,1);
      npr_P2 = freq_table(k,i,2)-freq_table(k,i,1);
      P2_tail=fisherextest(pr_P1,npr_P1,pr_P2,npr_P2);
      if (P2_tail<=0.05)&&(T_sig_dif(i)==0)
        T_sig_dif(i) = 1;
      elseif (P2_tail<=0.05)&&(T_sig_dif(i)==1)
        T_sig_dif(i) = 1;
      endif
    endfor
  endfor
endfor
% display the traits that did not exhibit any statistical significant
% difference between any population pair.
printf('Traits excluded due to statistically non significant\n');
printf('differences among populations.\n');
printf('%s\n',disp(T_index(find(T_sig_dif==0))));

% remove the traits with non statistical significant difference between
% any population pair from the frequency table and the relevant trait list
T_index(find(T_sig_dif==0)) = [];
freq_table(:,find(T_sig_dif==0),:) = [];

printf('The following traits are included in the final frequency table\n');
printf('%s\n',disp(T_index));
printf('The following table contains number of occurences for each,\n');
printf('trait across all populations.\n');
printf('%s\n',disp(freq_table(:,:,1)));
printf('The following table contains number of observations for each,\n');
printf('trait across all populations.\n');
printf('%s\n',disp(freq_table(:,:,2)));
