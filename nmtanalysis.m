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
function frequencies = nmtanalysis(filename,nlt)
% frequencies = nmtanalysis(filename,number_of_unilateral_traits)
%
% This function takes the .csv file containing the trait scores
% for a given population along with the number of unilateral traits
% and calculates their frequencies along with laterality, sexual
% dimorphism and inter trait correlation. The results are based on
% Fischer's exact test.
% All results are returned in a single data structure. Additionally,
% all results are also stored in a single spreadsheet file with the
% name of the population provided.
%
% The name of the population should correspond to the filename of the
% csv file without the extension, which is appended by the function.
% The number of the unilateral traits should be a non-negative integer.
%
% The trait scores should be 1 or 0 regarding presence or absence for
% each observation. For missing values NaN should be used, when no
% observation of a particular trait is available for an individual.
% The construction of the csv table should comply with the following
% rules:
%   The first column should list the names of the individuals
%   The second column should contain sex information of individuals,
%       0 for males, 1 for females, NaN for unknown
%   The consecutive columns should list unilateral traits if available,
%       followed by bilateral traits (left side should be listed first)
%   The first row should contain the names of the traits
%       (only letters and underscore should be used)
%
% The following example demostrates the required structure of a csv file
%
%   Individual,Sex,METO,APIC,INCA,PCTB,OMBL,OMBR,ASTL,ASTR
%   ind_name1,1,0,0,NaN,0,0,0,0,1
%   ind_name2,0,0,1,1,NaN,1,0,1,0
%   ind_name3,0,1,1,1,NaN,1,0,NaN,0
% 
%
% The function requires the io-package installed as well as the
% fisherextest.m and cellrange.m functions, which are necessary
% for Fisher's exact independence test and handling the indexing of
% the tables in speadsheets.

warning('off','Octave:divide-by-zero');
% Add .csv extension to the population name
filename = char(strcat(filename,".csv"));
% Read csv file and calculate its dimensions and load data into matrix
temp = csvread(filename);
row_len = length(temp(1,:));
col_len = length(temp(:,1));

% Remove first row and first column
trait_data = temp(2:end,2:end);

% Read string of traits from file
Traits = textread(filename, "%s", 1);
Trait_list = ostrsplit(Traits{1}, ",");
for i=3:row_len
  Trait{i-2,1} = Trait_list{:,i};
endfor

% Change the filename extension to .ods while retaining the initial
% name. The .ods format will be used for exporting results of the analysis.
filename = [strtrunc(filename,length(filename)-4) '.ods'];

% Measure population trait frequencies from data matrix

for i=2:row_len-1
  max_population = length(trait_data(:,i));
  nan_population = sum(isnan(trait_data(:,i)));
  trait_present = length(nonzeros(trait_data(:,i))) - nan_population;
  actual_population = length(trait_data(:,i)) - sum(isnan(trait_data(:,i)));
  percentage = (trait_present/actual_population);
  freq_table(i-1,:) = [trait_present actual_population percentage];
endfor

% Write frequency results in ods file on a sheet named 'Population Trait Frequencies'

header1 = {'Trait' 'Trait Presence' 'Actual Population' 'Percentage'};
celran = cellrange(1,1,length(header1(1,:)),length(header1(:,1)));
odswrite(filename, header1, 'Population Trait Frequencies',celran,'oct');
celran = cellrange(1,2,length(Trait(1,:)),length(Trait(:,1)));
odswrite(filename, Trait, 'Population Trait Frequencies',celran,'oct');
celran = cellrange(2,2,length(freq_table(1,:)),length(freq_table(:,1)));
odswrite(filename, freq_table, 'Population Trait Frequencies',celran,'oct');



% Perform trait laterality analysis. Calculate trait presence frequencies
% between sides and check for independence.
% Calculate the number and index of lateral traits.
lt = (length(trait_data(1,:)) - (nlt + 1))/2;
for i=1:lt
  xL_Tpr = 0;
  xL_Tab = 0;
  xR_Tpr = 0;
  xR_Tab = 0;
  not_available = 0;
  for j=1:col_len-1
    %Check for both left side and trait presence
    if (trait_data(j,(i*2)+nlt)==1)
      xL_Tpr = xL_Tpr +1;
    %Check for both left side and trait absence
    elseif (trait_data(j,(i*2)+nlt)==0)
      xL_Tab = xL_Tab +1;
    else
      not_available = not_available + 1;
    endif
    %Check for both right side and trait presence
    if (trait_data(j,(i*2)+1+nlt)==1)
      xR_Tpr = xR_Tpr +1;
    %Check for both right side and trait absence
    elseif (trait_data(j,(i*2)+1+nlt)==0)
      xR_Tab = xR_Tab +1;
    %Consider trait not available and exclude the individual
    else
      not_available = not_available + 1;
    endif
  endfor
  % Change the names of the bilateral traits by removing the 'L' at the end
  % that indicates side and by adding '_Lat' for clarification
  LTrait{i,1} = [strtrunc((Trait_list{:,(i*2)+1+nlt}),length(Trait_list{:,(i*2)+1+nlt})-1) '_Lat'];
  % constuct the contigency table for laterality independance and run
  % Fischer's exact test
  x = [xL_Tpr xL_Tab; xR_Tpr xR_Tab];
  L_P2_tail(i)=fisherextest(x(1), x(3), x(2), x(4));
  % save the individual counts for each case along with the 2-tailed p-value
  % of Fischer's exact test
  lat_freq_table(i,:) = [xL_Tpr xL_Tab xR_Tpr xR_Tab not_available L_P2_tail(i)];
endfor

% Write laterality results in ods file on a sheet named 'Laterality'

header2 = {'Trait' 'Left Side Trait Present' 'Left Side Trait Absent' ...
          'Right Side Trait Present' 'Right Side Trait Absent' 'Anavailable' ...
          'Fischer exact test'};
celran = cellrange(1,1,length(header2(1,:)),length(header2(:,1)));
odswrite(filename, header2, 'Laterality',celran,'oct');
celran = cellrange(1,2,length(LTrait(1,:)),length(LTrait(:,1)));
odswrite(filename, LTrait, 'Laterality',celran,'oct');
celran = cellrange(2,2,length(lat_freq_table(1,:)),length(lat_freq_table(:,1)));
odswrite(filename, lat_freq_table, 'Laterality',celran,'oct');


% Calculate the unilateral trait data matrix by merging the pairs 
% of lateral traits into single individual count traits.
for i=1:nlt+1
  uni_trait_data(:,i) = trait_data(:,i);
  UTrait{i,1} = Trait_list{:,i+2};
  U_Trait{1,i} = Trait_list{:,i+2};
endfor
for i=1:lt
  % Change the names of the bilateral traits by removing the '_Lat' at the end
  % that indicates side and by adding '_M' to declare merged traits
  UTrait{i+nlt,1} = [strtrunc(LTrait{i,1}, length(LTrait{i,1})-4) '_M']; 
  % apply the individual count scheme
  for j=1:col_len-1
    % Check if trait present on either side
    if (trait_data(j,(i*2)+nlt)==1) || (trait_data(j,(i*2)+1+nlt)==1)
      uni_trait_data(j,i+nlt+1) = 1;
    % Check if trait absent on both sides
    elseif (trait_data(j,(i*2)+nlt)==0) && (trait_data(j,(i*2)+nlt)==0)
      uni_trait_data(j,i+nlt+1) = 0;
    else
      uni_trait_data(j,i+nlt+1) = NaN;
    endif
  endfor
endfor


% Measure population trait frequencies from the merged unilateral trait matrix

for i=1:(length(uni_trait_data(1,:))-1)
  max_population = length(uni_trait_data(:,i+1));
  nan_population = sum(isnan(uni_trait_data(:,i+1)));
  uni_trait_present = length(nonzeros(uni_trait_data(:,i+1))) - sum(isnan(uni_trait_data(:,i+1)));
  uni_actual_population = length(uni_trait_data(:,i+1)) - sum(isnan(uni_trait_data(:,i+1)));
  uni_percentage = (uni_trait_present/uni_actual_population);
  uni_freq_table(i,:) = [uni_trait_present uni_actual_population uni_percentage];
endfor

% Write results in ods file on a sheet named 'Population Merged Trait Frequencies'
celran = cellrange(1,1,length(header1(1,:)),length(header1(:,1)));
odswrite(filename, header1, 'Population Merged Trait Frequencies',celran,'oct');
celran = cellrange(1,2,length(UTrait(1,:)),length(UTrait(:,1)));
odswrite(filename, UTrait, 'Population Merged Trait Frequencies',celran,'oct');
celran = cellrange(2,2,length(uni_freq_table(1,:)),length(uni_freq_table(:,1)));
odswrite(filename, uni_freq_table, 'Population Merged Trait Frequencies',celran,'oct');


% Perform trait sexual dimorphism analysis. Calculate trait presence frequencies
% between genders and check for independence.
% Calculate the number and index of unified traits.
ut = (length(uni_trait_data(1,:)))-1;
for i=1:ut
  xM_Tpr = 0;
  xM_Tab = 0;
  xF_Tpr = 0;
  xF_Tab = 0;
  not_available = 0;
  for j=1:col_len-1
    %Check for both male and trait presence
    if (uni_trait_data(j,1)==0) && (uni_trait_data(j,i+1)==1)
      xM_Tpr = xM_Tpr +1;
    %Check for both male and trait absence
    elseif (uni_trait_data(j,1)==0) && (uni_trait_data(j,i+1)==0)
      xM_Tab = xM_Tab +1;
    %Check for both female trait presence
    elseif (uni_trait_data(j,1)==1) && (uni_trait_data(j,i+1)==1)
      xF_Tpr = xF_Tpr +1;
    %Check for both female and trait absence
    elseif (uni_trait_data(j,1)==1) && (uni_trait_data(j,i+1)==0)
      xF_Tab = xF_Tab +1;
    %Consider trait not available and exclude the individual
    else
      not_available = not_available + 1;
    endif
  endfor
  if (i<=nlt)
    UTrait{i,1} = Trait_list{:,i+2};
    U_Trait{1,i} = Trait_list{:,i+2};
  elseif (i>nlt)
    % Change the names of the bilateral traits by removing the 'L' at the end
    % that indicates side and by adding '_M' to declare merged traits
    UTrait{i,1} = [strtrunc((Trait_list{:,((i-nlt)*2)+1+nlt}), ...
                  length(Trait_list{:,((i-nlt)*2)+1+nlt})-1) '_M'];
  endif
  % constuct the contigency table for laterality independance and run
  % Fischer's exact test
  x = [xM_Tpr xM_Tab; xF_Tpr xF_Tab];
  S_P2_tail(i)=fisherextest(x(1), x(3), x(2), x(4));
  % save the individual counts for each case along with the 2-tailed p-value
  % of Fischer's exact test
  sex_freq_table(i,:) = [xM_Tpr xM_Tab xF_Tpr xF_Tab not_available S_P2_tail(i)];
endfor
% Calculate per sex frequencies and number of actual population
for i=1:ut
  sex_percentage(i,:) = [((sex_freq_table(i,1))/((sex_freq_table(i,1)+...
          (sex_freq_table(i,2))))) ((sex_freq_table(i,1)+...
          (sex_freq_table(i,2)))) ((sex_freq_table(i,3))/((sex_freq_table(i,3)+...
          (sex_freq_table(i,4))))) ((sex_freq_table(i,3)+...
          (sex_freq_table(i,4))))];
          endfor
          
% Write sex dimorphism results in ods file on a sheet named 'Sex Dimorphism'
header3 = {'Trait' 'Mail Trait Present' 'Male Trait Absent' ...
          'Female Trait Present' 'Female Trait Absent' 'Anavailable' ...
          'Fischer exact test' 'Male presence' 'Male Population' ...
          'Female presence' 'Female Population'};
celran = cellrange(1,1,length(header3(1,:)),length(header3(:,1)));
odswrite(filename, header3, 'Sex Dimorphism',celran,'oct');
celran = cellrange(1,2,length(UTrait(1,:)),length(UTrait(:,1)));
odswrite(filename, UTrait, 'Sex Dimorphism',celran,'oct');
celran = cellrange(2,2,length(sex_freq_table(1,:)),length(sex_freq_table(:,1)));
odswrite(filename, sex_freq_table, 'Sex Dimorphism',celran,'oct');
celran = cellrange(8,2,length(sex_percentage(1,:)),length(sex_percentage(:,1)));
odswrite(filename, sex_percentage, 'Sex Dimorphism',celran,'oct');



% Perform inter trait independence test. Calculate trait presence frequencies
% in respect to other traits and check for independence.
% Calculate the number and index of unified traits.
for i=1:ut-1
  for k=i+1:ut
    xTipr_Tkpr = 0;
    xTiab_Tkpr = 0;
    xTipr_Tkab = 0;
    xTiab_Tkab = 0;
    not_available = 0;
    for j=1:col_len-1
      %Check for both i & k trait presence
      if (uni_trait_data(j,i+1)==1) && (uni_trait_data(j,k+1)==1)
        xTipr_Tkpr = xTipr_Tkpr +1;
      %Check for both i trait absence and k trait presence
      elseif (uni_trait_data(j,i+1)==0) && (uni_trait_data(j,k+1)==1)
        xTiab_Tkpr = xTiab_Tkpr +1;
      %Check for both i trait presence and k trait absence
      elseif (uni_trait_data(j,i+1)==1) && (uni_trait_data(j,k+1)==0)
        xTipr_Tkab = xTipr_Tkab +1;
      %Check for both i & k trait absence
      elseif (uni_trait_data(j,i+1)==0) && (uni_trait_data(j,k+1)==0)
        xTiab_Tkab = xTiab_Tkab +1;
      %Consider trait not available and exclude the individual
      else
        not_available = not_available + 1;
      endif
    endfor
    % constuct the contigency table for inter trait independance and run
    % Fischer's exact test
    x = [xTipr_Tkpr xTiab_Tkpr; xTipr_Tkab xTiab_Tkab];
    TI_P2_tail(i,k-1)=fisherextest(x(1), x(3), x(2), x(4));
    [TI_pval(i,k-1), TI_chisq(i,k-1), TI_df(i,k-1)] = chisquare_test_independence (x);
    % save the individual counts for each case along with the 2-tailed p-value
    % of Fischer's exact test
    trait_freq_table(i,k-1,:) = [xTipr_Tkpr xTiab_Tkpr xTipr_Tkab xTiab_Tkab ...
                                not_available TI_P2_tail(i,k-1)];
  endfor
endfor
% Write results in ods file on a sheet named 'Independence'
Tag = {'Fischer p'};
% The following saves the upper triangular matrix of Fischer's
% exact test p-values for all traits
celran = cellrange(2,1,length(UTrait(:,1))-1,length(UTrait(1,:)));
a = UTrait(2:end)';
odswrite(filename, a, 'Independence',celran,'oct');
celran = cellrange(1,2,length(UTrait(1,:)),length(UTrait(:,1))-1);
b = UTrait(1:end-1);
odswrite(filename, b, 'Independence',celran,'oct');
celran = cellrange(1,1,1,1);
odswrite(filename, Tag , 'Independence', celran,'oct');
celran = cellrange(2,2,length(TI_P2_tail(1,:)),length(TI_P2_tail(:,1)));
odswrite(filename, trait_freq_table(:,:,6), 'Independence',celran,'oct');



% Save frequency table for each analysis in a structure 
% for the function to return.
frequencies(1).fre = freq_table;
frequencies(1).lat = lat_freq_table;
frequencies(1).uni = uni_freq_table;
frequencies(1).sex = sex_freq_table;
frequencies(1).int = trait_freq_table;
endfunction
