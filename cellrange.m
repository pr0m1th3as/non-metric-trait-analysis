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
function cellrange = cellrange(startcol,startrow,len_row,len_col)
% function cellrange = cellrange(startcol,startrow,len_row,len_col)
%
% This function takes the integer coordinates of the first cell
% and the respective row length and column length of a table
% and returns the coresponding cell range of a spreadsheet as a string.
% The first argument defines the column and the second the row.
%
% All four input arguments are required as shown in following example.
%
% cellrange(2,3,5,6)
% ans = B3:F8
%
% The first cell of the table is expected in the column range A-Z.
% The length of a row cannot exceed 1024 allowing the last column
% to be AMJ in compliance with Libreoffice Calc limitations.

number1 = num2str(startrow);
number2 = num2str(startrow+len_col-1);
end_row = startcol+len_row-1;
if (startcol > 26)
  cellrange = 0;
  error ("cellref: expecting first column to start within A-Z range");
  elseif ((startcol+len_row-1) > 1024)
    cellrange = 0;
    error ("cellref: the maximum number of column cannot exceed 1024");
    elseif (end_row <= 26)
    cellrange = strcat(char(64+startcol), number1, ":", char(64+end_row), number2);
    elseif (end_row > 26)&&(end_row <= 702)
      if (rem(end_row,26) != 0)
      first_letter = floor(end_row/26);
      second_letter = rem(end_row,26);
      elseif (rem(end_row,26) == 0)
      first_letter = floor(end_row/26)-1;
      second_letter = 26;
      endif
      cellrange = strcat(char(64+startcol), number1, ":", char(64+first_letter), ...
                  char(64+second_letter), number2);
    elseif (end_row > 702)&&(end_row <= 1024)
      remain = end_row - 702;
      if (rem(remain,26) != 0)
      second_letter = floor(remain/26)+1;
      third_letter = rem(remain,26);
      elseif (rem(remain,26) == 0)
      second_letter = floor(remain/26);
      third_letter = 26;
      endif
      cellrange = strcat(char(64+startcol), number1, ":", "A", ...
                  char(64+second_letter), char(64+third_letter), number2);
endif
endfunction




