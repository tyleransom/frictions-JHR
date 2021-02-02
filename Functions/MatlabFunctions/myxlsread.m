
% Wrapper for the xlsread function so that it works on a non-pc platform
% using the open document format. Do not specify the extension on the input
% file, as .xls or .ods is automatically appended depending on the
% platform. Other arugments and outputs are the same as for the xlsread
% command.
%
% xxx - currently all input arguments except for mode must be specified 
%   (one could add nargin checking for optional inputs).
%
% Requires import of the odfdom library (java). Modify the static
% javaclasspath so that it includes odfdom.jar
% http://odftoolkit.org/projects/odfdom/pages/Home
%
% odfdom is dependent on xerces-j libraries as well for parsing xml, but 
% as of matlab 2008b these are already included automatically by matlab.
%
% Introduced mode input which if set to something other than 'basic' or
% empty string will force compatibility with column and row headers being
% ignored (same as PC implementation for xls files) for num output.
% Default is compatibility mode.

function [num str] = myxlsread(infile,wksht,rng,mode)

% call function as normal with xls on a pc using excel COM interface.
if ispc
  infile = [infile '.xls'];
  [num str] = xlsread(infile,wksht,rng);
  return;
end

% on a mac, use the odfdom interface to read the ods file
num = []; str = {}; infile = [infile '.ods'];

% mode switch default is for compatilibity with PC xlsread
if nargin < 4, mode = ''; end
if isempty(mode), mode = ''; end
switch mode
  case {'basic' ''}
    ignore_hdrs = true;
  otherwise
    ignore_hdrs = false;
end

import org.odftoolkit.odfdom.doc.*;
import org.odftoolkit.odfdom.doc.table.*;
import org.odftoolkit.odfdom.doc.office.*;

try
  odsDoc = OdfSpreadsheetDocument.loadDocument(infile);
  odsSpreadSheet = odsDoc.getContentRoot();
  odsTables = odsSpreadSheet.getChildNodes();
catch ME
  error(['Can not read ' infile ' ' ME.identifier]);
end

% iterate over the tables until we find the worksheet that we want.
len = odsTables.getLength(); found = 0;
for i=1:len
  child = odsTables.item(i-1);
  if strcmpi(child.getLocalName(),'table')
    odsTable = OdfTable.getInstance(child);
    if strcmp(odsTable.getTableName(),wksht), found = 1; break; end
  end
end
if ~found, return; end

% get the range of cells specified
i = strfind(rng,':'); brng = rng(1:i(1)-1); erng = rng(i(1)+1:end);
odsCells = odsTable.getCellRangeByPosition(brng,erng);

% iterate over the range and assign string and numeric values to outputs.
% only assign "float" and "string" types.
m = odsCells.getRowNumber(); n = odsCells.getColumnNumber(); 
str = cell(m,n); num = nan(m,n); firstnum = [inf inf]; lastcell = [0 0];
for r=1:m
  for c=1:n
    odsCell = odsCells.getCellByPosition(c-1,r-1);
    type = char(odsCell.getValueType());
    if isempty(type)
      str{r,c} = ''; num(r,c) = NaN;
    else
      switch lower(type)
        case 'float'
          num(r,c) = double(odsCell.getDoubleValue());
          str{r,c} = '';
          if r < firstnum(1), firstnum(1) = r; end
          if c < firstnum(2), firstnum(2) = c; end
        case 'string'
          str{r,c} = char(odsCell.getStringValue());
          num(r,c) = NaN;
        otherwise
          str{r,c} = ''; num(r,c) = NaN;
      end
      if r > lastcell(1), lastcell(1) = r; end
      if c > lastcell(2), lastcell(2) = c; end
    end % if cell is empty
  end % for each column
end % for each row

% for strict compatibility with PC xlsread (default).
% use anything other than 'basic' or '' for mode to return full num matrix
%   which has NaNs in all string locations for range being read.
if ignore_hdrs
  if all(isfinite(firstnum))
    num = num(firstnum(1):end,firstnum(2):end); 
  else
    num = [];
  end
end

% trim the read to include only non-empty cells
num = num(1:lastcell(1),1:lastcell(2));
str = str(1:lastcell(1),1:lastcell(2));

odsDoc.close(); % possible leaking memory without close?
