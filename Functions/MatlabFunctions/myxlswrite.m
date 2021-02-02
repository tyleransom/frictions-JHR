
% Wrapper for the xlswrite function so that it works on a non-pc platform
% using the open document format. Do not specify the extension on the input
% file, as .xls or .ods is automatically appended depending on the
% platform. Other arugments and outputs are the same as for the xlswrite
% command.
%
% xxx - currently all input arguments must be specified 
%   (one could add nargin checking for optional inputs).
%
% Requires import of the odfdom library (java). Modify the static
% javaclasspath so that it includes odfdom.jar
% http://odftoolkit.org/projects/odfdom/pages/Home
%
% odfdom is dependent on xerces-j libraries as well for parsing xml, but 
% as of matlab 2008b these are already included automatically by matlab.
%
function [status msg] = myxlswrite(outfile,A,wksht,rng)

% call function as normal with xls on a pc using excel COM interface.
if ispc
  outfile = [outfile '.xls'];
  [status msg] = myxlswrite(outfile,A,wksht,rng);
  return;
end

% on a mac, use the odfdom interface to write the ods file
status = true; msg = ''; outfile = [outfile '.ods'];

import org.odftoolkit.odfdom.doc.*;
import org.odftoolkit.odfdom.doc.table.*;
import org.odftoolkit.odfdom.doc.office.*;

if exist(outfile,'file')
  try
    odsDoc = OdfSpreadsheetDocument.loadDocument(outfile);
  catch ME
    status = false; msg = ME.identifier; return;
  end
else
  % if the file doesn't exist, create a new spreadsheet document and add
  try
    odsDoc = OdfSpreadsheetDocument.newSpreadsheetDocument();
  catch ME
    status = false; msg = ME.identifier; return;
  end
end

try
  odsSpreadSheet = odsDoc.getContentRoot();
  odsTables = odsSpreadSheet.getChildNodes();
catch ME
  status = false; msg = ME.identifier; return;
end

% iterate over the tables until we find the worksheet that we want
len = odsTables.getLength(); found = 0;
for i=1:len
  child = odsTables.item(i-1);
  if strcmpi(child.getLocalName(),'table')
    odsTable = OdfTable.getInstance(child);
    if strcmp(odsTable.getTableName(),wksht), found = 1; break; end
  end
end

% add the worksheet if not found
if ~found
  try
    odsTable = OdfTable.newTable(odsDoc);
    odsTable.setTableName(wksht);
  catch ME
    status = false; msg = ME.identifier; return;
  end
end

% check number of rows and columns of data to be written
[wr wc] = size(A); 

% get the range of cells specified
i = strfind(rng,':'); 
if isempty(i)
  % only the beginning of the range was provided
  brng = rng;
  [br bc] = cellID_to_row_col(brng);
  
  % figure out the end range based on how much to write.
  er = br+wr-1; ec = bc+wc-1;
  erng = row_col_to_cellID(er,ec);
else
  brng = rng(1:i(1)-1); erng = rng(i(1)+1:end);
end

% get the new range which should include all cells to be written.
% seems that odfdom automatically adds if they are not present so that
%   appending directly in this function is not necessary.
odsCells = odsTable.getCellRangeByPosition(brng,erng);
m = odsCells.getRowNumber(); n = odsCells.getColumnNumber(); 

% error if not enough data provided for specified range.
if wr < m || wc < n
  status = 0; msg = 'Not enough data provided for range specified'; return;
end

% iterate over the range and write data to cells.
for r=1:m
  for c=1:n
    odsCell = odsCells.getCellByPosition(c-1,r-1);
    
    if iscell(A)
      towrite = A{r,c};
    else
      towrite = A(r,c);
    end

    if ischar(towrite)
      odsCell.setStringValue(java.lang.String(towrite));
    elseif ~isempty(towrite)
      odsCell.setDoubleValue(java.lang.Double(towrite));
    end
  end
end

% save and close the new spreadsheet
odsDoc.save(outfile); odsDoc.close();

% couldn't find these functions in the odfdom, maybe missed?
function [r c] = cellID_to_row_col(id)
c = 0;
for i=1:length(id)
  n = str2double(id(i));
  if isfinite(n), break; end
  c = 26*c + upper(id(i)) - 'A' + 1;
end
r = str2double(id(i:end));

function id = row_col_to_cellID(r,c)
id = '';
while c > 0
  c = c-1; % weird spreadsheet base, no zero place
  id = [char('A' + mod(c,26)) id];
  c = fix(c/26);
end
id = [id num2str(r)];
