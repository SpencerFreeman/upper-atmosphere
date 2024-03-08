function data = importfile(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  HOMEWORK3INPUTSLE020624 = IMPORTFILE(FILE) reads data from the first
%  worksheet in the Microsoft Excel spreadsheet file named FILE.
%  Returns the data as a table.
%
%  HOMEWORK3INPUTSLE020624 = IMPORTFILE(FILE, SHEET) reads from the
%  specified worksheet.
%
%  HOMEWORK3INPUTSLE020624 = IMPORTFILE(FILE, SHEET, DATALINES) reads
%  from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  data = importfile("C:\Users\spenc\Desktop\school shit 2\Spring_2024\upper_atmosphere\hw3\Homework_3_input_SLE_020624.xlsx", "Sheet1", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 20-Feb-2024 18:41:23

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 6);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["ICONDiskSolarZenithAngleDegrees", "BrightnessR", "MSISAltkm", "Ocm3", "ICONLimbAltituideKm", "BrightnessR1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Import the data
data = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    data = [data; tb]; %#ok<AGROW>
end

end