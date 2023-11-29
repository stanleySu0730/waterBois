folder_path = 'C:\Users\stanl\MatLab\waterBois';
file_name = 'cases_wastewater_vaccine.csv';
file_path = fullfile (folder_path,file_name);
data = readtable(file_path);
file_name3 = 'id_zipcode_income_density.csv';
file_path3 = fullfile(folder_path,file_name3);
info = readtable(file_path3); 
incomeArray = info.MedianIncome;
densityArray = info(:,4);
ageArray = info(:,5);

if any(data(:, 1) == 1)
    filteredData = data(data(:, 1) == 1, :);
else
    disp('No entries found with the first column equal to 1.');
end
