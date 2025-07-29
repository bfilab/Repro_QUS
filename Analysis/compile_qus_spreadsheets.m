%Compile spreadsheets

main_dir = 'D:\Andrew\Prolapse Model';
cd(main_dir);
directories = {dir('M*')}; %Find all RTG folders
new_table = table; %Start a new table to compile data

for d = 1:length(directories)
    samples = directories{d};
    for i = 1:length(samples)
        this_id = samples(i).name; %Find animal ID
        cd(this_id); %Go to folder and save current directory
        csv_files = dir('2025-07-08*.csv'); %Find csv files with QUS data

        for f = 1:length(csv_files)
            this_table = readtable(csv_files(f).name); %Extract table
            %this_table.ID(:) = {this_id}; %Ensure the correct ID
            new_table = [new_table; this_table]; %Add to compiled table
        end
        cd(main_dir);
    end
end

save_str = strcat(char(datetime(datetime,'Format','yyyy-MM-dd')),...
    '-first-frame.csv');
writetable(new_table,save_str); %Save compiled table