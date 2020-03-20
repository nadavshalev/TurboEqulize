function text_from_mat(targetName, delim)

fn = uigetfile('MultiSelect', 'on');

id = fopen(targetName,'a');

for f = fn
    fprintf(id, delim);
    fprintf(id, [char(f) flip(delim)]);
    txt = fileread(char(f));
    fprintf(id, '%s', txt);
end

fclose(id);

disp('Success!');