function mat_from_text(txtFile, resDirName, delim)

txt = fileread(txtFile);

mkdir(resDirName);
txt_sections = split(txt, delim);
for i=1:length(txt_sections)-1   
    files = split(txt_sections{i+1}, flip(delim));
    fileName = files{1};
    disp(fileName);
    fileData = files{2};
    
    id = fopen(['./' resDirName '/' fileName],'a');
    fprintf(id, '%s', fileData);
    fclose(id);
end

disp('Success!');
end