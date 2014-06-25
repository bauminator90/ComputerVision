function coordinates = readplyfile(filename)
    fileid=fopen(filename,'rt');
    pointcount='';
    indheader=0;
    count=1;

    data=fgetl(fileid);
    while ischar(data)
        %Read-in vertice-count
        if size(data,2)>= 14
            if strcmp(data(1,1:14),'element vertex')
                vertexcount=data(1,15:end);
                vertexcount=str2num(vertexcount);
            end
        end
        
        %stop reading when reading 'end_header'.
        if strcmp(data,'end_header ')|strcmp(data,'end_header')
            indheader=count;
            break;
        end
        data = fgetl(fileid);
        count=count+1;
    end
    fclose(fileid);

    fileid=fopen(filename,'rt');
    data=fgetl(fileid);
    count=1;
    coordinates=[];
    
    while ischar(data)
        if (count >= indheader+1) & (count <= indheader + vertexcount)
            format long;
            coordinates=[coordinates;str2num(data)];
        end
        data=fgetl(fileid);
        count=count+1;
    end
   
end