function Wing_VortexRing_BATCH
    [FileName,PathName] = uigetfile('*.xlsx');
    M_input = xlsread([PathName FileName],'ENTERA','B2:F11');
    M_output = zeros(size(M_input,1),6);
    for i=1:size(M_input,1)
        fprintf('-------------------------------\n');
        fprintf('CASE %i/%i\n',i,size(M_input,1));
        fprintf('-------------------------------\n');
        input_struct.b = M_input(i,1); % Wingspan
        input_struct.cr = M_input(i,2); % Chord at root
        input_struct.ct = M_input(i,3); % Chord at tip
        input_struct.s = M_input(i,4); % Sweep angle
        input_struct.CL_target = M_input(i,5); % Required CL   
        output_struct = Wing_VortexRing(input_struct);
        M_output(i,:) = [output_struct.Ar  output_struct.CLCDi output_struct.CL output_struct.CDi output_struct.alpha output_struct.epsilon];
    end
    xlswrite([PathName FileName],M_output,'ENTERA','G2');
    keyboard;
end