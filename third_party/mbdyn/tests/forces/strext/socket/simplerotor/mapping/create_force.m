function   create_force(forcef, NNodes, Steps, Scale, Ref)

fid = fopen(forcef, 'w');
if (fid < 0)
	disp(['Errors in the opening of the file ', forcef]);
	return;
end

for i = 1 : Steps
    fprintf(fid,'Step %d\n', i - 1);
    if Ref
        fprintf(fid, 'REF %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n', zeros(1,6)); 
    end
    for j = 1 : NNodes
        fprintf(fid, '%16.8e %16.8e %16.8e\n', Scale*randn(1,3));
    end
end
fclose(fid);