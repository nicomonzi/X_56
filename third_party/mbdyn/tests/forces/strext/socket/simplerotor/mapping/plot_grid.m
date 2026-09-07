function [POut, VOut] = plot_grid(gridf, movf, varargin)

locframe = 1;
new_plot = 0;

if nargin > 2
	PStep = varargin{1};
else
	PStep = 1;
end

fid = fopen(gridf, 'r');
if (fid < 0)
	disp(['Errors in the opening of the file ', gridf]);
	return;
end
l = 0;
NodeN = 0;
while ~(feof(fid))
	tline = fgetl(fid);
        if ~strcmp(tline(1), '#') &&  ~isempty(tline)
		l = l + 1;
		if (l == 1)
			NodeN = sscanf(tline, '%d', 1);
			Nodes = zeros(NodeN, 3); 	
		else if (l < NodeN + 2)
				Nodes(l-1, :) = sscanf(tline, '%f', 3);
			else if (l == NodeN + 2)
					PanN = sscanf(tline, '%d', 1);
					Conn = zeros(PanN, 4);
				else
					Conn(l - (NodeN + 2) , :) = sscanf(tline, '%d', 4);
				end
			end
		end
	end
end
fclose(fid);
%NNodes = zeros(3, NodeN);
%NVel = zeros(3, NodeN);
fid = fopen(movf, 'r');
if (fid < 0)
	disp(['Errors in the opening of the file ', movf]);
	return;
end
i = 1;			
while ~(feof(fid))	
	S = fscanf(fid, 'STEP %d  ITERATION %d \n')
	Step = S(1);
	Iteration = S(2);
	tline = fgetl(fid);
	tline = fgetl(fid);
	P0 = sscanf(tline, '%f', 3)
	tline = fgetl(fid);
	R  = sscanf(tline(2:end), '%f %f %f', [3, 3])
	tline = fgetl(fid);
	V0 = sscanf(tline, '%f', 3)
	tline = fgetl(fid);
	W0 = sscanf(tline, '%f', 3)
	tline = fgetl(fid);
	Nd = fscanf(fid, 'POS %d\n');
	if (Nd ~= NodeN)
		error(['Wrong node number for POS list at Step ',num2str(Step)]);
	end
	NNodes = (fscanf(fid, '%f', [3, NodeN]))';
	tline = fgetl(fid);
	Nd = fscanf(fid, 'VEL %d\n');
	if (Nd ~= NodeN)
		error(['Wrong node number for VEL list at Step ',num2str(Step)]);
	end 
	NVel = (fscanf(fid, '%f', [3, NodeN]))';
	tline = fgetl(fid);
    if ~locframe
        NVel = ones(NodeN, 1) * V0' + (R * NVel')';
        NNodes = ones(NodeN, 1) * P0' + (R * NNodes')'; 
    end
    if not(mod(Step, PStep))
		POut(:,:,i) = NNodes;
		VOut(:,:,i) = NVel;
		i = i + 1;
        h = patch('Faces', Conn, 'Vertices', NNodes, ...
            'FaceColor','y','EdgeColor','r','LineWidth',.8);
        axis equal;
        view(3);
		pause;
		if new_plot
			close(gcf);
		end
	end
end
fclose(fid);
return;
