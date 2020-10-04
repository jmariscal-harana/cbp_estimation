%%	Based on Pete's algorithm, based on Jordi's, based on Marie's


function WriteInFile(PATHS,PARAM,MESH,BCS,ICS,HIS)

% ========================
% Trying to only write either visco-elastic or elastic files.
% START HERE
% ========================


%% Open IN file
fileIn = fopen([PATHS.Aortic1D_files, PATHS.PatientName, '.in'],'w');


%%	1. PARAMETER LIST
varName		= fieldnames(PARAM);
nbP			= numel(varName);

%	HEADER
fprintf(fileIn,'%d\t%s\n',nbP,'Number of parameters');

% Print each parameter value to .in file
for k=1:nbP
	eval(['varVal(k) = PARAM.',varName{k},'{1};']);
	eval(['varType{k} = PARAM.',varName{k},'{2};']);
    fprintf(fileIn,[varType{k},'\t%s\n'],varVal(k),varName{k});
end


%%	2. MESH DEFINITION
%	HEADER
fprintf(fileIn,'%s %d \n','Mesh -- Expansion order --  Quadrature order -- Ndomains =',MESH.nbSeg);

%	Specify the format of mechanical properties: 'Beta', 'Eh' or 'Empirical', as explained in the Reference Manual	
switch MESH.MechProps
	case 'Beta'
		% Write values of Area, Eh and (and Gamma if visco-elastic) for each arterial segment (each domain)
		for e=1:MESH.nbSeg
			
			% Print header line stating: Eh, Area (and Gamma if visco-elastic)
			fprintf(fileIn,'%d %s %d %s\n',MESH.nbElemBySeg(e),'    nel domain ',e, 'Beta area');  %If prescribed with empirical PWV and Area

			% Make string for radius
			Radius = ['(',num2str(MESH.Rad_coef(e,1)),'*x+',num2str(MESH.Rad_coef(e,2)),')'];
			
			% Print individual values for each element in this segment
			x_temp = MESH.x_elem{e};
			% cycle through elements
			for j=1:MESH.nbElemBySeg(e)
				
				% identify distances (x values)
				l_start = x_temp(j);
				l_end   = x_temp(j+1);
				
				% print mesh settings
				fprintf(fileIn,'%3.4f %3.4f %d %d %s\n', l_start, l_end, MESH.p(e), MESH.q(e), '# x_prox x_dist p q');
				
				% print beta - IMPORTANT: print beta before area
				fprintf(fileIn,'%s %f %s %s %s %s %s\n','Beta = ', MESH.beta_coef(e) ,'/sqrt(PI*',Radius,'*', Radius,')');
				
				% print area
				fprintf(fileIn,'%s %s %s %s\n','Area = PI*', Radius,'*',Radius);
				
				% print gamma if this segment is visco-elastic
				if (MESH.L(e)>MESH.elmt_shortL_vw) && MESH.Viscoeslastic == 1
					error('Does this even make sense?')
					fprintf(fileIn,'%s %s\n','Gamma = ', Gamma); %Visco-elastic coefficient
				end
			end
		end
		
	case 'Eh'
		% Write values of Area, Eh and (and Gamma if visco-elastic) for each arterial segment (each domain)
		for e=1:MESH.nbSeg
			
			% Make string for radius
			Radius = ['(',num2str(MESH.Rad_coef(e,1)),'*x+',num2str(MESH.Rad_coef(e,2)),')'];
			
			% Make string for Eh
			Eh = [num2str(0.1),'*(',num2str(MESH.Eh_k1),'*exp(',num2str(MESH.Eh_k2),'*100*',Radius,') +',num2str(MESH.Eh_k3),')*',Radius];
			
			% Make string for Gamma
			Gamma = ['4/PI*(',num2str(MESH.gamma_b1),'/(2*100*',Radius,')+',num2str(MESH.gamma_b0),')/1000/(2*',Radius,')^2'];

			% Pete's version: ghigo_gamma = (2/3)*sqrt(pi)*0.017.*Eh./(pi*(r.^2));
			
% 			Gamma = ['(2/3)*sqrt(PI)*0.017*' Eh '/(PI*(' Radius '*' Radius '))'];
			
% 			Beta = [num2str(0.1*4/3/sqrt(pi)),'*(',num2str(MESH.Eh_k1),'*exp(',num2str(MESH.Eh_k2),'*0.01*',Radius,') +',num2str(MESH.Eh_k3),')/',Radius];
			
			% Print individual values for each element in this segment
			x_temp = MESH.x_elem{e};
			
			% Print header line stating: Eh, Area (and Gamma if visco-elastic)
			if (MESH.L(e)>MESH.elmt_shortL_vw) && MESH.Viscoeslastic == 1
				fprintf(fileIn,'%d %s %d %s\n',MESH.nbElemBySeg(e),'    nel Domain ',e, 'Eh Area Gamma');
			elseif (MESH.L(e)<=MESH.elmt_shortL_vw) || MESH.Viscoeslastic == 0
				fprintf(fileIn,'%d %s %d %s\n',MESH.nbElemBySeg(e),'    nel Domain ',e, 'Eh Area');
			else
				error('\n Unrecognised wall type')
			end	
			
			% cycle through elements
			for j=1:MESH.nbElemBySeg(e)
				
				% identify distances (x values)
				l_start = x_temp(j);
				l_end   = x_temp(j+1);
				
				% print mesh settings
				fprintf(fileIn,'%3.4f %3.4f %d %d %s\n', l_start, l_end, MESH.p(e), MESH.q(e), '# x_prox x_dist p q');
				
				% print area - IMPORTANT: print Area before Eh
				fprintf(fileIn,'%s %s %s %s\n','Area = PI*', Radius,'*',Radius);
				
				% print Eh
				fprintf(fileIn,'%s %s\n','Eh = ', Eh); %If prescribed with x-varying Area and Eh
				
				% print gamma if this segment is visco-elastic
				if (MESH.L(e)>MESH.elmt_shortL_vw) && MESH.Viscoeslastic == 1
					fprintf(fileIn,'%s %s\n','Gamma = ', Gamma); %Visco-elastic coefficient
				end
			end
		end
		
	case 'Empirical'
		% Write values of Area, Eh and (and Gamma if visco-elastic) for each arterial segment (each domain)
		for e=1:MESH.nbSeg
			
			% Print header line stating: Eh, Area (and Gamma if visco-elastic)
			fprintf(fileIn,'%d %s %d %s\n',MESH.nbElemBySeg(e),'    nel domain ',e, 'Empirical_I Area');  %If prescribed with empirical PWV and Area

			% Make string for radius
			Radius = ['(',num2str(MESH.Rad_coef(e,1)),'*x+',num2str(MESH.Rad_coef(e,2)),')'];
			
			% Print individual values for each element in this segment
			x_temp = MESH.x_elem{e};
			% cycle through elements
			for j=1:MESH.nbElemBySeg(e)
				
				% identify distances (x values)
				l_start = x_temp(j);
				l_end   = x_temp(j+1);
				
				% print mesh settings
				fprintf(fileIn,'%3.4f %3.4f %d %d %s\n', l_start, l_end, MESH.p(e), MESH.q(e), '# x_prox x_dist p q');
				
				% print area - IMPORTANT: print Area before a
				fprintf(fileIn,'%s %s %s %s\n','Ao = PI*', Radius,'*',Radius);
				
				% print a - viscoelastic coefficient from Raymond 2011 * 1.3 (Marie's suggestion)
				fprintf(fileIn,'%s %2.2f\n','a = ', MESH.PWV_a(e));
				
				% print gamma if this segment is visco-elastic
				if (MESH.L(e)>MESH.elmt_shortL_vw) && MESH.Viscoeslastic == 1
					error('Does this even make sense?')
					fprintf(fileIn,'%s %s\n','Gamma = ', Gamma); %Visco-elastic coefficient
				end
			end
		end
end


%%	3. BOUNDARY CONDITIONS
%	HEADER
fprintf(fileIn,'%s\n','Boundary conditions');

% cycle through each segment (domain)
for i=1:MESH.nbSeg

    %--- Inlet Boundary Conditions

    if i==1
%         %- Domain 1 : Specify inlet as reflective flow
		switch BCS.InflowType
			case 'F0'
        fprintf(fileIn,'%s \t%s\n','F 0',' # Domain 1'); %Frequency input (sum of harmonics in bcs file (F 0 or F 3))
        fprintf(fileIn,'%s\n','F  0');
			case 'q2'
		fprintf(fileIn,'%s \t%s\n','q 2',' # Domain 1'); %time-flow input
        fprintf(fileIn,'%s\n','q = 0');
		fprintf(fileIn,'%s\n','q 2');
        fprintf(fileIn,'%s\n','q = 0');
		
		% 		%- Domain 1 : Specify inlet as reflective flow
%         fprintf(fileIn,'%s \t%s\n','u  5',' # A lhs boundary Domain 1');
%         fprintf(fileIn,'%s\n','u = sin(3.1415*t/0.165)*step(0.165-t,0)');
% 		fprintf(fileIn,'%s \t%s\n','u  5',' # U lhs boundary');
%         fprintf(fileIn,'%s\n','u = sin(3.1415*t/0.165)*step(0.165-t,0)');
		end
        
    else
        %- all other domains: specify which segments they are connected to
        
        % properties of inlet
        inNode = MESH.node1(i);
        deg = MESH.degree(inNode);
        inConnect = MESH.EdgesByNode(inNode,:);
        
        % switch according to the type of junction at the inlet
        switch deg
            
            case 1 %no connectivity at inlet
                disp(['there is a problem at segment ',num2str(i)]);
                return;
                
            case 2 %junction
                ij = find(not(inConnect([1 2]) == i));
                fprintf(fileIn,'%s %d %d  \t %s%d\n','J',inConnect(ij),inConnect(ij), ' # Domain ',i);
                fprintf(fileIn,'%s %d %d\n','J',inConnect(ij),inConnect(ij));
                
            case 3 %bifurcation
                ib = find(not(inConnect == i));  % the segments which join to the current segment
                % determine type of bifurcation
                if (MESH.BifType(inNode)<0)
                    BifType = 'B';
                else
                    BifType = 'C';
                end
                fprintf(fileIn,'%s %d %d  \t %s%d\n',BifType,inConnect(ib(1)),inConnect(ib(2)), ' # Domain ',i);
                fprintf(fileIn,'%s %d %d\n',BifType,inConnect(ib(1)),inConnect(ib(2)));
        end
    end
    
    %--- Outlet Boundary Conditions
    
    % properties of outlet
    outNode = MESH.node2(i);
    deg = MESH.degree(outNode);
    outConnect = MESH.EdgesByNode(outNode,:);
    
    % switch according to the type of junction at the outlet
    switch deg
        
        case 1 %outlet WK
            % varies according to whether the terminal boundary conditions absorb reflections (purely reistance), or reflect (windkessel)
			
			% imaginary part
			if strcmp(BCS.WkType, 'W')
				fprintf(fileIn,'%s %1.4e\n', BCS.WkType, BCS.C_Wk(i)); %Windkessel Compliance
			else
				fprintf(fileIn,'%s %1.4e\n',BCS.WkType, BCS.R1_Wk(i)); % RT = Z0
			end
			
			% real part
			if isfield(BCS,'Scale_R1')
				fprintf(fileIn,'%s %1.4e %1.4f\n','W',BCS.R_Wk_tot(i), BCS.Scale_R1(i)); %Windkessel Resistance where R1 = Scale_R1*Z_0
			else
				fprintf(fileIn,'%s %1.4e\n','W',BCS.R_Wk_tot(i)); %Windkessel Resistance where R1 = Z0
			end
			
			% resistance for absorbing terminal
            if strcmp(BCS.WkType, 'R')
                fprintf(fileIn,'%s %1.4e\n','R',BCS.R1_Wk(i)); % RT = Z0
            end
            
        case 2 %junction
            ij = find(not(outConnect([1 2]) == i));
            fprintf(fileIn,'%s %d %d\n','J',outConnect(ij),outConnect(ij));
            fprintf(fileIn,'%s %d %d\n','J',outConnect(ij),outConnect(ij));
            
        case 3 %bifurcation
            ib = find(not(outConnect == i));  % the segments which join to the current segment
            % determine type of bifurcation
            if (MESH.BifType(outNode)<0)
                BifType = 'B';
            else
                BifType = 'C';
            end
            fprintf(fileIn,'%s %d %d\n',BifType,outConnect(ib(1)),outConnect(ib(2)));
            fprintf(fileIn,'%s %d %d\n',BifType,outConnect(ib(1)),outConnect(ib(2)));
    end
end


%%	4. INITIAL CONDITIONS
%	HEADER
fprintf(fileIn,'%s\n','Initial condition');

switch ICS.Type
	case 'A0 mesh'
		% cycle through each segment (domain)
		for i=1:MESH.nbSeg
			fprintf(fileIn,'%s\n','a = Ao');
			fprintf(fileIn,'%s\n','u = 0');
		end

	case 'A0 DBP'
		% cycle through each segment (domain)
		for e=1:MESH.nbSeg
			
			% Make string for radius
			Radius = ['(',num2str(MESH.Rad_coef(e,1)),'*x+',num2str(MESH.Rad_coef(e,2)),')'];
			
			% Make string for Eh
			Eh = [num2str(0.1),'*(',num2str(MESH.Eh_k1),'*exp(',num2str(MESH.Eh_k2),'*100*',Radius,') +',num2str(MESH.Eh_k3),')*',Radius];			
			
			% print initial area
			fprintf(fileIn,'%s %s %s %s %s %s %s %s %s %s %s\n','a = PI*',Radius,'*',Radius,'*(1-3/4*',Radius,'*',num2str(ICS.Pdiastolic),'/(',Eh,'))^2');
% 			else
% 				fprintf(fileIn_cond,'%s\n','a = Ao');
% 			end
			
			% print initial velocity
			fprintf(fileIn,'%s\n','u = 0');
		end
end


%%	5. HISTORY POINTS
%	HEADER
fprintf(fileIn,'%s\n','History Pts');
fprintf(fileIn,'%d %s\n',length(HIS.ArteryOutput),' #Number of Domains with history Points');

%- Outputs at inlet, middle and outlet of some segment
locs = cell(length(HIS.ArteryOutput),1);
for seg_no_el = 1 : length(HIS.ArteryOutput)
    curr_seg_no = HIS.ArteryOutput(seg_no_el);
    temp = [0, MESH.L(curr_seg_no)./2, MESH.L(curr_seg_no)];
    locs{seg_no_el} = unique(temp); clear temp
end

%- Outputs at inlet
% locs = cell(length(HIS.ArteryOutput),1);
% for seg_no_el = 1 : length(HIS.ArteryOutput)
%     temp = [0];
%     locs{seg_no_el} = unique(temp); clear temp
% end

% cycle through each segment at which to write output measurements
for seg_no_el=1:length(HIS.ArteryOutput)
    fprintf(fileIn,'%d %d\n',length(locs{seg_no_el}), HIS.ArteryOutput(seg_no_el));
    format_spec = [repmat('%1.4f ', [1,length(locs{seg_no_el})-1]), '%1.4f\n'];
    fprintf(fileIn,format_spec,locs{seg_no_el});
end

%- Pressure along aorta - output point every 0.005m
% his_aorta_art = [1 3 18 22 32 34 36 38 40 41 43 49 51 53 55];
% el = 0.005;
% for l=1:length(his_aorta_art)
%     his_loc = [[0:el:L(his_aorta_art(l))],L(his_aorta_art(l))];
%     his_ao{l} = his_loc;
% end
% fprintf(fileIn,'%s\n','History Pts');
% fprintf(fileIn,'%d %s\n',length(his_aorta_art),' #Number of Domains with history Points');
% for i=1:length(his_aorta_art)
%     fprintf(fileIn,'%d %d\n',length(his_ao{i}),his_aorta_art(i));
%     his_loc = his_ao{i};
%     for j = 1:length(his_ao{i})
%         fprintf(fileIn,'%1.4f ', his_loc(j));
%     end
%     fprintf(fileIn,'\n');
% end


fprintf(fileIn,'\n\n%s %s','oneDbio -a -R', [PATHS.Aortic1D_files,PATHS.PatientName]);


%% Close file
fclose(fileIn);

end