function pot_orbit_weights(recompute,Parray,permutations)

filename = sprintf('localdata/orbits/pot/weights.mat');
if isfile(filename)&&~recompute

    fprintf('orbit POT weights already exists. \n');

else

    fprintf('computing POT weights...\n')

    P = numel(Parray);
    R = size(permutations,2);
    w = cell(P,1);
    str = '';
    for i = 1:P

        p = Parray(i);
        w{i} = zeros(p,R);

        for r = 1:R

            ind = permutations(1:p,r);
            w{i}(:,r) = compute.unordered_pot_weights(ind);

        end

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',i,P);
        fprintf(str);

    end

    % save out data
    save(filename,'w');
    fprintf('saved results to `%s`\n',filename)

end
