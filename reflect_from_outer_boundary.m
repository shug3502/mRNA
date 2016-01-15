function new_pos = reflect_from_outer_boundary(prev_pos,pos,direction,A)
           t_intersect = solve_ellipsoid_intersection(prev_pos,direction,A);
            if isnan(t_intersect)
                error('t is nan');
            end
            t_intersect = max(t_intersect); %want the first point it intersects, need this positive
            pos_intersect = prev_pos + t_intersect*direction;
            unit_normal = -diag(A)'*pos_intersect; %inward pointing normal?
            unit_normal = unit_normal./norm(unit_normal); %normalise

            Lstar = 2*pos_intersect - pos;
            new_pos = 2*(unit_normal'*Lstar)*unit_normal - Lstar; %new reflected position

                %checks on new position
                if abs(imag(t_intersect))>0
                        fprintf('xpos %f \n',new_pos);
                end
%                if new_pos'*A*new_pos >1
%                        res = new_pos'*A*new_pos -1;
%                        fprintf('Not reflected properly: res %f, t %f\n',res,t_intersect);
%                end
            if new_pos'*A*new_pos -1 > 10^-6 %some tolerance
                new_pos = pos_intersect; %take position as intersection pt. NB Not treated properly for second reflection
            end

end

