function [pos_intersect,t_intersect] = stop_at_boundary_intersection(prev_pos,pos,direction,A)
           t_intersect = solve_ellipsoid_intersection(prev_pos,direction,A);
            if isnan(t_intersect)
                error('t is nan');
            end
            t_intersect = max(t_intersect); %want the first point it intersects, need this positive
            pos_intersect = prev_pos + t_intersect*direction;
end


