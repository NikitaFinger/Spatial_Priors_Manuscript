function v_int = int_track(x,v,x_int)
v_int(:,1) = interp1(x,v(:,1),x_int);
v_int(:,2) = interp1(x,v(:,2),x_int);
v_int(:,3) = interp1(x,v(:,3),x_int);
