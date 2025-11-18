#Constant forces from FBD
Fx = 10 #N
Fy = 20 #N
Fz = 30 #N
Mx = 5  #Nm
My = 15 #Nm
Mz = 25 #Nm


def Number_Of_Fasteners(w, D_2, material, N_min):
    # x-z plane defined same as in 4.1

    # defining constraints

    # edge constraints are given as a range in 4.4
    if material == " 1 ":
        edge_center_min = [2, 3]
    elif material == " 2 ":
        edge_center_min = [4, 5]
    
    # fastener spacing always the same
    center_center_min = 1.5

    N_max = [] # Number of fasteners is a range, given the edge distance factor is also a range
    w_new = [] # Same for new height (calculated with minimum Number of fasteners)
    for e in range(0, len(edge_center_min)):

        something = 1 + (w-2*edge_center_min[e]*D_2)/(center_center_min*D_2) # unrounded N_max
        N_max.append(1 + (w-2*edge_center_min[e]*D_2)//(center_center_min*D_2)) 
        w_new.append((N_min-1) * center_center_min + 2 * edge_center_min[e])

        print("For chosen edge margin:", edge_center_min[e], "Max N, given w:", edge_center_min[e], N_max[e], something, "New w, given minimum N:", N_min, w_new[e])
    
    return N_max, w_new


# test

print(Number_Of_Fasteners(10, 2, " 1 ", 2))

