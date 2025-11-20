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

    N_max_check = [] # Number of fasteners check
    for e in range(0, len(edge_center_min)):

        something = 1 + (w-2*edge_center_min[e]*D_2)/(center_center_min*D_2) # unrounded N_max
        N_m = 1 + (w-2*edge_center_min[e]*D_2)//(center_center_min*D_2)
        N_max_check.append(N_m)
    
        e_test = ( N_m - 1 - w/(center_center_min*D_2) ) * -center_center_min/2 # chooses highest N and its associated edge spacing
        if e_test > edge_center_min[e]:
            N_max = max(N_max_check)
            edge_spacing = e_test
            #print(edge_spacing, edge_center_min[e], N_max, "1") # for testing

        w_new = (N_min-1) * center_center_min + 2 * min(edge_center_min)
        #print(e_test, N_max, "2") # for testing

    print("Max N, given w:",N_max,",", edge_spacing,",", "New w, given minimum N:", w_new, N_min)
    
    return N_max, edge_spacing, center_center_min, w_new, N_min


# test

print(Number_Of_Fasteners(15, 2, " 1 ", 2))