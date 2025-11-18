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

    # edge constraints have a range
    if material == " 1 "
        edge_center_min = [2, 3]
    elif material == " 2 "
        edge_center_min = [4, 5]
    
    center_center_min = 1.5

    for e in edge_center_min:
        
        N_max = (w-2*edge_center_min*D_2)//(center_center_min*D_2)
    print(N_max, (w-2*edge_center_min*D_2)/(center_center_min*D_2) )

    # if there is a minimum N, recalculate w to fit N fasteners with D_2 diameter
    if N_max < N_min:
        # 
    
    print("N_max, w:", N_max, w, "w_new given N")


    return #max N, margin for given D_2 that still meets constraints