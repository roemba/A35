import numpy as np

class parameters:
    chord = 0.547 #C_a	0.547	m
    span = 2.771 #l_a	2.771	m
    xlocation1 = 0.153 #x_1	0.153	m
    xlocation2 = 1.281 #x 2	1.281	m
    xlocation3 = 2.681 #x_3	2.681	m
    d12 = 28.0/100. #d_12	x_a	28.0	cm
    height = 22.5/100. #Height	h	22.5	cm
    skinthickness = 1.1/1000. #Skin_thickness	t_sk	1.1	mm
    sparthickness = 2.9/1000. #Spar_thickness	t sp	2.9	mm
    stiffenerthickness = 1.2/1000. #Stiffener_thickness	t_st	1.2	mm
    stiffenerheight = 1.5/100. #Stiffener_Height	h_st	1.5	cm
    stiffenerwidth = 2.0/100. #Stiffener_Width	w_st	2.0	cm
    stiffenernumber = 17. #Stiffener_Number	n_st	17	-
    verticaldisplacementhinge1 = (11.03)*0.0254 #Vertical_displacement_hinge_1	d_1	11.03	cm
    verticaldisplacementhinge3 = (16.42)*0.0254 #Vertical_displacement_hinge_3	d_3	16.42	cm
    maxupwarddeflection = np.radians(26.) #Max_upward_deflection	q	26	deg
    actuatorload = 9.17e3 #Actuator_load	P	9.17	kN
    aerodynamicload = 4.53e3 #Aerodynamic_load	q	4.53	kN/m
    youngsmodulus = 73.1e9 #Youngs modulus E 73.1 GPa
