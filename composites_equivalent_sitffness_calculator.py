# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:23:27 2023

@author: pawel
"""

import numpy as np 

class Q_matrix:
    def __init__(self, E11, E22, v12, G12):
        self.E11 = E11
        self.E22 = E22
        self.v12 = v12
        self.G12 = G12
    
    def calculated_Q(self):
        self.v21 = (self.v12)*(self.E22)/(self.E11)
        
        self.Q11 = (self.E11)/(1 - self.v12*self.v21)
        self.Q22 = (self.E22)/(1 - self.v12*self.v21)
        self.Q12 = (self.v12)*(self.E22)/(1 - self.v12*self.v21)
        self.Q66 = self.G12
        
        return self.Q11, self.Q22, self.Q12, self.Q66
    
    def transformed_Q_theta(self, theta):
        sin = np.sin(np.deg2rad(theta))        
        cos = np.cos(np.deg2rad(theta))

        self.Q11_theta = (self.Q11)*(cos**4) + 2*(self.Q12 + 2*self.Q66)*(cos**2)*(sin**2) + (self.Q22)*(sin**4)
        self.Q22_theta = (self.Q11)*(sin**4) + 2*(self.Q12 + 2*self.Q66)*(cos**2)*(sin**2) + (self.Q22)*(cos**4)
        
        
        self.Q16_theta = (self.Q11 - self.Q12 - 2*self.Q66)*(cos**3)*(sin) - (self.Q22 - self.Q12 - 2*self.Q66)*(cos)*(sin**3)
        self.Q26_theta = (self.Q11 - self.Q12 - 2*self.Q66)*(sin**3)*(cos) - (self.Q22 - self.Q12 - 2*self.Q66)*(sin)*(cos**3)
        
        self.Q12_theta = (self.Q12)*(cos**4 + sin**4) + (self.Q11 + self.Q22 - 4*self.Q66)*(cos**2)*(sin**2)
        self.Q66_theta = (self.Q11 + self.Q22 - 2*self.Q12 - 2*self.Q66 )*(cos**2)*(sin**2) + (self.Q66)*(cos**4 + sin**4)
        
        self.Q_matrix = [[self.Q11_theta, self.Q12_theta, self.Q16_theta],
                         [self.Q12_theta, self.Q22_theta, self.Q26_theta],
                         [self.Q16_theta, self.Q26_theta, self.Q66_theta]]

        return self.Q_matrix
    

class ABD_matrix:
    def __init__(self, Q_matrices, theta_list, thicknesses_list):
        if len(Q_matrices) != len(theta_list) != len(thicknesses_list):
            raise ValueError('Input lists must have the same length')        
        self.Q_matrices = Q_matrices
        self.theta_list = theta_list
        self.thicknesses_list = thicknesses_list

    def A_matrix(self):
        A_matrix_components = []
        for k in range(len(self.Q_matrices)):
            matrix_A = self.Q_matrices[k]
            scalar_A = self.thicknesses_list[k]
            result_A = [[scalar_A * element_A for element_A in row_A] for row_A in matrix_A]  
            A_matrix_components.append(result_A)
        
        A_matrix = np.array(A_matrix_components)
        self.A_matrix = np.sum(A_matrix, axis=0)
        return self.A_matrix


    def B_matrix(self):
        B_matrix_components = []
        for k in range(len(self.Q_matrices)):
            matrix_B = self.Q_matrices[k]
            total_thickness = sum(self.thicknesses_list)
            if k == 0:
                h_0 = total_thickness/2
                h_k = h_0 - self.thicknesses_list[k]
                z_k = (h_0 + h_k)/2
            else: 
                h_k_min1 = h_k - self.thicknesses_list[k]
                z_k = (h_k_min1 + h_k)/2
                h_k = h_k_min1
            scalar_B = self.thicknesses_list[k]*z_k
            result_B = [[scalar_B * element_B for element_B in row_B] for row_B in matrix_B]  
            B_matrix_components.append(result_B)
        
        B_matrix = np.array(B_matrix_components)
        self.B_matrix = np.sum(B_matrix, axis=0)
        return self.B_matrix
    
    def D_matrix(self):
        D_matrix_components = []
        for k in range(len(self.Q_matrices)):
            matrix_D = self.Q_matrices[k]
            total_thickness = sum(self.thicknesses_list)
            if k == 0:
                h_0 = total_thickness/2
                h_k = h_0 - self.thicknesses_list[k]
                z_k = (h_0 + h_k)/2
            else: 
                h_k_min1 = h_k - self.thicknesses_list[k]
                z_k = (h_k_min1 + h_k)/2
                h_k = h_k_min1
            scalar_D = ((self.thicknesses_list[k]**3)/12 + self.thicknesses_list[k]*z_k**2)
            result_D = [[scalar_D * element_D for element_D in row_D] for row_D in matrix_D]  
            D_matrix_components.append(result_D)
        
        D_matrix = np.array(D_matrix_components)
        self.D_matrix = np.sum(D_matrix, axis=0)
        return self.D_matrix 
    
    def equivalent_constants(self):
        matrix_A = self.A_matrix
        matrix_B = self.B_matrix
        matrix_D = self.D_matrix
        
        A11, A22, A12, A16, A26, A66 = matrix_A[0,0], matrix_A[1,1], matrix_A[0,1], \
            matrix_A[0,2], matrix_A[1,2], matrix_A[2,2]
            
        B11, B22, B12, B16, B26, B66 = matrix_B[0,0], matrix_B[1,1], matrix_B[0,1], \
            matrix_B[0,2], matrix_B[1,2], matrix_B[2,2]     
            
        D11, D22, D12, D16, D26, D66 = matrix_D[0,0], matrix_D[1,1], matrix_D[0,1], \
            matrix_D[0,2], matrix_D[1,2], matrix_D[2,2]    
            
        nominator = np.array([[A11, A12, A16,    B11, B12, B16],
                              [A12, A22, A26,    B12, B22, B26],
                              [A16, A26, A66,    B16, B26, B66],
                                 
                              [B11, B12, B16,    D11, D12, D16],
                              [B12, B22, B26,    D12, D22, D26],
                              [B16, B26, B66,    D16, D26, D66]])

        E_x_denominator = np.array([[A22, A26,    B12, B22, B26],
                                    [A26, A66,    B16, B26, B66],
                                 
                                    [B12, B16,    D11, D12, D16],
                                    [B22, B26,    D12, D22, D26],
                                    [B26, B66,    D16, D26, D66]])
        
        E_y_denominator = np.array([[A11, A16,    B11, B12, B16],
                                    [A16, A66,    B16, B26, B66],
                                 
                                    [B11, B16,    D11, D12, D16],
                                    [B12, B26,    D12, D22, D26],
                                    [B16, B66,    D16, D26, D66]])

        G_xy_denominator = np.array([[A11, A12,    B11, B12, B16],
                                     [A12, A22,    B12, B22, B26],
                                 
                                     [B11, B12,    D11, D12, D16],
                                     [B12, B22,    D12, D22, D26],
                                     [B16, B26,    D16, D26, D66]])
        
        
        v_xy_nominator = np.array([[A12, A26,    B12, B22, B26],
                                   [A16, A66,    B16, B26, B66],
                                 
                                   [B11, B16,    D11, D12, D16],
                                   [B12, B26,    D12, D22, D26],
                                   [B16, B66,    D16, D26, D66]])        



        v_xy_denominator = np.array([[A22, A26,    B12, B22, B26],
                                     [A26, A66,    B16, B26, B66],
                                 
                                     [B12, B16,    D11, D12, D16],
                                     [B22, B26,    D12, D22, D26],
                                     [B26, B66,    D16, D26, D66]])                
        
        v_yx_nominator = np.array([[A12, A16,    B11, B12, B16],
                                   [A16, A66,    B16, B26, B66],
                                 
                                   [B12, B16,    D11, D12, D16],
                                   [B22, B26,    D12, D22, D26],
                                   [B16, B66,    D16, D26, D66]])        
                
        v_yx_denominator = np.array([[A11, A16,    B12, B22, B26],
                                     [A16, A66,    B16, B26, B66],
                                 
                                     [B11, B16,    D11, D12, D16],
                                     [B12, B26,    D12, D22, D26],
                                     [B16, B66,    D16, D26, D66]])                
        
        self.E_x = (np.linalg.det(nominator))/((np.linalg.det(E_x_denominator)*(sum(self.thicknesses_list))))
        self.E_y = (np.linalg.det(nominator))/((np.linalg.det(E_y_denominator)*(sum(self.thicknesses_list))))
        self.G_xy = (np.linalg.det(nominator))/((np.linalg.det(G_xy_denominator)*(sum(self.thicknesses_list))))

        self.v_xy = np.linalg.det(v_xy_nominator)/np.linalg.det(v_xy_denominator)
        self.v_yx = np.linalg.det(v_yx_nominator)/np.linalg.det(v_yx_denominator)
        
        return self.E_x, self.E_y, self.G_xy, self.v_xy, self.v_yx
    

# list1 = [1, 2, 3, 4, 5]
# list2 = [10, 20, 30, 40, 50]

# result = np.sum(np.multiply(list1, list2))
# print(result)

E11 = 3500
E22 = 2800
v12 = 0.3
G12 = E11/(2*(1+v12))


thetas = [45, -45, 45, -45, -45, 45, -45, 45]
thickness = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]



Q_s = []

for i in range(len(thetas)):
    theta = thetas[i]
    q = 0
    q = Q_matrix(E11, E22, v12, G12)
    q.calculated_Q()
    Q_s.append(q.transformed_Q_theta(theta))



ABD = ABD_matrix(Q_s, thetas, thickness)

ABD.A_matrix()
ABD.B_matrix()
ABD.D_matrix()


# print(ABD.A_matrix())
# print(ABD.B_matrix())
# print(ABD.D_matrix())

print(ABD.equivalent_constants())

# q = Q_matrix(E11, E22, v12, G12)
# q.calculated_Q()

# theta = 0
# print(q.transformed_Q_theta(theta))
