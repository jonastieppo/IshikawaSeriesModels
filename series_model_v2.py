'''
Model to implement the closed solutions for wove-fibers composites obtained by Ishikawa (1982)
'''
# %%
import numpy as np
import matplotlib.pyplot as plt
from CLT import*

class Series_Models():
    '''
    Implemeting Ishikawa Solution for series model
    '''
    def __init__(self):
        print("Ishikawa Solution Started!")

    def DefineEngConstants(self,dictionary_with_values):
        '''
        Method to define the Engineering Constants. The input must
        come in a dictionary form.
        
        Mandatory shape:

        Eng_const = {

            "E1":11500,
            "E2":900,
            "G12":455,
            "nu_12":0.3
        }
        '''
        self.EngineeringConstants = dictionary_with_values


    def ShowEngConstants(self):
        '''
        Method to print the Inserted Eng. Constants
        '''
        print(self.EngineeringConstants)


    def CheckIfEngWereInserted(self):
        '''
        Method to check if the Engineering contants were Inserted correctly
        '''
        right_k_words = ['E1','E2','G12','v12','v23']
        c = 0
        for k in self.EngineeringConstants:
            c=c+1
            if not k in right_k_words:
                return "error: Input the Eng. Correctly. Please run ShowEngConstants and check if everything is right."
            if c>5:
                return "error: More than 4 eng. constats were inserted.Please run ShowEngConstants and check if everything is right." 

    def CalculateZinfo(self,h,n=2):
        '''
        Method to calculate the z  location of the i-esim laminate
        '''
        z_cood = [-h/2]
        
        z_cood = [-h/2]

        for each_ply in range(n):
            z_cood.append(z_cood[each_ply]+h/n)

        self.z_cood=z_cood
        
    def WritePlateStiffness(self,h):
        '''
        Method to write the A,B,D coeffiencts following Ishikawa
        presented equations.
        '''
        self.h=h
        Qbar_0=Qbar(0,self.EngineeringConstants)
        Qbar_90=Qbar(90,self.EngineeringConstants)

        self.CalculateZinfo(h)

        z_info1 = z_info(self.z_cood[0],self.z_cood[1])
        z_info2 = z_info(self.z_cood[1],self.z_cood[2])

        Ninfo1 = N_info(Qbar_90,z_info1)
        Ninfo2 = N_info(Qbar_0,z_info2)

        Laminate_info = {
                        'n1':Ninfo1,
                        'n2':Ninfo2
                        }
        A = A_submatriz(Laminate_info)
        B = B_submatriz(Laminate_info)
        D = D_submatriz(Laminate_info)

        ABD = MountABD(A,B,D)

        self.plate_stiffness = ABD
        C = self.plate_stiffness
        self.plate_compliance = np.linalg.inv(self.plate_stiffness)

    def IsoForce(self,n):
        '''
        Method to write the modified compliance based on Series Model.
        n is the repetition unity for the weave cloth. 
        For the plain-weave, n=2. For eigh-harness sating weave, n=8
        '''
        factor = (1-2/n)
        Matrix_fator = np.array([[1,1,1,factor,factor,factor],
                                 [1,1,1,factor,factor,factor],
                                 [1,1,1,factor,factor,factor],
                                 [factor,factor,factor,1,1,1],
                                 [factor,factor,factor,1,1,1],
                                 [factor,factor,factor,1,1,1]])
       
        new_abd_inv = np.multiply(Matrix_fator,self.plate_compliance)
        self.Stifness_lower_iso_force = np.linalg.inv(new_abd_inv)
        self.Compliance_upper_iso_force = np.linalg.inv(self.Stifness_lower_iso_force)
        self.Ex_lower,self.Ey_lower,self.Gxy_lower,self.vxy_lower = CalculateEngineeringConstants(self.Stifness_lower_iso_force,self.h)

    def IsoStrain(self,n):
        '''
        Method to write the modiefied stiffness. It represents
        the upper bonds of stiffness
        '''
        factor = (1-2/n)
        Matrix_fator = np.array([[1,1,1,factor,factor,factor],
                                 [1,1,1,factor,factor,factor],
                                 [1,1,1,factor,factor,factor],
                                 [factor,factor,factor,1,1,1],
                                 [factor,factor,factor,1,1,1],
                                 [factor,factor,factor,1,1,1]])

        new_abd = np.multiply(Matrix_fator,self.plate_stiffness)
        self.Ex_upper,self.Ey_upper,self.Gxy_upper,self.vxy_upper = CalculateEngineeringConstants(new_abd,self.h)       
        self.Stifness_upper_iso_strain = new_abd
        self.Compliance_lower_iso_strain = np.linalg.inv(self.Stifness_upper_iso_strain)

# %%
'''
Testing the class
'''
Model = Series_Models()

Eng_const = {

    "E1":11500,
    "E2":900,
    "G12":455,
    "v12":0.3,
    "v23":0.3
}


Model.DefineEngConstants(Eng_const)
Model.WritePlateStiffness(0.4)
# %%
'''
Plotting the therma A11
'''
nlist = np.linspace(2,1000,1000)
y1 = []
y2 = []
for i in nlist:
    Model.IsoForce(i)
    Model.IsoStrain(i)
    y1.append(Model.Stifness_upper_iso_strain[0,0])
    y2.append(Model.Stifness_lower_iso_force[0,0])
    # y3.append(Model.Compliance_upper[0,0])
    # y4.append(Model.Compliance_lower[0,0])
fig,ax = plt.subplots()

nlist = [1/i for i in nlist]
ax.plot(nlist,y1,label='$A_{11}$ upper Bound (Iso Strain)')
ax.plot(nlist,y2,label='$A_{11}$ lower Bound (Iso Force)')
plt.xlabel('$1/n$')
plt.ylabel('$A_{11}$')
plt.legend()

plt.savefig("A11.pdf")
# %%
'''
Plotting the compliance coupling therm b*11
'''
nlist = np.linspace(2,1000,1000)
y1 = []
y2 = []
for i in nlist:
    Model.IsoForce(i)
    Model.IsoStrain(i)
    y1.append(Model.Compliance_upper_iso_force[3,0])
    y2.append(Model.Compliance_lower_iso_strain[3,0])
    # y3.append(Model.Compliance_upper[0,0])
    # y4.append(Model.Compliance_lower[0,0])
fig,ax = plt.subplots()

nlist = [1/i for i in nlist]
ax.plot(nlist,y1,label='$b^*_{11}$ upper Bound (Iso Force)')
ax.plot(nlist,y2,label='$b^*_{11}$ lower Bound (Iso Strain)')
plt.xlabel('$1/n$')
plt.ylabel('$b^*_{11}$')
plt.legend()

plt.savefig("b11.pdf")

# %%
