# Log likelihood = -18.26563974403333290297
# Number of coding sequence sites = 1

# Substitution model parameters:
model=TwoParameterBinary(mu=0.193055454634,pi0=0.950000000000)

nonhomogeneous=general
nonhomogeneous.number_of_models=2
nonhomogeneous.stationarity = yes
model1=RELAX(1_Full.theta=0.500000000000,1_Full.theta1=0.500000000000,1_Full.theta2=0.500000000000,2_Full.theta=0.500000000000,2_Full.theta1=0.500000000000,2_Full.theta2=0.500000000000,3_Full.theta=0.500000000000,3_Full.theta1=0.500000000000,3_Full.theta2=0.500000000000,kappa=4.006417673813,theta1=0.500000000000,theta2=0.800000000000,p=0.100000000000,omega1=1.000000000000,omega2=1.000000000000,k=1.000000000000)
model1.nodes_id=0,8,1,10,2,12,3,14,4,16,17

model2=RELAX(kappa=RELAX.kappa_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,k=0.000000000001)
model2.nodes_id=7,9,11,13,15,5

