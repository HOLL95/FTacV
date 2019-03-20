import  matplotlib.pyplot as plt
import pints.plot
import numpy as np

length_list=[2e4]
dec_list=[32]
repeat_num=5
for lcv_1 in range(0, len(length_list)):
    for lcv_2 in range(0, len(dec_list)):
        for lcv_3 in range(0, repeat_num):
            desired_length=int(length_list[lcv_1])
            dec_amount=dec_list[lcv_2]
            filename=str(desired_length)+"_"+str(dec_amount)+"_"+str(lcv_3)+".exp"
            chains=np.load(filename)
            mean_params=np.zeros(len(chains[0, 0,:]))
            for i in range(0, len(chains[0, 0,:])):
                mean_params[i]= np.mean(chains[:,5000:,i])
            print mean_params

            print filename
            pints.plot.trace(chains)
            plt.show()
"""
[2.40666629e-01, 1.33659203e+04, 8.32138327e+00, 6.48986896e-0,  6.99137283e-10 ]
 [2.89424312e-01, 1.42829102e+04, 1.16728294e+00, 2.95794563e-02, 1.54470074e-10 ]
 [2.37426876e-01, 1.40219177e+04, 1.19240444e+01, 2.29657002e-01, 1.42238592e-09 ]
 [2.41984331e-01, 1.03875610e+00, 2.04067494e+01, 5.02952291e-03, 3.85177644e-10 ]
[2.34130799e-01, 4.02039142e+01, 2.46522095e+02, 3.91211938e-04, 2.59101097e-10 ]
 [2.77009135e-01, 1.01677918e+04, 1.22365956e+01, 1.17780590e+00, 1.21045676e-09 ]
 [2.85687934e-01, 1.39493817e+04, 5.61819373e-01, 4.73712354e-02, 1.46875670e-10 ]
 [2.85658802e-01, 1.31128704e+04, 5.60503182e-01, 4.73852493e-02, 1.46876400e-10 ]
 [2.74746341e-01 1.57993170e+03, 1.46569159e+02, 3.42081253e-03, 6.72772614e-10]

 """
