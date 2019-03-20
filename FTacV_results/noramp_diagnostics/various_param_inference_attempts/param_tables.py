import numpy as np
def latex_dict_plotter(dict):
    keys=dict.keys()
    num_rows=len(keys)
    top_row=' & '.join(keys)
    print top_row +" \\\\"
    print "\\hline"
    val_list=[str(param_list[x]) for x in keys]
    bottom_row=' & '.join(val_list)
    print bottom_row +" \\\\"
    print "\\hline"
param_list={
            '$E_{start}$': -0.85, #(starting dc voltage - V)
            '$E_{reverse}$': -0.1,    #  (reverse dc voltage - V)
            '$\\omega$':8.8,#8.88480830076,  #    (frequency Hz)
            '$\\DeltaE$': "150e-3",   #(ac voltage amplitude - V) freq_range[j],#
            'v': "24.97e-3",   #       (scan rate s^-1)
            'area': 0.03, #(electrode surface area cm^2)
            'Ru': 20.0,  #     (uncompensated resistance ohms)
            'C_{dl}': 0.000134, #(capacitance parameters)
            '$\\Gamma$': 6.5e-12,          # (surface coverage per unit area)
            '$E_0$': -0.4,      #       (reversible potential V)
            '$k_0$': 10.0, #(reaction rate s-1)
            '$\\alpha$': 0.5,
            'Sampling frequency' : (1.0/200),
        }
latex_dict_plotter(param_list)
