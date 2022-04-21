# Modelling ventral and dorsal BG circuits
# Ventral Striatum Modifies Activity of Basal Ganglia Circuit: A Computational Model


# With this code BG circuits is modelled; single neurons, synaptic currents, nucleus and circuits levels. 
# The results are given with membrane potentials, synaptic currents, raster plots, 
# frequency analysis  (power spectrum density, frequency time plot) and local field potentials.
# Single neurons are modelled using Izhikevich neuron model and modified Izhikevich neuron model.


from brian2 import *
start_scope()

# Parameters
#neuron parameters
Vthr = 30 * mvolt
EL = -75 * mV

#synaptic parameters
tau_s = 1 * ms
we = 0.1 *amp/mV
wi = 0.1 *amp/mV
Vi = -90 * mV
Ve = 0 * mV
dly=(3+rand())*ms

### synaptic weigths
w_cse = 0.20 *we        
w_tse = 0.25 *we
w_vse = 3.00 *we
w_vsi = 3.00 *wi
w_sse = 3.00 *we
w_ssi = 3.00 *wi


par_percent=10

number_of_neurons_in_ACA_pyramid = 900
number_of_neurons_in_ACA_in = 100

number_of_neurons_in_PFC_pyramid = 900
number_of_neurons_in_PFC_in = 100

number_of_neurons_in_MC_pyramid = 900
number_of_neurons_in_MC_in = 100


number_of_neurons_in_nacc=450
number_of_neurons_in_msnd1_core_121 = 100
number_of_neurons_in_msnd2_core_122 = 100
number_of_neurons_in_msnd1_shell_123 = 100
number_of_neurons_in_msnd2_shell_124 = 100
number_of_neurons_in_nacc_in = 50
number_of_neurons_in_msnd1_caudate_221 = 200
number_of_neurons_in_msnd2_caudate_222 = 200
number_of_neurons_in_caudate_in = 50


number_of_neurons_in_bg=300
number_of_neurons_in_GPe=100

number_of_neurons_in_thl = 100
number_of_neurons_in_vta = 100




#############   ---- Poisson Groups ------ #####################

PG_ACA_pyramid_1011 = PoissonGroup(number_of_neurons_in_ACA_pyramid, 5 * Hz)  #nominal value: 5*Hz
PG_ACA_in_1012 = PoissonGroup(number_of_neurons_in_ACA_in, 50 * Hz)

PG_BG_1040 = PoissonGroup(number_of_neurons_in_bg, 1 * Hz)

PG_Ventral_THL_1051 = PoissonGroup(number_of_neurons_in_thl, 5 * Hz)


PG_PFC_pyramid_2011 = PoissonGroup(number_of_neurons_in_PFC_pyramid, 5 * Hz)  #nominal value: 5*Hz
PG_PFC_in_2012 = PoissonGroup(number_of_neurons_in_PFC_in, 50 * Hz)    #nominal value: 50*Hz

PG_BG_2040 = PoissonGroup(number_of_neurons_in_bg, 1 * Hz)

PG_Dorsal_THL_2051 = PoissonGroup(number_of_neurons_in_thl, 5 * Hz)



print('Equations')


eqs_dyn = """
dv/dt=(0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u/ms+I*mV/(amp*ms)+Is*mV/(amp*ms) : volt
du/dt=a*(b*v-u)/ms                                           : volt
I : amp
Is=ge*(Ve-v)+gi*(Vi-v) :  amp
dge/dt=-ge/tau_e	: amp/volt
dgi/dt=-gi/tau_i	: amp/volt
a : 1
b : 1
c : volt
d : volt
tau_e : second
tau_i : second
"""


tau_glu=2*ms
tau_DA=1.5*ms
tau_i_msn=tau_s

eqs_msn = """
dv/dt=(0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u/ms+I*mV/(amp*ms)+Is*mV/(amp*ms) : volt
du/dt=a*(b*v+k*mV-u)/ms                                           : volt
I : amp
I_Glu=g_glu*(Ve-v) :  amp
I_DA=g_DA*(V_DA-v) :  amp
I_Ach=g_Ach*(Vi-v)  :  amp
I_GABA=g_GABA*(Vi-v) :  amp
Is=I_Glu+I_DA+I_Ach+I_GABA :  amp
dg_glu/dt=-g_glu/tau_glu	: amp/volt
dg_DA/dt=-g_DA/tau_DA	: amp/volt
dg_Ach/dt=-g_Ach/tau_i_msn	: amp/volt
dg_GABA/dt=-g_GABA/tau_i_msn	: amp/volt
a : 1
b : 1
c : volt
d : volt
k : 1
V_DA : volt
"""


eqs_reset = '''
v = c
u = u+d
'''



ACA_Pyramid_111 = NeuronGroup(number_of_neurons_in_ACA_pyramid, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_ACA_pyramid):
    ACA_Pyramid_111.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    ACA_Pyramid_111.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    ACA_Pyramid_111.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    ACA_Pyramid_111.d[i] = 8*(100-par_percent+2*par_percent*rand())/100* mvolt 
    ACA_Pyramid_111.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    ACA_Pyramid_111.u[i] = (-14.5*((100-par_percent+2*par_percent*rand())/100))*mvolt
ACA_Pyramid_111.tau_e = tau_s
ACA_Pyramid_111.tau_i = tau_s
    

ACA_in_112 = NeuronGroup(number_of_neurons_in_ACA_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_ACA_in):
    ACA_in_112.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    ACA_in_112.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    ACA_in_112.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    ACA_in_112.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    ACA_in_112.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    ACA_in_112.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
ACA_in_112.tau_e = tau_s
ACA_in_112.tau_i = tau_s


###################################################
####### ----------- Ventral BG Groups -----#######
###################################################

####### ----------- NAcc -----#######

msnd1_core_121 =  NeuronGroup(number_of_neurons_in_msnd1_core_121, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd1_core_121):
    msnd1_core_121.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd1_core_121.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd1_core_121.c[i] = -62*((100-par_percent+2*par_percent*rand())/100) * mvolt         # NV:-52  
    msnd1_core_121.d[i] = 0.6*((100-par_percent+2*par_percent*rand())/100) * mvolt          # NV: 1.9
    msnd1_core_121.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd1_core_121.u[i] = 35*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd1_core_121.k[i] = 35*((100-2*par_percent+4*par_percent*rand())/100)   #parametre araligi: 0.02 - 0.07
msnd1_core_121.V_DA = 0*mV



msnd2_core_122= NeuronGroup(number_of_neurons_in_msnd2_core_122, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd2_core_122):
    msnd2_core_122.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd2_core_122.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd2_core_122.c[i] = -60*((100-par_percent+2*par_percent*rand())/100) * mvolt        # NV:-52
    msnd2_core_122.d[i] = 0.6*((100-par_percent+2*par_percent*rand())/100) * mvolt         #NV: 1.9
    msnd2_core_122.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd2_core_122.u[i] = 25*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd2_core_122.k[i] = 20*((100-2*par_percent+4*par_percent*rand())/100)
msnd2_core_122.V_DA = -90*mV


msnd1_shell_123= NeuronGroup(number_of_neurons_in_msnd1_shell_123, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd1_shell_123):
    msnd1_shell_123.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd1_shell_123.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd1_shell_123.c[i] = -56*((100-par_percent+2*par_percent*rand())/100) * mvolt      # NV:-52  
    msnd1_shell_123.d[i] = 0.4*((100-par_percent+2*par_percent*rand())/100) * mvolt     #NV: 1.9
    msnd1_shell_123.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd1_shell_123.u[i] = 35*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd1_shell_123.k[i] = 35*((100-2*par_percent+4*par_percent*rand())/100)
msnd1_shell_123.V_DA = 0*mV    
    


msnd2_shell_124= NeuronGroup(number_of_neurons_in_msnd2_shell_124, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd2_shell_124):
    msnd2_shell_124.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd2_shell_124.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd2_shell_124.c[i] = -55*((100-par_percent+2*par_percent*rand())/100) * mvolt        #NV: -52
    msnd2_shell_124.d[i] = 0.4*((100-par_percent+2*par_percent*rand())/100) * mvolt       #NV: 1.9   
    msnd2_shell_124.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd2_shell_124.u[i] = 25*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd2_shell_124.k[i] = 20*((100-2*par_percent+4*par_percent*rand())/100)
msnd2_shell_124.V_DA = -90*mV



nacc_in_125= NeuronGroup(number_of_neurons_in_nacc_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_nacc_in):
    nacc_in_125.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    nacc_in_125.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    nacc_in_125.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    nacc_in_125.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt
    nacc_in_125.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    nacc_in_125.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
nacc_in_125.tau_e = tau_s
nacc_in_125.tau_i = tau_s


####### ----------- Ventral Pallidal -----#######

Ventral_GPe_141 = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    Ventral_GPe_141.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    Ventral_GPe_141.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    Ventral_GPe_141.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Ventral_GPe_141.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Ventral_GPe_141.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Ventral_GPe_141.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100)*mV
Ventral_GPe_141.tau_e = tau_s
Ventral_GPe_141.tau_i = tau_s*2


Ventral_GPi_142 = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    Ventral_GPi_142.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    Ventral_GPi_142.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    Ventral_GPi_142.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Ventral_GPi_142.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Ventral_GPi_142.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Ventral_GPi_142.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100)*mV
Ventral_GPi_142.tau_e = tau_s
Ventral_GPi_142.tau_i = tau_s*2


Ventral_STN_143 = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    Ventral_STN_143.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    Ventral_STN_143.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    Ventral_STN_143.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Ventral_STN_143.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Ventral_STN_143.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Ventral_STN_143.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
Ventral_STN_143.tau_e = tau_s
Ventral_STN_143.tau_i = tau_s



VTA_DA_131 = NeuronGroup(number_of_neurons_in_vta, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_vta):
    VTA_DA_131.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    VTA_DA_131.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    VTA_DA_131.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    VTA_DA_131.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    VTA_DA_131.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    VTA_DA_131.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
VTA_DA_131.tau_e = tau_s
VTA_DA_131.tau_i = tau_s



Ventral_THL_151 = NeuronGroup(number_of_neurons_in_thl, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_thl):
    Ventral_THL_151.a[i] = 0.03*((100-par_percent+2*par_percent*rand())/100)
    Ventral_THL_151.b[i] = 0.25*((100-par_percent+2*par_percent*rand())/100)
    Ventral_THL_151.c[i] = -52*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Ventral_THL_151.d[i] = 0.01*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Ventral_THL_151.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Ventral_THL_151.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100) *mV
Ventral_THL_151.tau_e = tau_s
Ventral_THL_151.tau_i = tau_s*10



ACA_Pyramid_111.I = 0*amp
ACA_in_112.I = 0*amp

msnd1_core_121.I = 0*amp
msnd2_core_122.I = 0*amp
msnd1_shell_123.I = 0*amp
msnd2_shell_124.I = 0*amp
nacc_in_125.I = 0*amp

Ventral_GPe_141.I = 0*amp
Ventral_GPi_142.I = 0*amp
Ventral_STN_143.I = 0*amp

VTA_DA_131.I = 0*amp
Ventral_THL_151.I = 0*amp



###################################################
####### ----------- Dorsal BG Groups -----#######
###################################################

PFC_Pyramid_211 = NeuronGroup(number_of_neurons_in_PFC_pyramid, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_PFC_pyramid):
    PFC_Pyramid_211.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    PFC_Pyramid_211.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    PFC_Pyramid_211.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    PFC_Pyramid_211.d[i] = 8*(100-par_percent+2*par_percent*rand())/100* mvolt 
    PFC_Pyramid_211.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    PFC_Pyramid_211.u[i] = (-14.5*((100-par_percent+2*par_percent*rand())/100))*mvolt
PFC_Pyramid_211.tau_e = tau_s
PFC_Pyramid_211.tau_i = tau_s
    
    


PFC_in_212 = NeuronGroup(number_of_neurons_in_PFC_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_ACA_in):
    PFC_in_212.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    PFC_in_212.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    PFC_in_212.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    PFC_in_212.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    PFC_in_212.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    PFC_in_212.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
PFC_in_212.tau_e = tau_s
PFC_in_212.tau_i = tau_s


####### ----------- Caudate -----#######

msnd1_caudate_221 =  NeuronGroup(number_of_neurons_in_msnd1_caudate_221, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd1_caudate_221):
    msnd1_caudate_221.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd1_caudate_221.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd1_caudate_221.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt         # NV:-52  
    msnd1_caudate_221.d[i] = 0.6*((100-par_percent+2*par_percent*rand())/100) * mvolt          # NV: 1.9
    msnd1_caudate_221.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd1_caudate_221.u[i] = 35*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd1_caudate_221.k[i] = 35*((100-2*par_percent+4*par_percent*rand())/100)   #parametre araligi: 0.02 - 0.07
msnd1_caudate_221.V_DA = 0*mV



msnd2_caudate_222= NeuronGroup(number_of_neurons_in_msnd2_caudate_222, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd2_caudate_222):
    msnd2_caudate_222.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd2_caudate_222.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd2_caudate_222.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        # NV:-52
    msnd2_caudate_222.d[i] = 0.5*((100-par_percent+2*par_percent*rand())/100) * mvolt         #NV: 1.9
    msnd2_caudate_222.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd2_caudate_222.u[i] = 25*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd2_caudate_222.k[i] = 20*((100-2*par_percent+4*par_percent*rand())/100)
msnd2_caudate_222.V_DA = -90*mV


caudate_in_223= NeuronGroup(number_of_neurons_in_caudate_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_nacc_in):
    caudate_in_223.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    caudate_in_223.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    caudate_in_223.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    caudate_in_223.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt
    caudate_in_223.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    caudate_in_223.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
caudate_in_223.tau_e = tau_s
caudate_in_223.tau_i = tau_s


####### ----------- Dorsal Pallidal -----#######


Dorsal_GPe_241 = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    Dorsal_GPe_241.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_GPe_241.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_GPe_241.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Dorsal_GPe_241.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Dorsal_GPe_241.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_GPe_241.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100)*mV
Dorsal_GPe_241.tau_e = tau_s
Dorsal_GPe_241.tau_i = tau_s*2




Dorsal_GPi_242 = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    Dorsal_GPi_242.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_GPi_242.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_GPi_242.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Dorsal_GPi_242.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Dorsal_GPi_242.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_GPi_242.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100)*mV
Dorsal_GPi_242.tau_e = tau_s
Dorsal_GPi_242.tau_i = tau_s*2


Dorsal_STN_243 = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    Dorsal_STN_243.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_STN_243.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_STN_243.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Dorsal_STN_243.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Dorsal_STN_243.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_STN_243.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
Dorsal_STN_243.tau_e = tau_s
Dorsal_STN_243.tau_i = tau_s


SNc_DA_231 = NeuronGroup(number_of_neurons_in_vta, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_vta):
    SNc_DA_231.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    SNc_DA_231.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    SNc_DA_231.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    SNc_DA_231.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    SNc_DA_231.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    SNc_DA_231.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
SNc_DA_231.tau_e = tau_s
SNc_DA_231.tau_i = tau_s




Dorsal_THL_251 = NeuronGroup(number_of_neurons_in_thl, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_thl):
    Dorsal_THL_251.a[i] = 0.03*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_THL_251.b[i] = 0.25*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_THL_251.c[i] = -52*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    Dorsal_THL_251.d[i] = 0.01*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    Dorsal_THL_251.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    Dorsal_THL_251.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100) *mV

Dorsal_THL_251.tau_e = tau_s
Dorsal_THL_251.tau_i = tau_s*5



MC_Pyramid_311 = NeuronGroup(number_of_neurons_in_MC_pyramid, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_MC_pyramid):
    MC_Pyramid_311.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    MC_Pyramid_311.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    MC_Pyramid_311.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    MC_Pyramid_311.d[i] = 8*(100-par_percent+2*par_percent*rand())/100* mvolt 
    MC_Pyramid_311.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    MC_Pyramid_311.u[i] = (-14.5*((100-par_percent+2*par_percent*rand())/100))*mvolt
MC_Pyramid_311.tau_e = tau_s
MC_Pyramid_311.tau_i = tau_s
    
    

MC_in_312 = NeuronGroup(number_of_neurons_in_MC_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_MC_in):
    MC_in_312.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    MC_in_312.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    MC_in_312.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    MC_in_312.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    MC_in_312.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    MC_in_312.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
MC_in_312.tau_e = tau_s
MC_in_312.tau_i = tau_s


PFC_Pyramid_211.I = 0*amp
PFC_in_212.I = 0*amp

msnd1_caudate_221.I = 0*amp
msnd2_caudate_222.I = 0*amp
caudate_in_223.I = 0*amp

Dorsal_GPe_241.I = 0*amp
Dorsal_GPi_242.I = 0*amp
Dorsal_STN_243.I = 0*amp

SNc_DA_231.I = 0*amp
Dorsal_THL_251.I = 0*amp

MC_Pyramid_311.I = 0*amp
MC_in_312.I = 0*amp


###################################################
###### Synapses #############
#############################
print('Synapses')

################ -----  Ventral Synapses -------------- ######
## 111 ACA_Pyramid_111

S01_1011_111 = Synapses(PG_ACA_pyramid_1011, ACA_Pyramid_111,  'w :siemens', delay=dly, on_pre='ge += w')
S01_1011_111.connect(True, p = 0.25)

S02_112_111 = Synapses(ACA_in_112, ACA_Pyramid_111,  delay=dly, on_pre='gi += wi')
S02_112_111.connect(True, p = 0.25)

#thalamocortical synapses

S55_151_111 = Synapses(Ventral_THL_151, ACA_Pyramid_111,  'w :siemens', delay=dly, on_pre='ge += w')
S55_151_111.connect(True, p = 0.25)
S55_151_111.w=we*.25

# 112 IN

S03_1012_112 = Synapses(PG_ACA_in_1012, ACA_in_112,  delay=dly, on_pre='ge += 5*we')
S03_1012_112.connect(True, p = 0.25)

S04_111_112 = Synapses(ACA_Pyramid_111, ACA_in_112, delay=dly, on_pre='ge += we')
S04_111_112.connect(True, p = 0.25)

# 2 NAcc
# 121 Core MSND1

S05_111_121 = Synapses(ACA_Pyramid_111, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_glu += w')
S05_111_121.connect(True, p = 0.25)
S05_111_121.w=w_cse

S06_125_121 = Synapses(nacc_in_125, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_Ach += w')
S06_125_121.connect(True, p = 0.25)
S06_125_121.w=wi

S07_151_121 = Synapses( Ventral_THL_151, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_glu += w')
S07_151_121.connect(True, p = 0.25)
S07_151_121.w=w_tse

S08_131_121 = Synapses(VTA_DA_131, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_DA += w')
S08_131_121.connect(True, p = 0.25)
S08_131_121.w=w_vse

##### Colateral inhibitions from MSNs

S09_121_121 = Synapses(msnd1_core_121, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S09_121_121.connect(True, p = 0.05)
S09_121_121.w=wi

S10_122_121 = Synapses(msnd2_core_122, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S10_122_121.connect(True, p = 0.25)
S10_122_121.w=wi*2

S11_123_121 = Synapses(msnd1_shell_123, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S11_123_121.connect(True, p = 0.05)
S11_123_121.w=wi

S12_124_121 = Synapses(msnd2_shell_124, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S12_124_121.connect(True, p = 0.25)
S12_124_121.w=wi*2

# 122 Core MSND2

S13_111_122 = Synapses(ACA_Pyramid_111, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_glu += w')
S13_111_122.connect(True, p = 0.25)
S13_111_122.w=w_cse

S14_125_122 = Synapses(nacc_in_125, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S14_125_122.connect(True, p = 0.25)
S14_125_122.w=wi

S15_151_122 = Synapses( Ventral_THL_151, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_glu += w')
S15_151_122.connect(True, p = 0.25)
S15_151_122.w=w_tse

S16_131_122 = Synapses(VTA_DA_131, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_DA += w')
S16_131_122.connect(True, p = 0.25)
S16_131_122.w=w_vsi

##### Colateral inhibitions from MSNs

S17_121_122 = Synapses(msnd1_core_121, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S17_121_122.connect(True, p = 0.25)
S17_121_122.w=wi*2

S18_122_122 = Synapses(msnd2_core_122, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S18_122_122.connect(True, p = 0.05)
S18_122_122.w=wi

S19_123_122 = Synapses(msnd1_shell_123, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S19_123_122.connect(True, p = 0.25)
S19_123_122.w=wi*2

S20_124_122 = Synapses(msnd2_shell_124, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S20_124_122.connect(True, p = 0.05)
S20_124_122.w=wi

# 123 Shell MSND1

S21_111_123 = Synapses(ACA_Pyramid_111, msnd1_shell_123,  'w :siemens', delay=dly, on_pre='g_glu += w')
S21_111_123.connect(True, p = 0.25)
S21_111_123.w=w_cse

S22_125_123 = Synapses(nacc_in_125, msnd1_shell_123,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S22_125_123.connect(True, p = 0.25)
S22_125_123.w=wi

S23_151_123 = Synapses(Ventral_THL_151, msnd1_shell_123,  'w :siemens', delay=dly, on_pre='g_glu += w')
S23_151_123.connect(True, p = 0.25)
S23_151_123.w=w_tse

S24_131_123 = Synapses(VTA_DA_131, msnd1_shell_123,  'w :siemens', delay=dly, on_pre='g_DA += w')
S24_131_123.connect(True, p = 0.25)
S24_131_123.w=w_vse

##### Colateral inhibitions from MSNs

S25_121_123 = Synapses(msnd1_core_121, msnd1_shell_123, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S25_121_123.connect(True, p = 0.05)
S25_121_123.w=wi

S26_122_123 = Synapses(msnd2_core_122, msnd1_shell_123, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S26_122_123.connect(True, p = 0.25)
S26_122_123.w=wi*2

S27_123_123 = Synapses(msnd1_shell_123, msnd1_shell_123, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S27_123_123.connect(True, p = 0.05)
S27_123_123.w=wi

S28_124_123 = Synapses(msnd2_shell_124, msnd1_shell_123, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S28_124_123.connect(True, p = 0.25)
S28_124_123.w=wi*2

# 124 Shell MSND2

S29_111_124 = Synapses(ACA_Pyramid_111, msnd2_shell_124,  'w :siemens', delay=dly, on_pre='g_glu += w')
S29_111_124.connect(True, p = 0.25)
S29_111_124.w=w_cse

S30_125_124 = Synapses(nacc_in_125, msnd2_shell_124,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S30_125_124.connect(True, p = 0.25)
S30_125_124.w=wi

S31_151_124 = Synapses(Ventral_THL_151, msnd2_shell_124,  'w :siemens', delay=dly, on_pre='g_glu += w')
S31_151_124.connect(True, p = 0.25)
S31_151_124.w=w_tse

S32_131_124 = Synapses(VTA_DA_131, msnd2_shell_124,  'w :siemens', delay=dly, on_pre='g_DA += w')
S32_131_124.connect(True, p = 0.25)
S32_131_124.w=w_vsi

##### Colateral inhibitions from MSNs

S33_121_124 = Synapses(msnd1_core_121, msnd2_shell_124, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S33_121_124.connect(True, p = 0.25)
S33_121_124.w=wi*2

S34_122_124 = Synapses(msnd2_core_122, msnd2_shell_124, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S34_122_124.connect(True, p = 0.05)
S34_122_124.w=wi

S35_123_124 = Synapses(msnd1_shell_123, msnd2_shell_124, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S35_123_124.connect(True, p = 0.25)
S35_123_124.w=wi*2

S36_124_124 = Synapses(msnd2_shell_124, msnd2_shell_124, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S36_124_124.connect(True, p = 0.05)
S36_124_124.w=wi

# 125 IN

S37_111_125 = Synapses(ACA_Pyramid_111, nacc_in_125,  delay=dly, on_pre='ge += 0.25*we')
S37_111_125.connect(True, p = 0.20)

S38_121_125 = Synapses(msnd1_core_121, nacc_in_125, delay=dly, on_pre='gi += wi')
S38_121_125.connect(True, p = 0.25)

S39_122_125 = Synapses(msnd2_core_122, nacc_in_125, delay=dly, on_pre='gi += wi')
S39_122_125.connect(True, p = 0.25)

S40_123_125 = Synapses(msnd1_shell_123, nacc_in_125, delay=dly, on_pre='gi += wi')
S40_123_125.connect(True, p = 0.25)

S41_124_125 = Synapses(msnd2_shell_124, nacc_in_125, delay=dly, on_pre='gi += wi')
S41_124_125.connect(True, p = 0.25)

# 141 Ventral_GPe_141

S42_1040_141 = Synapses(PG_BG_1040, Ventral_GPe_141,  delay=dly, on_pre='ge += 2*we')
S42_1040_141.connect(True, p = 0.25)

S43_122_141 = Synapses(msnd2_core_122, Ventral_GPe_141,  delay=dly, on_pre='gi += wi')
S43_122_141.connect(True, p = 0.25)

S44_124_141 = Synapses(msnd2_shell_124, Ventral_GPe_141,  delay=dly, on_pre='gi += wi')
S44_124_141.connect(True, p = 0.25)

S45_143_141 = Synapses(Ventral_STN_143, Ventral_GPe_141,  delay=dly, on_pre='ge += we')
S45_143_141.connect(True, p = 0.25)

# 142 Ventral_GPi_142

S46_1040_142 = Synapses(PG_BG_1040, Ventral_GPi_142,  delay=dly, on_pre='ge += 2*we')
S46_1040_142.connect(True, p = 0.25)

S47_121_142 = Synapses(msnd1_core_121, Ventral_GPi_142,  delay=dly, on_pre='gi += wi')
S47_121_142.connect(True, p = 0.25)

S48_123_142 = Synapses(msnd1_shell_123, Ventral_GPi_142,  delay=dly, on_pre='gi += wi')
S48_123_142.connect(True, p = 0.25)

S49_141_142 = Synapses(Ventral_GPe_141, Ventral_GPi_142,  delay=dly, on_pre='gi += wi')
S49_141_142.connect(True, p = 0.25)

S93_143_142 = Synapses(Ventral_STN_143, Ventral_GPi_142,  delay=dly, on_pre='ge += .1*wi')
S93_143_142.connect(True, p = 0.25)

# 143 Ventral_STN_143

S50_111_143 = Synapses(ACA_Pyramid_111, Ventral_STN_143,  delay=dly, on_pre='ge += 0.5*we')
S50_111_143.connect(True, p = 0.25)

S51_141_143 = Synapses(Ventral_GPe_141, Ventral_STN_143,  delay=dly, on_pre='gi += wi')
S51_141_143.connect(True, p = 0.25)

S94_142_143 = Synapses(Ventral_GPi_142, Ventral_STN_143,  delay=dly, on_pre='gi += .1*wi')
S94_142_143.connect(True, p = 0.25)

# 151 Ventral_THL_151

S52_1051_151 = Synapses(PG_Ventral_THL_1051, Ventral_THL_151,  delay=dly, on_pre='ge += 0.20*we')    #agirlik 0.20 ye indirildi. 20200429 re
S52_1051_151.connect(True, p = 0.25)

S53_142_151 = Synapses(Ventral_GPi_142, Ventral_THL_151, delay=dly, on_pre='gi += 2.5*wi')
S53_142_151.connect(True, p = 0.5)

# 131 VTA

S54_1011_131 = Synapses(PG_ACA_pyramid_1011, VTA_DA_131,  'w :siemens', delay=dly, on_pre='ge += w')
S54_1011_131.connect(True, p = 0.25)
S54_1011_131.w=we*0.5


S92_211_131 = Synapses(PFC_Pyramid_211, VTA_DA_131,  'w :siemens', delay=dly, on_pre='ge += w')
S92_211_131.connect(True, p = 0.25)
S92_211_131.w=we*0.5  




###################################################
################ -----  Dorsal Synapses -------------- ######
## 211 PFC_Pyramid_211

S56_2011_211 = Synapses(PG_PFC_pyramid_2011, PFC_Pyramid_211,  'w :siemens', delay=dly, on_pre='ge += w')
S56_2011_211.connect(True, p = 0.25)
S56_2011_211.w=we*0

S57_212_211 = Synapses(PFC_in_212, PFC_Pyramid_211,  delay=dly, on_pre='gi += wi')
S57_212_211.connect(True, p = 0.25)

#thalamocortical synapses

S58_151_211 = Synapses(Ventral_THL_151, PFC_Pyramid_211,  'w :siemens', delay=dly, on_pre='ge += w')
S58_151_211.connect(True, p = 0.5)                                       # olasilik 0.5 e cikarildi re 20200429
S58_151_211.w=we*0.25

S59_251_211 = Synapses(Dorsal_THL_251, PFC_Pyramid_211,  'w :siemens', delay=dly, on_pre='ge += w')
S59_251_211.connect(True, p = 0.25)
S59_251_211.w=we*0.25

# 212 IN

S60_2012_212 = Synapses(PG_PFC_in_2012, PFC_in_212,  delay=dly, on_pre='ge += 5*we')
S60_2012_212.connect(True, p = 0.25)

S61_211_212 = Synapses(PFC_Pyramid_211, PFC_in_212, delay=dly, on_pre='ge += we')
S61_211_212.connect(True, p = 0.25)




# 2 Caudate
# 221 Caudate MSND1
dorsal_w_factor=1.1
S62_211_221 = Synapses(PFC_Pyramid_211, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_glu += w')
S62_211_221.connect(True, p = 0.25)
S62_211_221.w=w_cse*dorsal_w_factor

S63_223_221 = Synapses(caudate_in_223, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_Ach += w')
S63_223_221.connect(True, p = 0.25)
S63_223_221.w=wi

S64_251_221 = Synapses( Dorsal_THL_251, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_glu += w')
S64_251_221.connect(True, p = 0.25)
S64_251_221.w=w_tse*dorsal_w_factor

S65_231_221 = Synapses(SNc_DA_231, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_DA += w')
S65_231_221.connect(True, p = 0.25)
S65_231_221.w=w_vse*dorsal_w_factor

##### Colateral inhibitions from MSNs

S66_221_221 = Synapses(msnd1_caudate_221, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S66_221_221.connect(True, p = 0.05)
S66_221_221.w=wi

S67_222_221 = Synapses(msnd2_caudate_222, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S67_222_221.connect(True, p = 0.25)
S67_222_221.w=wi*2

# 222 Caudate MSND2

S68_211_222 = Synapses(PFC_Pyramid_211, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_glu += w')
S68_211_222.connect(True, p = 0.25)
S68_211_222.w=w_cse*dorsal_w_factor

S69_223_222 = Synapses(caudate_in_223, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S69_223_222.connect(True, p = 0.25)
S69_223_222.w=wi

S70_251_222 = Synapses(Dorsal_THL_251, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_glu += w')
S70_251_222.connect(True, p = 0.25)
S70_251_222.w=w_tse*dorsal_w_factor

S71_231_222 = Synapses(SNc_DA_231, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_DA += w')
S71_231_222.connect(True, p = 0.25)
S71_231_222.w=w_vsi*dorsal_w_factor

##### Colateral inhibitions from MSNs

S72_221_222 = Synapses(msnd1_caudate_221, msnd2_caudate_222, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S72_221_222.connect(True, p = 0.25)
S72_221_222.w=wi*2

S73_222_222 = Synapses(msnd2_caudate_222, msnd2_caudate_222, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S73_222_222.connect(True, p = 0.05)
S73_222_222.w=wi

# 223 Caudate IN

S74_211_223 = Synapses(PFC_Pyramid_211, caudate_in_223,  delay=dly, on_pre='ge += 0.25*we')
S74_211_223.connect(True, p = 0.20)

S75_221_223 = Synapses(msnd1_caudate_221, caudate_in_223, delay=dly, on_pre='gi += wi')
S75_221_223.connect(True, p = 0.25)

S76_222_223 = Synapses(msnd2_caudate_222, caudate_in_223, delay=dly, on_pre='gi += wi')
S76_222_223.connect(True, p = 0.25)

# 31 Dorsal_GPe_241

S77_2040_241 = Synapses(PG_BG_2040, Dorsal_GPe_241,  delay=dly, on_pre='ge += 2*we')
S77_2040_241.connect(True, p = 0.25)

S78_222_241 = Synapses(msnd2_caudate_222, Dorsal_GPe_241,  delay=dly, on_pre='gi += wi')
S78_222_241.connect(True, p = 0.25)

S79_243_241 = Synapses(Dorsal_STN_243, Dorsal_GPe_241,  delay=dly, on_pre='ge += we')
S79_243_241.connect(True, p = 0.25)

# 32 Dorsal_GPi_242

S80_2040_242 = Synapses(PG_BG_2040, Dorsal_GPi_242,  delay=dly, on_pre='ge += 2*we')
S80_2040_242.connect(True, p = 0.25)

S81_221_242 = Synapses(msnd1_caudate_221, Dorsal_GPi_242,  delay=dly, on_pre='gi += wi')
S81_221_242.connect(True, p = 0.25)

S82_241_242 = Synapses(Dorsal_GPe_241, Dorsal_GPi_242,  delay=dly, on_pre='gi += wi')
S82_241_242.connect(True, p = 0.25)

# 33 Dorsal_STN_243

S83_211_243 = Synapses(PFC_Pyramid_211, Dorsal_STN_243,  delay=dly, on_pre='ge += 0.5*we')
S83_211_243.connect(True, p = 0.25)

S84_241_243 = Synapses(Dorsal_GPe_241, Dorsal_STN_243,  delay=dly, on_pre='gi += wi')
S84_241_243.connect(True, p = 0.25)

# 41 Dorsal_THL_251

S85_2051_251 = Synapses(PG_Dorsal_THL_2051, Dorsal_THL_251,  delay=dly, on_pre='ge += 0.15*we')
S85_2051_251.connect(True, p = 0.25)

S86_242_251 = Synapses(Dorsal_GPi_242, Dorsal_THL_251, delay=dly, on_pre='gi += wi')
S86_242_251.connect(True, p = 0.25)

# 231 SNc

S87_2011_231 = Synapses(PG_PFC_pyramid_2011, SNc_DA_231,  'w :siemens', delay=dly, on_pre='ge += w')
S87_2011_231.connect(True, p = 0.25)
S87_2011_231.w=we 

###############################################################################
############## ---------------- Motor Cortex -----------------------  #########
###############################################################################

#### MC 311

S88_312_311 = Synapses(MC_in_312, MC_Pyramid_311,  delay=dly, on_pre='gi += wi')
S88_312_311.connect(True, p = 0.25)

#thalamocortical synapses

S89_251_311 = Synapses(Dorsal_THL_251, MC_Pyramid_311,  'w :siemens', delay=dly, on_pre='ge += w')
S89_251_311.connect(True, p = 0.25)
S89_251_311.w=we


# 312 IN

S90_2012_312 = Synapses(PG_PFC_in_2012, MC_in_312,  delay=dly, on_pre='ge += 5*we')
S90_2012_312.connect(True, p = 0.25)

S91_311_312 = Synapses(MC_Pyramid_311, MC_in_312, delay=dly, on_pre='ge += we')
S91_311_312.connect(True, p = 0.25)

#==============================================================================
#==============================================================================

import time
init_time=time.time()

# Simulation will begin now after defining everything needed to for neurons and the neural structures.
# But first 100ms of the simulation is just for settling all neurons to equilibrium. 
# So this part will not be used in analysis.
######### First 100ms #########
################################
duration1=100*ms

##Ventral Poisson
S01_1011_111.w=we*0
S54_1011_131.w=we*0
##Dorsal Poisson
S56_2011_211.w=we*0         #Poisson to PFC
S87_2011_231.w=we*0         #Poisson to SNc


print("----------------------------------------------")
print('100 ms initial conditions delay')
print("sim_time="+str(duration1))
run(duration1, report='text')


##Ventral Poisson
S01_1011_111.w=we*.75      ############# ACA Cortex
S54_1011_131.w=we     ############# VTA DA
##Dorsal Poisson
S56_2011_211.w=we*0      ############# PFC Cortex
S87_2011_231.w=we*0.1      ############# SNc DA


###############################################################################
##################### ---- Monitors ---- ######################################
###############################################################################
print('Monitors')

############ Ventral ####################
trace_ACA_Pyramid_111 = StateMonitor(ACA_Pyramid_111, 'v', record=9)
spikes_ACA_Pyramid_111 = SpikeMonitor(ACA_Pyramid_111)

trace_ACA_in_112 = StateMonitor(ACA_in_112, 'v', record=9)
spikes_ACA_in_112 = SpikeMonitor(ACA_in_112)

trace_msnd1_core_121 = StateMonitor(msnd1_core_121, 'v', record=True)
spikes_msnd1_core_121 = SpikeMonitor(msnd1_core_121)
trace_I_s_msnd1_core = StateMonitor(msnd1_core_121, 'Is', record=True)

trace_msnd2_core_122 = StateMonitor(msnd2_core_122, 'v', record=True)
spikes_msnd2_core_122 = SpikeMonitor(msnd2_core_122)
trace_I_s_msnd2_core = StateMonitor(msnd2_core_122, 'Is', record=True)

trace_msnd1_shell_123 = StateMonitor(msnd1_shell_123, 'v', record=True)
spikes_msnd1_shell_123 = SpikeMonitor(msnd1_shell_123)
trace_I_s_msnd1_shell = StateMonitor(msnd1_shell_123, 'Is', record=True)

trace_msnd2_shell_124 = StateMonitor(msnd2_shell_124, 'v', record=True)
spikes_msnd2_shell_124 = SpikeMonitor(msnd2_shell_124)
trace_I_s_msnd2_shell = StateMonitor(msnd2_shell_124, 'Is', record=True)

trace_nacc_in_125 = StateMonitor(nacc_in_125, 'v', record=True)
spikes_nacc_in_125 = SpikeMonitor(nacc_in_125)
trace_I_s_nacc_in = StateMonitor(nacc_in_125, 'Is', record=True)

trace_Ventral_GPe_141 = StateMonitor(Ventral_GPe_141, 'v', record=9)
spikes_Ventral_GPe_141 = SpikeMonitor(Ventral_GPe_141)

trace_Ventral_GPi_142 = StateMonitor(Ventral_GPi_142, 'v', record=9)
spikes_Ventral_GPi_142 = SpikeMonitor(Ventral_GPi_142)

trace_Ventral_STN_143 = StateMonitor(Ventral_STN_143, 'v', record=9)
spikes_Ventral_STN_143 = SpikeMonitor(Ventral_STN_143)

trace_VTA_DA_131 = StateMonitor(VTA_DA_131, 'v', record=9)
spikes_VTA_DA_131 = SpikeMonitor(VTA_DA_131)

trace_Ventral_THL_151 = StateMonitor(Ventral_THL_151, 'v', record=9)
spikes_Ventral_THL_151 = SpikeMonitor(Ventral_THL_151)

###############################################################################
########### Dorsal monitors ###################

trace_PFC_Pyramid_211 = StateMonitor(PFC_Pyramid_211, 'v', record=9)
spikes_PFC_Pyramid_211 = SpikeMonitor(PFC_Pyramid_211)

trace_MC_Pyramid_311 = StateMonitor(MC_Pyramid_311, 'v', record=9)
spikes_MC_Pyramid_311 = SpikeMonitor(MC_Pyramid_311)

trace_PFC_in_212 = StateMonitor(PFC_in_212, 'v', record=9)
spikes_PFC_in_212 = SpikeMonitor(PFC_in_212)

trace_msnd1_caudate_221 = StateMonitor(msnd1_caudate_221, 'v', record=True)
spikes_msnd1_caudate_221 = SpikeMonitor(msnd1_caudate_221)
trace_I_s_msnd1_caudate = StateMonitor(msnd1_caudate_221, 'Is', record=True)

trace_msnd2_caudate_222 = StateMonitor(msnd2_caudate_222, 'v', record=True)
spikes_msnd2_caudate_222 = SpikeMonitor(msnd2_caudate_222)
trace_I_s_msnd2_caudate = StateMonitor(msnd2_caudate_222, 'Is', record=True)

trace_caudate_in_223 = StateMonitor(caudate_in_223, 'v', record=True)
spikes_caudate_in_223 = SpikeMonitor(caudate_in_223)
trace_I_s_caudate_in = StateMonitor(caudate_in_223, 'Is', record=True)

trace_Dorsal_GPe_241 = StateMonitor(Dorsal_GPe_241, 'v', record=9)
spikes_Dorsal_GPe_241 = SpikeMonitor(Dorsal_GPe_241)

trace_Dorsal_GPi_242 = StateMonitor(Dorsal_GPi_242, 'v', record=9)
spikes_Dorsal_GPi_242 = SpikeMonitor(Dorsal_GPi_242)

trace_Dorsal_STN_243 = StateMonitor(Dorsal_STN_243, 'v', record=9)
spikes_Dorsal_STN_243 = SpikeMonitor(Dorsal_STN_243)

trace_SNc_DA_231 = StateMonitor(SNc_DA_231, 'v', record=9)
spikes_SNc_DA_231 = SpikeMonitor(SNc_DA_231)

trace_Dorsal_THL_251 = StateMonitor(Dorsal_THL_251, 'v', record=9)
spikes_Dorsal_THL_251 = SpikeMonitor(Dorsal_THL_251)

###############################################################################
###############################################################################
print("----------------------------------------------")
print("start")

number_of_scenario=0

#########################  ----- Scenarios   --------------  ####################

##-------------------- Scenario 0 --------------------------------##
## It is used for testing purposes.
# The scenario is used to test whether the code works.


if number_of_scenario==0:

    
    print("----------------------------------------------")
    print('Scenario 0')
    run(250*ms,report='text')



##-------------------- Scenario 1 --------------------------------##

###### Changing cortex and VTA input, investigate the relation of stimulus and reward

elif number_of_scenario==1:

    print("----------------------------------------------")
    print('Scenario 1')
    
    print('Scenario 1: D1 Resting')  
    duration1=100*ms
    duration2=500*ms
    duration3=1000*ms
    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0.1
    S92_211_131.w=we*0
    S87_2011_231.w=we*0.1      
            
    
    S05_111_121.w=w_cse*0        # from ACA_Pyramid_111 to msnd1 core
    S07_151_121.w=w_tse*0.1        # from Ventral_THL_151 to     msnd1 core
    S08_131_121.w=w_vse*0        # from VTA to     msnd1 core
    
    
    S13_111_122.w=w_cse*0        # from ACA_Pyramid_111 to msnd2 core
    S15_151_122.w=w_tse*0.1        # from Ventral_THL_151 to     msnd2 core
    S16_131_122.w=w_vsi*0        # from VTA to     msnd2 core
    
    
    S21_111_123.w=w_cse*0        # from ACA_Pyramid_111 to msnd1 shell
    S23_151_123.w=w_tse*0.1        # from Ventral_THL_151 to     msnd1 shell
    S24_131_123.w=w_vse*0        # from VTA to     msnd1 shell
    
    
    S29_111_124.w=w_cse*0        # from ACA_Pyramid_111 to msnd2 shell
    S31_151_124.w=w_tse*0.1        # from Ventral_THL_151 to     msnd2 shell
    S32_131_124.w=w_vsi*0        # from VTA to     msnd2 shell
    
    
    
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')
    
    print('Scenario 1: D2 Only cortex')    

    S01_1011_111.w=we*1
    S54_1011_131.w=we*0.1
    S92_211_131.w=we*0 
    S87_2011_231.w=we*0.1
    
    S05_111_121.w=w_cse        # from ACA_Pyramid_111 to msnd1 core
    S07_151_121.w=w_tse        # from Ventral_THL_151 to     msnd1 core
    S08_131_121.w=w_vse        # from VTA to     msnd1 core
    
    
    S13_111_122.w=w_cse        # from ACA_Pyramid_111 to msnd2 core
    S15_151_122.w=w_tse        # from Ventral_THL_151 to     msnd2 core
    S16_131_122.w=w_vsi        # from VTA to     msnd2 core
    
    
    S21_111_123.w=w_cse        # from ACA_Pyramid_111 to msnd1 shell
    S23_151_123.w=w_tse        # from Ventral_THL_151 to     msnd1 shell
    S24_131_123.w=w_vse        # from VTA to     msnd1 shell
    
    
    S29_111_124.w=w_cse        # from ACA_Pyramid_111 to msnd2 shell
    S31_151_124.w=w_tse        # from Ventral_THL_151 to     msnd2 shell
    S32_131_124.w=w_vsi        # from VTA to     msnd2 shell  
    
    
    print("sim_time="+str(duration3))
    run(duration3, report='text')    
    
    
    print('Scenario 1: D3 Resting')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')
    
    
    print('Scenario 1: D4 Only VTA')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we
    S92_211_131.w=we*0.5
    S87_2011_231.w=we*0.1
    
    S05_111_121.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd1 core
    S13_111_122.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd2 core    
    S21_111_123.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd1 shell    
    S29_111_124.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd2 shell    
    
    
    print("sim_time="+str(duration3))
    run(duration3, report='text')


    
    print('Scenario 1: D5 Resting')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')


    print('Scenario 1: D6 Cortex and VTA')    

    S01_1011_111.w=we*0.5
    S54_1011_131.w=we
    S87_2011_231.w=we*0.1
    
    
    S05_111_121.w=w_cse        # from ACA_Pyramid_111 to msnd1 core
    S13_111_122.w=w_cse        # from ACA_Pyramid_111 to msnd2 core    
    S21_111_123.w=w_cse        # from ACA_Pyramid_111 to msnd1 shell    
    S29_111_124.w=w_cse        # from ACA_Pyramid_111 to msnd2 shell    

    
    print("sim_time="+str(duration3))
    run(duration3, report='text')
    
    
    print('Scenario 1: D7 Resting')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')



##-------------------- Scenario 2 --------------------------------##

###### Changing cortex, VTA and SNc input, investigate the relation of stimulus and reward

elif number_of_scenario==2:

    print("----------------------------------------------")
    print('Scenario 2')
    
    print('Scenario 2: D1 Resting')  
    duration1=100*ms
    duration2=500*ms
    duration3=2000*ms
    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0.1
    S92_211_131.w=we*0
    S87_2011_231.w=we*0.1      
            
    
    S05_111_121.w=w_cse*0        # from ACA_Pyramid_111 to msnd1 core
    S07_151_121.w=w_tse*0.1        # from Ventral_THL_151 to     msnd1 core
    S08_131_121.w=w_vse*0        # from VTA to     msnd1 core
    
    
    S13_111_122.w=w_cse*0        # from ACA_Pyramid_111 to msnd2 core
    S15_151_122.w=w_tse*0.1        # from Ventral_THL_151 to     msnd2 core
    S16_131_122.w=w_vsi*0        # from VTA to     msnd2 core
    
    
    S21_111_123.w=w_cse*0        # from ACA_Pyramid_111 to msnd1 shell
    S23_151_123.w=w_tse*0.1        # from Ventral_THL_151 to     msnd1 shell
    S24_131_123.w=w_vse*0        # from VTA to     msnd1 shell
    
    
    S29_111_124.w=w_cse*0        # from ACA_Pyramid_111 to msnd2 shell
    S31_151_124.w=w_tse*0.1        # from Ventral_THL_151 to     msnd2 shell
    S32_131_124.w=w_vsi*0        # from VTA to     msnd2 shell
    
    
    
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')
    

    print('Scenario 2: D2 Only cortex')    

    S01_1011_111.w=we*1
    S54_1011_131.w=we*0.2
    S92_211_131.w=we*0 
    S87_2011_231.w=we*0.2
    
    S05_111_121.w=w_cse        # from ACA_Pyramid_111 to msnd1 core
    S07_151_121.w=w_tse        # from Ventral_THL_151 to     msnd1 core
    S08_131_121.w=w_vse        # from VTA to     msnd1 core
    
    
    S13_111_122.w=w_cse        # from ACA_Pyramid_111 to msnd2 core
    S15_151_122.w=w_tse        # from Ventral_THL_151 to     msnd2 core
    S16_131_122.w=w_vsi        # from VTA to     msnd2 core
    
    
    S21_111_123.w=w_cse        # from ACA_Pyramid_111 to msnd1 shell
    S23_151_123.w=w_tse        # from Ventral_THL_151 to     msnd1 shell
    S24_131_123.w=w_vse        # from VTA to     msnd1 shell
    
    
    S29_111_124.w=w_cse        # from ACA_Pyramid_111 to msnd2 shell
    S31_151_124.w=w_tse        # from Ventral_THL_151 to     msnd2 shell
    S32_131_124.w=w_vsi        # from VTA to     msnd2 shell  
    
    
    print("sim_time="+str(duration3))
    run(duration3, report='text')    
    
    
    print('Scenario 2: D3 Resting')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')
    
    
    print('Scenario 2: Only VTA')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we
    S92_211_131.w=we*0.5
    S87_2011_231.w=we*0.2
    S56_2011_211.w=we*0.5    #Poisson to PFC
    
    S05_111_121.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd1 core
    S13_111_122.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd2 core    
    S21_111_123.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd1 shell    
    S29_111_124.w=w_cse*0.1        # from ACA_Pyramid_111 to msnd2 shell    
    
    
    print("sim_time="+str(duration3))
    run(duration3, report='text')

    print('Scenario 2: D5 Resting')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    S56_2011_211.w=we*0    #Poisson to PFC
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')


    print('Scenario 2: D6 Cortex, VTA and SNc')    

    S01_1011_111.w=we*0.5
    S54_1011_131.w=we
    S87_2011_231.w=we
    
    
    S05_111_121.w=w_cse        # from ACA_Pyramid_111 to msnd1 core
    S13_111_122.w=w_cse        # from ACA_Pyramid_111 to msnd2 core    
    S21_111_123.w=w_cse        # from ACA_Pyramid_111 to msnd1 shell    
    S29_111_124.w=w_cse        # from ACA_Pyramid_111 to msnd2 shell    

    
    print("sim_time="+str(duration3))
    run(duration3, report='text')
    
    
    print('Scenario 2: D7 Resting')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')

    
else:
    print('Invalid scenario number!!!')


    
print("----------------------------------------------")
filename=str(time.time())
final_time=time.time()
print("Simulation Time:",str(final_time-init_time))

#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####


print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->  : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->  : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->  : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->  : " +str(spikes_nacc_in_125.num_spikes/number_of_neurons_in_nacc_in))
print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->  : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
print("Ventral_THL_151 : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))


print("----------------------------------------------")

print("PFC_Pyramid_211 : "+str(spikes_PFC_Pyramid_211.num_spikes)+"  ---->  : " +str(spikes_PFC_Pyramid_211.num_spikes/number_of_neurons_in_PFC_pyramid))
print("Crtx IN         : "+str(spikes_PFC_in_212.num_spikes)+"  ---->  : " +str(spikes_PFC_in_212.num_spikes/number_of_neurons_in_PFC_in))
print("MSND1 Caudate   : "+str(spikes_msnd1_caudate_221.num_spikes)+"  ---->  : " +str(spikes_msnd1_caudate_221.num_spikes/number_of_neurons_in_msnd1_caudate_221))
print("MSND2 Caudate   : "+str(spikes_msnd2_caudate_222.num_spikes)+"  ---->  : " +str(spikes_msnd2_caudate_222.num_spikes/number_of_neurons_in_msnd2_caudate_222))
print("Caudate IN      : "+str(spikes_caudate_in_223.num_spikes)+"  ---->  : " +str(spikes_caudate_in_223.num_spikes/number_of_neurons_in_caudate_in))
print("Dorsal_GPe_241  : "+str(spikes_Dorsal_GPe_241.num_spikes)+"  ---->  : " +str(spikes_Dorsal_GPe_241.num_spikes/number_of_neurons_in_GPe))
print("Dorsal_GPi_242  : "+str(spikes_Dorsal_GPi_242.num_spikes)+"  ---->  : " +str(spikes_Dorsal_GPi_242.num_spikes/number_of_neurons_in_GPe))
print("Dorsal_STN_243  : "+str(spikes_Dorsal_STN_243.num_spikes)+"  ---->  : " +str(spikes_Dorsal_STN_243.num_spikes/number_of_neurons_in_GPe))
print("SNc             : "+str(spikes_SNc_DA_231.num_spikes)+"  ---->  : " +str(spikes_SNc_DA_231.num_spikes/number_of_neurons_in_vta))
print("Dorsal_THL_251  : "+str(spikes_Dorsal_THL_251.num_spikes)+"  ---->  : " +str(spikes_Dorsal_THL_251.num_spikes/number_of_neurons_in_thl))

print("----------------------------------------------")


#########################################
############# --- Figures --- #########
#########################################

figure()
subplot(411)
plot(trace_ACA_Pyramid_111.t / ms, trace_ACA_Pyramid_111[9].v / mV)
ylabel('ACA_Pyramid_111')


subplot(412)
plot(trace_ACA_in_112.t / ms, trace_ACA_in_112[9].v / mV)
ylabel('Crtx IN')


subplot(413)
plot(trace_VTA_DA_131.t / ms, trace_VTA_DA_131[9].v / mV)
ylabel('VTA')


subplot(414)
plot(trace_Ventral_THL_151.t / ms, trace_Ventral_THL_151[9].v / mV, linewidth=1)
xlabel('time (ms)')
ylabel('Ventral_THL_151')
tight_layout()
#savefig('example1.pdf', dpi=600)

figure()
subplot(411)
plot(trace_msnd1_core_121.t / ms, trace_msnd1_core_121[9].v / mV)
ylabel('MSND1 Core')


subplot(412)
plot(trace_msnd2_core_122.t / ms, trace_msnd2_core_122[9].v / mV)
ylabel('MSND2 Core')




subplot(413)
plot(trace_msnd1_shell_123.t / ms, trace_msnd1_shell_123[9].v / mV)
ylabel('MSND1 Shell')



subplot(414)
plot(trace_msnd2_shell_124.t / ms, trace_msnd2_shell_124[9].v / mV)
xlabel('time, ms')
ylabel('MSND2 Shell')






figure()
subplot(411)
plot(trace_Ventral_GPe_141.t / ms, trace_Ventral_GPe_141[9].v / mV)
ylabel('Ventral_GPe_141')

subplot(412)
plot(trace_Ventral_GPi_142.t / ms, trace_Ventral_GPi_142[9].v / mV)
ylabel('Ventral_GPi_142')

subplot(413)
plot(trace_Ventral_STN_143.t / ms, trace_Ventral_STN_143[9].v / mV)
ylabel('Ventral_STN_143')

subplot(414)
plot(trace_nacc_in_125.t / ms, trace_nacc_in_125[9].v / mV)
xlabel('time (ms)')
ylabel('NAcc IN')

#==============================================================================
# 
#==============================================================================
#==============================================================================
#
if (number_of_scenario==1 or number_of_scenario==2):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure()
    plt.subplots_adjust(top=0.99,bottom=0.05,left=0.035,right=0.99,hspace=0.2,wspace=0.18)
    
    plt.subplot(541)
    plt.plot(spikes_ACA_Pyramid_111.t/ms, spikes_ACA_Pyramid_111.i, '.k')
    plt.axis([100,4400,0,900])
    #plt.axis('off')
    plt.ylabel('ACA');
    


    
    plt.subplot(542)
    plt.plot(spikes_VTA_DA_131.t/ms, spikes_VTA_DA_131.i, '.k')
    plt.ylabel('VTA DA');
    plt.axis([100,4400,0,100])
    
    plt.subplot(543)
    plt.plot(spikes_PFC_Pyramid_211.t/ms, spikes_PFC_Pyramid_211.i, '.k')
    plt.axis([100,4400,0,900])
    plt.ylabel('PFC');
    
    
    
    plt.subplot(544)
    plt.plot(spikes_SNc_DA_231.t/ms, spikes_SNc_DA_231.i, '.k')
    plt.ylabel('SNc DA');
    plt.axis([100,4400,0,100])
    
    
    
    plt.subplot(545)
    plt.plot(spikes_msnd1_core_121.t/ms, spikes_msnd1_core_121.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('MSND1 Core');
    
    plt.subplot(546)
    plt.plot(spikes_msnd2_core_122.t/ms, spikes_msnd2_core_122.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('MSND2 Core');
    
    
    
    plt.subplot(547)
    plt.plot(spikes_msnd1_caudate_221.t/ms, spikes_msnd1_caudate_221.i, '.k')
    plt.axis([100,4400,0,200])
    plt.ylabel('MSND1 Caudate');
    
    plt.subplot(548)
    plt.plot(spikes_msnd2_caudate_222.t/ms, spikes_msnd2_caudate_222.i, '.k')
    plt.axis([100,4400,0,200])
    plt.ylabel('MSND2 Caudate');
    
    
    
    plt.subplot(549)
    plt.plot(spikes_msnd1_shell_123.t/ms, spikes_msnd1_shell_123.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('MSND1 Shell');
    
    plt.subplot(5,4,10)
    plt.plot(spikes_msnd2_shell_124.t/ms, spikes_msnd2_shell_124.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('MSND2 Shell');
    
    plt.subplot(5,4,11)
    plt.plot(spikes_nacc_in_125.t/ms, spikes_nacc_in_125.i, '.k')
    plt.axis([100,4400,0,50])
    plt.ylabel('NAcc IN');
    
    plt.subplot(5,4,12)
    plt.plot(spikes_caudate_in_223.t/ms, spikes_caudate_in_223.i, '.k')
    plt.axis([100,4400,0,50])
    plt.ylabel('Caudate IN');
    
    
    plt.subplot(5,4,13)
    plt.plot(spikes_Ventral_STN_143.t/ms, spikes_Ventral_STN_143.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('Ventral_STN');
    
    plt.subplot(5,4,14)
    plt.plot(spikes_Ventral_THL_151.t/ms, spikes_Ventral_THL_151.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('Ventral_THL');
    
    plt.subplot(5,4,15)
    plt.plot(spikes_Dorsal_STN_243.t/ms, spikes_Dorsal_STN_243.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('Dorsal_STN');
    
    plt.subplot(5,4,16)
    plt.plot(spikes_Dorsal_THL_251.t/ms, spikes_Dorsal_THL_251.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('Dorsal_THL');
    
    plt.subplot(5,4,17)
    plt.plot(spikes_Ventral_GPi_142.t/ms, spikes_Ventral_GPi_142.i, '.k')
    plt.axis([100,4400,0,100])
    plt.ylabel('Ventral_GPi');
    plt.xlabel('zaman, ms')
    
    plt.subplot(5,4,18)
    plt.plot(spikes_Ventral_GPe_141.t/ms, spikes_Ventral_GPe_141.i, '.k')
    plt.axis([100,4400,0,100])
    plt.xlabel('zaman, ms')
    plt.ylabel('Ventral_GPe');
    
    
    plt.subplot(5,4,19)
    plt.plot(spikes_Dorsal_GPi_242.t/ms, spikes_Dorsal_GPi_242.i, '.k')
    plt.axis([100,4400,0,100])
    plt.xlabel('zaman, ms')
    plt.ylabel('Dorsal GPi');
    
    plt.subplot(5,4,20)
    plt.plot(spikes_Dorsal_GPe_241.t/ms, spikes_Dorsal_GPe_241.i, '.k')
    plt.axis([100,4400,0,100])
    plt.xlabel('zaman, ms')
    plt.ylabel('Dorsal GPe');
        
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    figurre = plt.gcf() # get current figure
    figurre.set_size_inches(18, 10)
    plt.savefig('raster_plot_sc01'+filename+'.pdf',dpi=600)
    plt.savefig('raster_plot_sc01'+filename+'.png',dpi=600)
    plt.show()
    
    

#==============================================================================
plt.figure()
plt.subplots_adjust(top=0.99,bottom=0.05,left=0.035,right=0.99,hspace=0.2,wspace=0.18)
 
plt.subplot(511)
hist_cortex = hist(spikes_ACA_Pyramid_111.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('ACA');
 
plt.subplot(512)
hist_VTA_DA_131 = hist(spikes_VTA_DA_131.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('VTA');
 
plt.subplot(513) 
hist_msnd1_core_121 = hist(spikes_msnd1_core_121.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
hist_msnd1_shell_123 = hist(spikes_msnd1_shell_123.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('MSND1');
 #
plt.subplot(514)
hist_msnd2_core_122 = hist(spikes_msnd2_core_122.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
hist_msnd2_shell_124 = hist(spikes_msnd2_shell_124.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('MSND2');
 
plt.subplot(515)
hist_Ventral_THL_151 = hist(spikes_Ventral_THL_151.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('v THL');
xlabel('time, ms')


    
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
figurre = plt.gcf() # get current figure
figurre.set_size_inches(18, 10)
plt.savefig('firing_rate_v'+filename+'.pdf',dpi=600)
plt.savefig('firing_rate_v'+filename+'.png',dpi=600)







plt.figure()
plt.subplots_adjust(top=0.99,bottom=0.05,left=0.035,right=0.99,hspace=0.2,wspace=0.18)

plt.subplot(511)
hist_cortex_pfc = hist(spikes_PFC_Pyramid_211.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('PFC');
 
plt.subplot(512)
hist_SNc_DA_231 = hist(spikes_SNc_DA_231.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('SNc');
 
plt.subplot(513)
hist_msnd1_caudate_221 = hist(spikes_msnd1_caudate_221.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('MSND1');
 #
plt.subplot(514)
hist_msnd2_caudate_222 = hist(spikes_msnd2_caudate_222.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('MSND2');
 
 
plt.subplot(515)
hist_Dorsal_THL_251 = hist(spikes_Dorsal_THL_251.t/ms, 50)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))
ylabel('d THL');
xlabel('time, ms')


    
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
figurre = plt.gcf() # get current figure
figurre.set_size_inches(18, 10)
plt.savefig('firing_rate_d'+filename+'.pdf',dpi=600)
plt.savefig('firing_rate_d'+filename+'.png',dpi=600)
plt.show()


# 
# 






#==============================================================================
#############---------LFP-------#########
#==============================================================================
#  LFPs are obtained calculating the distance of each neuron to electrodes and the total synaptic current for each neuron.        
        
### NAcc neuron dimension #################
xN_msnd1_core=5*rand(number_of_neurons_in_msnd1_core_121)
yN_msnd1_core=10*rand(number_of_neurons_in_msnd1_core_121)
zN_msnd1_core=10*rand(number_of_neurons_in_msnd1_core_121)
xN_msnd2_core=5*rand(number_of_neurons_in_msnd2_core_122)
yN_msnd2_core=10*rand(number_of_neurons_in_msnd2_core_122)
zN_msnd2_core=10*rand(number_of_neurons_in_msnd2_core_122)
xN_msnd1_shell=5+5*rand(number_of_neurons_in_msnd1_shell_123)
yN_msnd1_shell=10*rand(number_of_neurons_in_msnd1_shell_123)
zN_msnd1_shell=10*rand(number_of_neurons_in_msnd1_shell_123)
xN_msnd2_shell=5+5*rand(number_of_neurons_in_msnd2_shell_124)
yN_msnd2_shell=10*rand(number_of_neurons_in_msnd2_shell_124)
zN_msnd2_shell=10*rand(number_of_neurons_in_msnd2_shell_124)
xN_nacc_in=10*rand(number_of_neurons_in_nacc_in)
yN_nacc_in=10*rand(number_of_neurons_in_nacc_in)
zN_nacc_in=10*rand(number_of_neurons_in_nacc_in)

### Caudate neuron dimension #################
xN_msnd1_caudate=10*rand(number_of_neurons_in_msnd1_caudate_221)
yN_msnd1_caudate=10*rand(number_of_neurons_in_msnd1_caudate_221)
zN_msnd1_caudate=10*rand(number_of_neurons_in_msnd1_caudate_221)
xN_msnd2_caudate=10*rand(number_of_neurons_in_msnd2_caudate_222)
yN_msnd2_caudate=10*rand(number_of_neurons_in_msnd2_caudate_222)
zN_msnd2_caudate=10*rand(number_of_neurons_in_msnd2_caudate_222)
xN_caudate_in=10*rand(number_of_neurons_in_caudate_in)
yN_caudate_in=10*rand(number_of_neurons_in_caudate_in)
zN_caudate_in=10*rand(number_of_neurons_in_caudate_in)


array_dim=.1

xE11=5*rand()
yE11=10*rand()
zE11=10*rand()

xE21=5+5*rand()
yE21=10*rand()
zE21=10*rand()


xE3=3+4*rand()
yE3=3+4*rand()
zE3=3+4*rand()


print('Electrode 1 (Core): ('+str(xE11)+','+str(yE11)+','+str(zE11)+')')
print('Electrode 2 (Shell): ('+str(xE21)+','+str(yE21)+','+str(zE21)+')')
print('Electrode 3 (Caudate): ('+str(xE3)+','+str(yE3)+','+str(zE3)+')')


LFP_Is_vE1=zeros(len(trace_I_s_msnd1_core.Is[0]))*amp
LFP_Is_vE2=zeros(len(trace_I_s_msnd1_core.Is[0]))*amp
LFP_Is_vE3=zeros(len(trace_I_s_msnd1_caudate.Is[0]))*amp

for i_lfp in range(number_of_neurons_in_msnd1_core_121):
    d_msnd1_core_E11=((xE11-xN_msnd1_core[i_lfp])*(xE11-xN_msnd1_core[i_lfp])+(yE11-yN_msnd1_core[i_lfp])*(yE11-yN_msnd1_core[i_lfp])+(zE11-zN_msnd1_core[i_lfp])*(zE11-zN_msnd1_core[i_lfp]))
    d_msnd2_core_E11=((xE11-xN_msnd2_core[i_lfp])*(xE11-xN_msnd2_core[i_lfp])+(yE11-yN_msnd2_core[i_lfp])*(yE11-yN_msnd2_core[i_lfp])+(zE11-zN_msnd2_core[i_lfp])*(zE11-zN_msnd2_core[i_lfp]))
    d_msnd1_shell_E11=((xE11-xN_msnd1_shell[i_lfp])*(xE11-xN_msnd1_shell[i_lfp])+(yE11-yN_msnd1_shell[i_lfp])*(yE11-yN_msnd1_shell[i_lfp])+(zE11-zN_msnd1_shell[i_lfp])*(zE11-zN_msnd1_shell[i_lfp]))
    d_msnd2_shell_E11=((xE11-xN_msnd2_shell[i_lfp])*(xE11-xN_msnd2_shell[i_lfp])+(yE11-yN_msnd2_shell[i_lfp])*(yE11-yN_msnd2_shell[i_lfp])+(zE11-zN_msnd2_shell[i_lfp])*(zE11-zN_msnd2_shell[i_lfp]))
    d_msnd1_core_E21=((xE21-xN_msnd1_core[i_lfp])*(xE21-xN_msnd1_core[i_lfp])+(yE21-yN_msnd1_core[i_lfp])*(yE21-yN_msnd1_core[i_lfp])+(zE21-zN_msnd1_core[i_lfp])*(zE21-zN_msnd1_core[i_lfp]))
    d_msnd2_core_E21=((xE21-xN_msnd2_core[i_lfp])*(xE21-xN_msnd2_core[i_lfp])+(yE21-yN_msnd2_core[i_lfp])*(yE21-yN_msnd2_core[i_lfp])+(zE21-zN_msnd2_core[i_lfp])*(zE21-zN_msnd2_core[i_lfp]))
    d_msnd1_shell_E21=((xE21-xN_msnd1_shell[i_lfp])*(xE21-xN_msnd1_shell[i_lfp])+(yE21-yN_msnd1_shell[i_lfp])*(yE21-yN_msnd1_shell[i_lfp])+(zE21-zN_msnd1_shell[i_lfp])*(zE21-zN_msnd1_shell[i_lfp]))
    d_msnd2_shell_E21=((xE21-xN_msnd2_shell[i_lfp])*(xE21-xN_msnd2_shell[i_lfp])+(yE21-yN_msnd2_shell[i_lfp])*(yE21-yN_msnd2_shell[i_lfp])+(zE21-zN_msnd2_shell[i_lfp])*(zE21-zN_msnd2_shell[i_lfp]))
    if i_lfp<50:  
        d_nacc_in_E11=((xE11-xN_nacc_in[i_lfp])*(xE11-xN_nacc_in[i_lfp])+(yE11-yN_nacc_in[i_lfp])*(yE11-yN_nacc_in[i_lfp])+(zE11-zN_nacc_in[i_lfp])*(zE11-zN_nacc_in[i_lfp]))
        d_nacc_in_E21=((xE21-xN_nacc_in[i_lfp])*(xE21-xN_nacc_in[i_lfp])+(yE21-yN_nacc_in[i_lfp])*(yE21-yN_nacc_in[i_lfp])+(zE21-zN_nacc_in[i_lfp])*(zE21-zN_nacc_in[i_lfp]))
        LFP_Is_vE1=LFP_Is_vE1+trace_I_s_msnd1_core[i_lfp].Is/d_msnd1_core_E11+trace_I_s_msnd2_core[i_lfp].Is/d_msnd2_core_E11+trace_I_s_msnd1_shell[i_lfp].Is/d_msnd1_shell_E11+trace_I_s_msnd2_shell[i_lfp].Is/d_msnd2_shell_E11+trace_I_s_nacc_in[i_lfp].Is/d_nacc_in_E11
        LFP_Is_vE2=LFP_Is_vE2+trace_I_s_msnd1_core[i_lfp].Is/d_msnd1_core_E21+trace_I_s_msnd2_core[i_lfp].Is/d_msnd2_core_E21+trace_I_s_msnd1_shell[i_lfp].Is/d_msnd1_shell_E21+trace_I_s_msnd2_shell[i_lfp].Is/d_msnd2_shell_E21+trace_I_s_nacc_in[i_lfp].Is/d_nacc_in_E21
    else:         
        LFP_Is_vE1=LFP_Is_vE1+trace_I_s_msnd1_core.Is[i_lfp]/d_msnd1_core_E11+trace_I_s_msnd2_core.Is[i_lfp]/d_msnd2_core_E11+trace_I_s_msnd1_shell.Is[i_lfp]/d_msnd1_shell_E11+trace_I_s_msnd2_shell.Is[i_lfp]/d_msnd2_shell_E11
        LFP_Is_vE2=LFP_Is_vE2+trace_I_s_msnd1_core.Is[i_lfp]/d_msnd1_core_E21+trace_I_s_msnd2_core.Is[i_lfp]/d_msnd2_core_E21+trace_I_s_msnd1_shell.Is[i_lfp]/d_msnd1_shell_E21+trace_I_s_msnd2_shell.Is[i_lfp]/d_msnd2_shell_E21
    


for i_lfp in range(number_of_neurons_in_msnd1_caudate_221):
    d_msnd1_caudate_E3=((xE3-xN_msnd1_caudate[i_lfp])*(xE3-xN_msnd1_caudate[i_lfp])+(yE3-yN_msnd1_caudate[i_lfp])*(yE3-yN_msnd1_caudate[i_lfp])+(zE3-zN_msnd1_caudate[i_lfp])*(zE3-zN_msnd1_caudate[i_lfp]))
    d_msnd2_caudate_E3=((xE3-xN_msnd2_caudate[i_lfp])*(xE3-xN_msnd2_caudate[i_lfp])+(yE3-yN_msnd2_caudate[i_lfp])*(yE3-yN_msnd2_caudate[i_lfp])+(zE3-zN_msnd2_caudate[i_lfp])*(zE3-zN_msnd2_caudate[i_lfp]))
    if i_lfp<50:  
        d_caudate_in_E3=((xE3-xN_caudate_in[i_lfp])*(xE3-xN_caudate_in[i_lfp])+(yE3-yN_caudate_in[i_lfp])*(yE3-yN_caudate_in[i_lfp])+(zE3-zN_caudate_in[i_lfp])*(zE3-zN_caudate_in[i_lfp]))
        LFP_Is_vE3=LFP_Is_vE3+trace_I_s_msnd1_caudate[i_lfp].Is/d_msnd1_caudate_E3+trace_I_s_msnd2_caudate[i_lfp].Is/d_msnd2_caudate_E3+trace_I_s_caudate_in[i_lfp].Is/d_caudate_in_E3
    else:
        LFP_Is_vE3=LFP_Is_vE3+trace_I_s_msnd1_caudate[i_lfp].Is/d_msnd1_caudate_E3+trace_I_s_msnd2_caudate[i_lfp].Is/d_msnd2_caudate_E3
        
        

figure()
subplot(211)
plot(trace_I_s_msnd1_core.t / ms, trace_I_s_msnd1_core[0].Is / amp)
plot(trace_I_s_msnd1_core.t / ms, trace_I_s_msnd1_core[1].Is / amp)
plot(trace_I_s_msnd1_caudate.t / ms, trace_I_s_msnd1_caudate[1].Is / amp)
ylabel('MSN Is')


subplot(212)
plot(trace_I_s_msnd1_core.t / ms, LFP_Is_vE1 /volt)
plot(trace_I_s_msnd1_core.t / ms, LFP_Is_vE2 /volt)
plot(trace_I_s_msnd1_caudate.t / ms, LFP_Is_vE3 /volt)
ylabel('LFP_Is_E123')
xlabel('time (ms)')
tight_layout()
savefig('LFP_Is_E123'+filename+'.pdf', dpi=600)

show()

#==============================================================================
#############------LFP (end)----#########
#==============================================================================
