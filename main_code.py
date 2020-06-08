# This code is written by Rahmi Elibol at December 2019 for Ventral to Dorsal article.



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
w_cse = 0.20 *we        #for lfp old value 0.20 re 20200429
w_tse = 0.25 *we        #for lfp old value 0.25 re 20200429
w_vse = 3.00 *we
w_vsi = 3.00 *wi
w_sse = 3.00 *we
w_ssi = 3.00 *wi




par_percent=10


#Regular Spike (RS) parameters
#a = 0.02 / ms
#b = 0.2 / ms  
#c = -65 * mvolt        
#d = 8 * mvolt 

#Fast Spike (FS) parameters
#a = 0.1 / ms
#b = 0.2 / ms  
#c = -65 * mvolt        
#d = 2 * mvolt 



#Chattering (CH) parameters
#a = 0.02 / ms
#b = 0.2 / ms  
#c = -50 * mvolt        
#d = 2 * mvolt 




#number_of_neurons_in_cortex=900
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


    #nominal value: 50*Hz


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


###################################################
####### ----------- Ventral BG Groups -----#######
###################################################


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



####### ----------- Ventral THL -----#######


##      0.03      0.25    -60     4        0;...      % rebound spike
##      0.03      0.25    -52     0        0;...      % rebound burst
## THL icin rebaund burst alindi.

Ventral_THL_151 = NeuronGroup(number_of_neurons_in_thl, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)



##rebound burst
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




####### ----------- Ventral Pallidal -----#######



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



####### ----------- Dorsal THL -----#######


##      0.03      0.25    -60     4        0;...      % rebound spike
##      0.03      0.25    -52     0        0;...      % rebound burst
## THL icin rebaund burst alindi.

Dorsal_THL_251 = NeuronGroup(number_of_neurons_in_thl, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)



##rebound burst
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
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################



###### Synapses #############
#############################
print('Synapses')

################ -----  Ventral Synapses -------------- ######
# 1  korteks
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

print("----------------------------------------------")
print("NAcc Core MSND1 hucresine gelen sinaps sayilari:")


S05_111_121 = Synapses(ACA_Pyramid_111, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_glu += w')
S05_111_121.connect(True, p = 0.25)
S05_111_121.w=w_cse

print("piramitten gelen sinaps sayisi: "+str(S05_111_121.N)+" ortalama:"+str(S05_111_121.N/100))

S06_125_121 = Synapses(nacc_in_125, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_Ach += w')
S06_125_121.connect(True, p = 0.25)
S06_125_121.w=wi

print("IN gelen sinaps sayisi: "+str(S06_125_121.N)+" ortalama:"+str(S06_125_121.N/100))


S07_151_121 = Synapses( Ventral_THL_151, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_glu += w')
S07_151_121.connect(True, p = 0.25)
S07_151_121.w=w_tse

print("Ventral_THL_151 gelen sinaps sayisi: "+str(S07_151_121.N)+" ortalama:"+str(S07_151_121.N/100))

S08_131_121 = Synapses(VTA_DA_131, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_DA += w')
S08_131_121.connect(True, p = 0.25)
S08_131_121.w=w_vse

print("VTA_DA_131 gelen sinaps sayisi: "+str(S08_131_121.N)+" ortalama:"+str(S08_131_121.N/100))

##### Colateral inhibitions from MSNs

S09_121_121 = Synapses(msnd1_core_121, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S09_121_121.connect(True, p = 0.05)
S09_121_121.w=wi
print("MSND1C gelen sinaps sayisi: "+str(S09_121_121.N)+" ortalama:"+str(S09_121_121.N/100))


S10_122_121 = Synapses(msnd2_core_122, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S10_122_121.connect(True, p = 0.25)
S10_122_121.w=wi*2
print("MSND2C gelen sinaps sayisi: "+str(S10_122_121.N)+" ortalama:"+str(S10_122_121.N/100))

S11_123_121 = Synapses(msnd1_shell_123, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S11_123_121.connect(True, p = 0.05)
S11_123_121.w=wi
print("MSND1S gelen sinaps sayisi: "+str(S11_123_121.N)+" ortalama:"+str(S11_123_121.N/100))


S12_124_121 = Synapses(msnd2_shell_124, msnd1_core_121, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S12_124_121.connect(True, p = 0.25)
S12_124_121.w=wi*2

print("MSND2S gelen sinaps sayisi: "+str(S12_124_121.N)+" ortalama:"+str(S12_124_121.N/100))
print("Toplam ortalama sinaps sayilari :" +str((S05_111_121.N+S06_125_121.N+S07_151_121.N+S08_131_121.N+S09_121_121.N+S10_122_121.N+S11_123_121.N+S12_124_121.N)/100))


print("----------------------------------------------")
print("NAcc Core MSND2 hucresine gelen sinaps sayilari:")


# 122 Core MSND2

S13_111_122 = Synapses(ACA_Pyramid_111, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_glu += w')
S13_111_122.connect(True, p = 0.25)
S13_111_122.w=w_cse
print("piramitten gelen sinaps sayisi: "+str(S13_111_122.N)+" ortalama:"+str(S13_111_122.N/100))

S14_125_122 = Synapses(nacc_in_125, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S14_125_122.connect(True, p = 0.25)
S14_125_122.w=wi
print("IN gelen sinaps sayisi: "+str(S14_125_122.N)+" ortalama:"+str(S14_125_122.N/100))

S15_151_122 = Synapses( Ventral_THL_151, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_glu += w')
S15_151_122.connect(True, p = 0.25)
S15_151_122.w=w_tse
print("Ventral_THL_151 gelen sinaps sayisi: "+str(S15_151_122.N)+" ortalama:"+str(S15_151_122.N/100))

S16_131_122 = Synapses(VTA_DA_131, msnd2_core_122,  'w :siemens', delay=dly, on_pre='g_DA += w')
S16_131_122.connect(True, p = 0.25)
S16_131_122.w=w_vsi

print("VTA_DA_131 gelen sinaps sayisi: "+str(S16_131_122.N)+" ortalama:"+str(S16_131_122.N/100))


##### Colateral inhibitions from MSNs

S17_121_122 = Synapses(msnd1_core_121, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S17_121_122.connect(True, p = 0.25)
S17_121_122.w=wi*2
print("MSND1C gelen sinaps sayisi: "+str(S17_121_122.N)+" ortalama:"+str(S17_121_122.N/100))


S18_122_122 = Synapses(msnd2_core_122, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S18_122_122.connect(True, p = 0.05)
S18_122_122.w=wi
print("MSND2C gelen sinaps sayisi: "+str(S18_122_122.N)+" ortalama:"+str(S18_122_122.N/100))


S19_123_122 = Synapses(msnd1_shell_123, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S19_123_122.connect(True, p = 0.25)
S19_123_122.w=wi*2


print("MSND1S gelen sinaps sayisi: "+str(S19_123_122.N)+" ortalama:"+str(S19_123_122.N/100))


S20_124_122 = Synapses(msnd2_shell_124, msnd2_core_122, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S20_124_122.connect(True, p = 0.05)
S20_124_122.w=wi


print("MSND2S gelen sinaps sayisi: "+str(S20_124_122.N)+" ortalama:"+str(S20_124_122.N/100))
print("Toplam ortalama sinaps sayilari :" +str((S13_111_122.N+S14_125_122.N+S15_151_122.N+S16_131_122.N+S17_121_122.N+S18_122_122.N+S19_123_122.N+S20_124_122.N)/100))




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



##### Ventral_GPe_141 ve Ventral_GPi_142 baglantilari yeniden duzenlenecek...
#### Poisson dan *5 geliyordu....

 




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
S54_1011_131.w=we*0.5     #0.4 odul yok #1.0 odul var


S92_211_131 = Synapses(PFC_Pyramid_211, VTA_DA_131,  'w :siemens', delay=dly, on_pre='ge += w')
S92_211_131.connect(True, p = 0.25)
S92_211_131.w=we*0.5  





#==============================================================================
#==============================================================================
#==============================================================================
# # STDP ve Hebbian kodlari
# stdp parameters
#==============================================================================
# 
# taupre = 20*ms
# taupost = taupre
# 
# w_base=10
# gmax = 100
# dApre = .01
# dApost = -dApre * taupre / taupost * 1.05
# dApost *= gmax
# dApre *= gmax
# 
# 
# 
# 

# 
# #agirliklar run komutunun ustunde tanimli
# 

# 
# 
# dApre = .1
# #This synapse is Hebbian learning
# S111= Synapses(pyramid, pyramid, 
#              '''w : 1
#                 dApre/dt = -Apre / taupre : 1 (event-driven)''',
#              on_pre='''ge += w*Hz
#                     Apre += dApre ''',
#              on_post='''w = clip(w + Apre, 0, gmax)''', #   10 yerine ----> Apre#9.05.2018 nss ve re
#              ) 
# S111.connect()
# 
# #initial value for weights
# 
# S111.w= 5*rand()
# 
# print("ilk Hebbian agirliklar:")
# print(S111.w[300,30])
# 
# dApre = .01
# 
# 
# 
# 
# 
# 
# 
# 
# 
# S131 = Synapses(pyramid, msnd1, 
#              '''w : 1
#                 dApre/dt = -Apre / taupre : 1 (event-driven)
#                 dApost/dt = -Apost / taupost : 1 (event-driven)''',
#              on_pre='''ge += w*Hz
#                     Apre += dApre
#                     w = clip(w + Apost, 0, gmax)''',
#              on_post='''Apost += dApost
#                      w = clip(w + Apre, 0, gmax)''',
#              )
# S131.connect()
# 
# #initial value for weights
# 
# 
# for i in range(500):
#     for j in range(100):
#         S131.w[i,j] = w_base+20*rand()
# 
# #print("ilk agirliklar:")
# #print(S131.w[300,:])
# 
# 
# 
# 
# 
# 
# S141 = Synapses(pyramid, msnd2, 
#              '''w : 1
#                 dApre/dt = -Apre / taupre : 1 (event-driven)
#                 dApost/dt = -Apost / taupost : 1 (event-driven)''',
#              on_pre='''ge += w*Hz
#                     Apre += dApre
#                     w = clip(w + Apost, 0, gmax)''',
#              on_post='''Apost += dApost
#                      w = clip(w + Apre, 0, gmax)''',
#              )
# S141.connect()
# 
# #initial value for weights
# 
# 
# for i in range(500):
#     for j in range(100):
#         S141.w[i,j] = w_base+10*rand()
# 
# #print("ilk agirliklar:")
# #print(S141.w[299,:])
# 


###################################################


################ -----  Dorsal Synapses -------------- ######
# 1  korteks
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



########## Onemli Not ################      20202904    dorsal_w_carpani
#####      w_cse    Her iki grup icin   1.1 ile carpildi
#####      w_tse    Her iki grup icin   1.1 ile carpildi
#####      w_vse    (D1 ler icin)        1.1 ile carpildi
#####      w_vsi    (D2 ler icin)        1.1 ile carpildi

dorsal_w_carpani=1.1


# 2 Caudate
# 221 Caudate MSND1
print("----------------------------------------------")
print("Caudate MSND1 hucresine gelen sinaps sayilari:")

S62_211_221 = Synapses(PFC_Pyramid_211, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_glu += w')
S62_211_221.connect(True, p = 0.25)
S62_211_221.w=w_cse*dorsal_w_carpani
print("piramitten gelen sinaps sayisi: "+str(S62_211_221.N)+" ortalama:"+str(S62_211_221.N/200))

S63_223_221 = Synapses(caudate_in_223, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_Ach += w')
S63_223_221.connect(True, p = 0.25)
S63_223_221.w=wi
print("IN gelen sinaps sayisi: "+str(S63_223_221.N)+" ortalama:"+str(S63_223_221.N/200))

S64_251_221 = Synapses( Dorsal_THL_251, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_glu += w')
S64_251_221.connect(True, p = 0.25)
S64_251_221.w=w_tse*dorsal_w_carpani
print("Dorsal_THL_251 gelen sinaps sayisi: "+str(S64_251_221.N)+" ortalama:"+str(S64_251_221.N/200))

S65_231_221 = Synapses(SNc_DA_231, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_DA += w')
S65_231_221.connect(True, p = 0.25)
S65_231_221.w=w_vse*dorsal_w_carpani
print("VTA gelen sinaps sayisi: "+str(S65_231_221.N)+" ortalama:"+str(S65_231_221.N/200))

##### Colateral inhibitions from MSNs

S66_221_221 = Synapses(msnd1_caudate_221, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S66_221_221.connect(True, p = 0.05)
S66_221_221.w=wi
print("MSND1C gelen sinaps sayisi: "+str(S66_221_221.N)+" ortalama:"+str(S66_221_221.N/200))

S67_222_221 = Synapses(msnd2_caudate_222, msnd1_caudate_221, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S67_222_221.connect(True, p = 0.25)
S67_222_221.w=wi*2
print("MSND2C gelen sinaps sayisi: "+str(S67_222_221.N)+" ortalama:"+str(S67_222_221.N/200))

print("Toplam ortalama sinaps sayilari :" +str((S62_211_221.N+S63_223_221.N+S64_251_221.N+S65_231_221.N+S66_221_221.N+S67_222_221.N)/200))





# 222 Caudate MSND2

S68_211_222 = Synapses(PFC_Pyramid_211, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_glu += w')
S68_211_222.connect(True, p = 0.25)
S68_211_222.w=w_cse*dorsal_w_carpani

S69_223_222 = Synapses(caudate_in_223, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S69_223_222.connect(True, p = 0.25)
S69_223_222.w=wi

S70_251_222 = Synapses(Dorsal_THL_251, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_glu += w')
S70_251_222.connect(True, p = 0.25)
S70_251_222.w=w_tse*dorsal_w_carpani

S71_231_222 = Synapses(SNc_DA_231, msnd2_caudate_222,  'w :siemens', delay=dly, on_pre='g_DA += w')
S71_231_222.connect(True, p = 0.25)
S71_231_222.w=w_vsi*dorsal_w_carpani


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




##### Dorsal_GPe_241 ve Dorsal_GPi_242 baglantilari yeniden duzenlenecek...
#### Poisson dan *5 geliyordu....

 




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
S86_242_251.connect(True, p = 0.5)





# 231 SNc

S87_2011_231 = Synapses(PG_PFC_pyramid_2011, SNc_DA_231,  'w :siemens', delay=dly, on_pre='ge += w')
S87_2011_231.connect(True, p = 0.25)
S87_2011_231.w=we     #0.4 odul yok #1.0 odul var


###############################################################################
###############################################################################
############## ---------------- Motor Cortex -----------------------  #########
###############################################################################
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



###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================









import time
init_time=time.time()


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
#monge_ACA_Pyramid_111 = StateMonitor(ACA_Pyramid_111, 'ge', record=True)
#mongi_ACA_Pyramid_111 = StateMonitor(ACA_Pyramid_111, 'gi', record=True)
spikes_ACA_Pyramid_111 = SpikeMonitor(ACA_Pyramid_111)


trace_ACA_in_112 = StateMonitor(ACA_in_112, 'v', record=9)
#monge_ACA_in_112 = StateMonitor(ACA_in_112, 'ge', record=True)
#mongi_ACA_in_112 = StateMonitor(ACA_in_112, 'gi', record=True)
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
#monge_ACA_in_112 = StateMonitor(ACA_in_112, 'ge', record=True)
#mongi_ACA_in_112 = StateMonitor(ACA_in_112, 'gi', record=True)
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
########### Excitation currents ##########
###############################################################################
print("----------------------------------------------")
print("start")

number_of_scenario=2


#########################  ----- Scenarios   --------------  ####################

##-------------------- Scenario 0 --------------------------------##
#### 10 ms sureli kodun calisip calismadigini teset etmek icin kullanilir. 


if number_of_scenario==0:

    
    print("----------------------------------------------")
    print('Senaryo 0')
    run(250*ms,report='text')



##-------------------- Scenario 1 --------------------------------##

##### Korteks ve VTA girisleri degistirilerek odul ve uyaran iliskisi icin kullanilir. 

elif number_of_scenario==1:

    print("----------------------------------------------")
    print('Senaryo 1')
    
    print('Senaryo 1: D1 Dinlenim')  
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
    
    print("D1 Dinlenimden sonra")
    print("ACA_Pyramid_111 :"+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->   : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->       : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->   : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->   : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->      : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->       : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL     : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))
    
    
    print('Senaryo 1: D2 Korteks var')    

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
    
    print("D2 Korteksten sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->       : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->   : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->   : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->      : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->       : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL     : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))    
    
    print('Senaryo 1: D3 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')
    
    
    print('Senaryo 1: D4 VTA var')    
    
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

    
    print("D4 VTA dan sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN : "+str(spikes_ACA_in_112.num_spikes)+"  ---->  : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->  : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->  : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN : "+str(spikes_nacc_in_125.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->  : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))    
    
    
    print('Senaryo 1: D5 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')


    print('Senaryo 1: D6 Korteks ve VTA var')    

    S01_1011_111.w=we*0.5
    S54_1011_131.w=we
    S87_2011_231.w=we*0.1
    
    
    S05_111_121.w=w_cse        # from ACA_Pyramid_111 to msnd1 core
    S13_111_122.w=w_cse        # from ACA_Pyramid_111 to msnd2 core    
    S21_111_123.w=w_cse        # from ACA_Pyramid_111 to msnd1 shell    
    S29_111_124.w=w_cse        # from ACA_Pyramid_111 to msnd2 shell    

    
    print("sim_time="+str(duration3))
    run(duration3, report='text')
    
    
    
    print("D6 Korteks + VTA dan sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN : "+str(spikes_ACA_in_112.num_spikes)+"  ---->  : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->  : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->  : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN : "+str(spikes_nacc_in_125.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->  : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))
    
    print('Senaryo 1: D7 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')



##-------------------- Scenario 2 --------------------------------##

##### Korteks, VTA ve SNc girisleri degistirilerek odul ve uyaran iliskisi icin kullanilir. 

elif number_of_scenario==2:

    print("----------------------------------------------")
    print('Senaryo 2')
    
    print('Senaryo 2: D1 Dinlenim')  
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
    
    print("D1 Dinlenimden sonra")
    print("ACA_Pyramid_111 :"+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->   : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->       : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->   : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->   : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->      : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->       : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL     : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))
    
    
    print('Senaryo 2: D2 Korteks var')    

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
    
    print("D2 Korteksten sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->       : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->   : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->   : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->      : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->       : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL     : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))    
    
    print('Senaryo 2: D3 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')
    
    
    print('Senaryo 2: D4 VTA  var (PFC)')    
    
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

    
    print("D4 VTA dan sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN : "+str(spikes_ACA_in_112.num_spikes)+"  ---->  : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->  : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->  : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN : "+str(spikes_nacc_in_125.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->  : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))    
    
    
    print('Senaryo 2: D5 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    S56_2011_211.w=we*0    #Poisson to PFC
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')


    print('Senaryo 2: D6 Korteks, VTA ve SNc var')    

    S01_1011_111.w=we*0.5
    S54_1011_131.w=we
    S87_2011_231.w=we
    
    
    S05_111_121.w=w_cse        # from ACA_Pyramid_111 to msnd1 core
    S13_111_122.w=w_cse        # from ACA_Pyramid_111 to msnd2 core    
    S21_111_123.w=w_cse        # from ACA_Pyramid_111 to msnd1 shell    
    S29_111_124.w=w_cse        # from ACA_Pyramid_111 to msnd2 shell    

    
    print("sim_time="+str(duration3))
    run(duration3, report='text')
    
    
    
    print("D6 Korteks + VTA + SNc den sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN : "+str(spikes_ACA_in_112.num_spikes)+"  ---->  : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->  : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->  : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN : "+str(spikes_nacc_in_125.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->  : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))
    
    print('Senaryo 2: D7 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')




##-------------------- Scenario 3 --------------------------------##

##### Korteks var D2 lerin aktivitesi icin kullanilir



elif number_of_scenario==3:

    print("----------------------------------------------")
    print('Senaryo 3')
    
    print('Senaryo 3: D1 Dinlenim')  
    duration1=100*ms
    duration2=500*ms
    duration3=2000*ms
    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0.1
    S92_211_131.w=we*0
    S87_2011_231.w=we*0.1
    S56_2011_211.w=we*0      
            
    
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
    
    print("D1 Dinlenimden sonra")
    print("ACA_Pyramid_111 :"+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->   : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->       : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->   : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->   : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->      : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->       : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL     : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))
    
    
    print('Senaryo 3: D2 Korteks var')    

    S01_1011_111.w=we*1
    S54_1011_131.w=we*0.1
    S92_211_131.w=we*0 
    S87_2011_231.w=we*0.1
    S56_2011_211.w=we*1

    
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
    
    print("D2 Korteksten sonra")
    print("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid))
    print("Crtx IN         : "+str(spikes_ACA_in_112.num_spikes)+"  ---->       : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in))
    print("MSND1 Core      : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->   : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121))
    print("MSND2 Core      : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->   : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122))
    print("MSND1 Shell     : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123))
    print("MSND2 Shell     : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124))
    print("NAcc IN         : "+str(spikes_nacc_in_125.num_spikes)+"  ---->      : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in))
    print("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe))
    print("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe))
    print("VTA             : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->       : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta))
    print("Ventral THL     : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl))    
    
    print('Senaryo 3: D3 Dinlenim')    
    
    S01_1011_111.w=we*0
    S54_1011_131.w=we*0
    S87_2011_231.w=we*0
    S56_2011_211.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')



#################################   scenario 3 end ####################################
#######################################################################################






    
else:
    print('Senaryo falan yok!!!!')


    
print("----------------------------------------------")
filename=str(time.time())
final_time=time.time()
print("benzetim suresi:",str(final_time-init_time))

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
#############synchronization_measure
#########################################
#
#
#v_values=matrix(trace.v/mV)
#(rows_v,columns_v)=shape(v_values)
#every_t_time=range(0,columns_v)
#
#
#
#temp_time = np.zeros(columns_v)
#for i in range(columns_v):
#	for j in range(rows_v):
#		temp_time[i]+=v_values[j,i]
#average_time=temp_time/rows_v
#print(average_time)
#
#general_variation=np.var(average_time)
#print ('general_variation : ' + str(general_variation))
#
#
#individual_variation = np.zeros(rows_v)
#for i in range(rows_v):
#	individual_variation[i]=np.var(v_values[i])
#	print ('individual_variation['  + str(i) + '] :' + str(individual_variation[i]))
#
#
#print ('average(individual_variation): ' + str(average(individual_variation)))
#sync_measure=(general_variation/average(individual_variation))
#print ('sync_measure : ' + str(sync_measure))
#
#
#########################################
#############synchronization_measure
#########################################


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
# ###################Is=ge*(Ve-v)+gi*(Vi-v) :  volt/second
# figure()
# subplot(511)
# for  i in range(1):
# 	plot(monge_ACA_Pyramid_111.t / ms, monge_ACA_Pyramid_111[i].ge*(Ve-trace_ACA_Pyramid_111[i].v)+mongi_ACA_Pyramid_111[i].gi*(Vi-trace_ACA_Pyramid_111[i].v))
# #xlabel('time (ms)')
# #ylabel('ACA_Pyramid_111')
# 
# #
# subplot(512)
# 
# for  i in range(1):
#     plot(monge_ACA_in_112.t / ms, monge_ACA_in_112[i].ge*(Ve-trace_ACA_in_112[i].v)+mongi_ACA_in_112[i].gi*(Vi-trace_ACA_in_112[i].v))
# #xlabel('time (ms)')
# #ylabel('crtx in')
# #
# #
# #
# subplot(513)
# for  i in range(1):
#     plot(monge_msnd1_core_121.t / ms, monge_msnd1_core_121[i].ge*(Ve-trace_msnd1_core_121[i].v)+mongi_msnd1_core_121[i].gi*(Vi-trace_msnd1_core_121[i].v))
# 
# #xlabel('time (ms)')
# #ylabel('msnd1_core_121')
# #
# #
# #
# subplot(514)
# for  i in range(1):
#     plot(monge_msnd2_core_122.t / ms, monge_msnd2_core_122[i].ge*(Ve-trace_msnd2_core_122[i].v)+mongi_msnd2_core_122[i].gi*(Vi-trace_msnd2_core_122[i].v))
# #xlabel('time (ms)')
# #ylabel('msnd2_core_122')
# #
# #
# #
# #
# subplot(515)
# for  i in range(1):
#     plot(monge_nacc_in_125.t / ms, monge_nacc_in_125[i].ge*(Ve-trace_nacc_in_125[i].v)+mongi_nacc_in_125[i].gi*(Vi-trace_nacc_in_125[i].v))
# 
# #xlabel('time (ms)')
# #ylabel('str in')
# 
# 
#==============================================================================
#==============================================================================
# 
# 
# figure()
# subplot(411)
# plot(spikes_ACA_Pyramid_111.t/ms, spikes_ACA_Pyramid_111.i, '.k')
# 
# ylabel('ACA_Pyramid_111');
# 
# subplot(412)
# plot(spikes_ACA_in_112.t/ms, spikes_ACA_in_112.i, '.k')
# 
# ylabel('Crtx IN');
# 
# 
# subplot(413)
# plot(spikes_VTA_DA_131.t/ms, spikes_VTA_DA_131.i, '.k')
# ylabel('VTA');
# 
# 
# subplot(414)
# plot(spikes_Ventral_THL_151.t/ms, spikes_Ventral_THL_151.i, '.k')
# xlabel('Time (ms)')
# ylabel('Ventral_THL_151');
# 
# 
# 
# figure()
# 
# subplot(511)
# plot(spikes_msnd1_core_121.t/ms, spikes_msnd1_core_121.i, '.k')
# ylabel('MSND1 Core');
# 
# subplot(512)
# plot(spikes_msnd2_core_122.t/ms, spikes_msnd2_core_122.i, '.k')
# ylabel('MSND2 Core');
# 
# 
# subplot(513)
# plot(spikes_msnd1_shell_123.t/ms, spikes_msnd1_shell_123.i, '.k')
# ylabel('MSND1 Shell');
# 
# subplot(514)
# plot(spikes_msnd2_shell_124.t/ms, spikes_msnd2_shell_124.i, '.k')
# ylabel('MSND2 Shell');
# 
# 
# subplot(515)
# plot(spikes_nacc_in_125.t/ms, spikes_nacc_in_125.i, '.k')
# xlabel('time, ms')
# ylabel('NAcc IN');
# 
# 
#==============================================================================

#==============================================================================
# figure()
# subplot(311)
# plot(spikes_Ventral_GPe_141.t/ms, spikes_Ventral_GPe_141.i, '.k')
# ylabel('Ventral_GPe_141');
# 
# subplot(312)
# plot(spikes_Ventral_GPi_142.t/ms, spikes_Ventral_GPi_142.i, '.k')
# ylabel('Ventral_GPi_142');
# 
# 
# subplot(313)
# plot(spikes_Ventral_STN_143.t/ms, spikes_Ventral_STN_143.i, '.k')
# ylabel('Ventral_STN_143');
# 
#==============================================================================

################################################################
################################################################
#
if (number_of_scenario==1 or number_of_scenario==2 or number_of_scenario==3):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure()
    plt.subplots_adjust(top=0.99,bottom=0.05,left=0.035,right=0.99,hspace=0.2,wspace=0.18)
    
    plt.subplot(541)
    plt.plot(spikes_ACA_Pyramid_111.t/ms, spikes_ACA_Pyramid_111.i, '.k')
    plt.axis([100,7000,0,900])
    #plt.axis('off')
    plt.ylabel('ACA');
    


    
    plt.subplot(542)
    plt.plot(spikes_VTA_DA_131.t/ms, spikes_VTA_DA_131.i, '.k')
    plt.ylabel('VTA DA');
    plt.axis([100,7000,0,100])
    
    plt.subplot(543)
    plt.plot(spikes_PFC_Pyramid_211.t/ms, spikes_PFC_Pyramid_211.i, '.k')
    plt.axis([100,7000,0,900])
    plt.ylabel('PFC');
    
    
    
    plt.subplot(544)
    plt.plot(spikes_SNc_DA_231.t/ms, spikes_SNc_DA_231.i, '.k')
    plt.ylabel('SNc DA');
    plt.axis([100,7000,0,100])
    
    
    
    plt.subplot(545)
    plt.plot(spikes_msnd1_core_121.t/ms, spikes_msnd1_core_121.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('MSND1 Core');
    
    plt.subplot(546)
    plt.plot(spikes_msnd2_core_122.t/ms, spikes_msnd2_core_122.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('MSND2 Core');
    
    
    
    plt.subplot(547)
    plt.plot(spikes_msnd1_caudate_221.t/ms, spikes_msnd1_caudate_221.i, '.k')
    plt.axis([100,7000,0,200])
    plt.ylabel('MSND1 Caudate');
    
    plt.subplot(548)
    plt.plot(spikes_msnd2_caudate_222.t/ms, spikes_msnd2_caudate_222.i, '.k')
    plt.axis([100,7000,0,200])
    plt.ylabel('MSND2 Caudate');
    
    
    
    plt.subplot(549)
    plt.plot(spikes_msnd1_shell_123.t/ms, spikes_msnd1_shell_123.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('MSND1 Shell');
    
    plt.subplot(5,4,10)
    plt.plot(spikes_msnd2_shell_124.t/ms, spikes_msnd2_shell_124.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('MSND2 Shell');
    
    plt.subplot(5,4,11)
    plt.plot(spikes_nacc_in_125.t/ms, spikes_nacc_in_125.i, '.k')
    plt.axis([100,7000,0,50])
    plt.ylabel('NAcc IN');
    
    plt.subplot(5,4,12)
    plt.plot(spikes_caudate_in_223.t/ms, spikes_caudate_in_223.i, '.k')
    plt.axis([100,7000,0,50])
    plt.ylabel('Caudate IN');
    
    
    plt.subplot(5,4,13)
    plt.plot(spikes_Ventral_STN_143.t/ms, spikes_Ventral_STN_143.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('Ventral_STN');
    
    plt.subplot(5,4,14)
    plt.plot(spikes_Ventral_THL_151.t/ms, spikes_Ventral_THL_151.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('Ventral_THL');
    
    plt.subplot(5,4,15)
    plt.plot(spikes_Dorsal_STN_243.t/ms, spikes_Dorsal_STN_243.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('Dorsal_STN');
    
    plt.subplot(5,4,16)
    plt.plot(spikes_Dorsal_THL_251.t/ms, spikes_Dorsal_THL_251.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('Dorsal_THL');
    
    plt.subplot(5,4,17)
    plt.plot(spikes_Ventral_GPi_142.t/ms, spikes_Ventral_GPi_142.i, '.k')
    plt.axis([100,7000,0,100])
    plt.ylabel('Ventral_GPi');
    plt.xlabel('zaman, ms')
    
    plt.subplot(5,4,18)
    plt.plot(spikes_Ventral_GPe_141.t/ms, spikes_Ventral_GPe_141.i, '.k')
    plt.axis([100,7000,0,100])
    plt.xlabel('zaman, ms')
    plt.ylabel('Ventral_GPe');
    
    
    plt.subplot(5,4,19)
    plt.plot(spikes_Dorsal_GPi_242.t/ms, spikes_Dorsal_GPi_242.i, '.k')
    plt.axis([100,7000,0,100])
    plt.xlabel('zaman, ms')
    plt.ylabel('Dorsal GPi');
    
    plt.subplot(5,4,20)
    plt.plot(spikes_Dorsal_GPe_241.t/ms, spikes_Dorsal_GPe_241.i, '.k')
    plt.axis([100,7000,0,100])
    plt.xlabel('zaman, ms')
    plt.ylabel('Dorsal GPe');
        
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    figurre = plt.gcf() # get current figure
    figurre.set_size_inches(18, 10)
    plt.savefig('raster_plot_sc01'+filename+'.pdf',dpi=100)
    plt.savefig('raster_plot_sc01'+filename+'.png',dpi=150)
    plt.show()
    
    
    
    

#hist_cortex = hist(spikes_ACA_Pyramid_111.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_nacc_in_125))/(number_of_neurons_in_nacc_in_125*defaultclock.dt))




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
# #==============================================================================
# #========     Firing Rate data ====================================================
# #==============================================================================
# 
# 
# fr_cortex=hist_cortex[0]
# 
# 
# fr_msnd1c=hist_msnd1_core_121[0]
# fr_msnd2c=hist_msnd2_core_122[0]
# fr_msnd1s=hist_msnd1_shell_123[0]
# fr_msnd2s=hist_msnd2_shell_124[0]
# 
# 
# fr_VTA_DA_131=hist_VTA_DA_131[0]
# 
# fr_Ventral_THL_151=hist_Ventral_THL_151[0]
# #==============================================================================
# # fr_nacc_in=hist_nacc_in_125[0]
# #==============================================================================
# 
# 
# fr_nacc=fr_msnd1c+fr_msnd2c+fr_msnd1s+fr_msnd2s
# fr_D1=fr_msnd1c+fr_msnd1s
# fr_D2=fr_msnd2c+fr_msnd2s
# 
# #==============================================================================
# #========     File Writing ====================================================
# #==============================================================================
# 

 
 
filename_fundamentals_informations="001_filename_fundamentals_informations_"+filename+".dat"
with open(filename_fundamentals_informations,"w") as results_filename_fundamentals_informations:
    results_filename_fundamentals_informations.write("benzetim suresi : "+str(final_time-init_time)+"\n")
    results_filename_fundamentals_informations.write("ACA_Pyramid_111 : "+str(spikes_ACA_Pyramid_111.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_ACA_pyramid)+"\n")
    results_filename_fundamentals_informations.write("Crtx IN : "+str(spikes_ACA_in_112.num_spikes)+"  ---->  : " +str(spikes_ACA_in_112.num_spikes/number_of_neurons_in_ACA_in)+"\n")
    results_filename_fundamentals_informations.write("MSND1 Core : "+str(spikes_msnd1_core_121.num_spikes)+"  ---->  : " +str(spikes_msnd1_core_121.num_spikes/number_of_neurons_in_msnd1_core_121)+"\n")
    results_filename_fundamentals_informations.write("MSND2 Core : "+str(spikes_msnd2_core_122.num_spikes)+"  ---->  : " +str(spikes_msnd2_core_122.num_spikes/number_of_neurons_in_msnd2_core_122)+"\n")
    results_filename_fundamentals_informations.write("MSND1 Shell : "+str(spikes_msnd1_shell_123.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell_123.num_spikes/number_of_neurons_in_msnd1_shell_123)+"\n")
    results_filename_fundamentals_informations.write("MSND2 Shell : "+str(spikes_msnd2_shell_124.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell_124.num_spikes/number_of_neurons_in_msnd2_shell_124)+"\n")
    results_filename_fundamentals_informations.write("NAcc IN : "+str(spikes_nacc_in_125.num_spikes)+"  ---->  : " +str(spikes_ACA_Pyramid_111.num_spikes/number_of_neurons_in_nacc_in)+"\n")
    results_filename_fundamentals_informations.write("Ventral_GPe_141 : "+str(spikes_Ventral_GPe_141.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPe_141.num_spikes/number_of_neurons_in_GPe)+"\n")
    results_filename_fundamentals_informations.write("Ventral_GPi_142 : "+str(spikes_Ventral_GPi_142.num_spikes)+"  ---->  : " +str(spikes_Ventral_GPi_142.num_spikes/number_of_neurons_in_GPe)+"\n")
    results_filename_fundamentals_informations.write("Ventral_STN_143 : "+str(spikes_Ventral_STN_143.num_spikes)+"  ---->  : " +str(spikes_Ventral_STN_143.num_spikes/number_of_neurons_in_GPe)+"\n")
    results_filename_fundamentals_informations.write("VTA : "+str(spikes_VTA_DA_131.num_spikes)+"  ---->  : " +str(spikes_VTA_DA_131.num_spikes/number_of_neurons_in_vta)+"\n")
    results_filename_fundamentals_informations.write("Ventral_THL_151 : "+str(spikes_Ventral_THL_151.num_spikes)+"  ---->  : " +str(spikes_Ventral_THL_151.num_spikes/number_of_neurons_in_thl)+"\n")
 
 
 
 
# =============================================================================
# filename_crtx="011_crtx_"+filename+".dat"
# with open(filename_crtx,"w") as results_crtx:
# for i in range(len(fr_cortex)):
#  		 	results_crtx.write(str(fr_cortex[i])+"\n")
# 
# 
# filename_str="012_str_"+filename+".dat"
# with open(filename_str,"w") as results_str:
# for i in range(len(fr_nacc)):
#  		 	results_str.write(str(fr_nacc[i])+"\n")
# =============================================================================
 
 
# =============================================================================
# 
# filename_msnd1c="001_trace_msnd1_core_121_"+filename+".dat"
# with open(filename_msnd1c,"w") as results_msnd1c:
# for i in range(len(trace_msnd1_core_121.t)):
#  		 	results_msnd1c.write(str(trace_msnd1_core_121.t[i]/ms)+","+str(trace_msnd1_core_121[9].v[i]/mV)+"\n")
# 
# 
# filename_msnd2c="002_trace_msnd2_core_122_"+filename+".dat"
# with open(filename_msnd2c,"w") as results_msnd2c:
# for i in range(len(trace_msnd2_core_122.t)):
#  		 	results_msnd2c.write(str(trace_msnd2_core_122.t[i]/ms)+","+str(trace_msnd2_core_122[9].v[i]/mV)+"\n")
# 
#  
# filename_msnd1s="003_trace_msnd1_shell_123_"+filename+".dat"
# with open(filename_msnd1s,"w") as results_msnd1s:
# for i in range(len(trace_msnd1_shell_123.t)):
#  		 	results_msnd1s.write(str(trace_msnd1_shell_123.t[i]/ms)+","+str(trace_msnd1_shell_123[9].v[i]/mV)+"\n")
#           
# filename_msnd2s="004_trace_msnd2_shell_124_"+filename+".dat"
# with open(filename_msnd2s,"w") as results_msnd2s:
# for i in range(len(trace_msnd2_shell_124.t)):
#  		 	results_msnd2s.write(str(trace_msnd2_shell_124.t[i]/ms)+","+str(trace_msnd2_shell_124[9].v[i]/mV)+"\n")              
#  
#  
# 
# filename_innacc="005_trace_in_nacc_"+filename+".dat"
# with open(filename_innacc,"w") as results_innacc:
# for i in range(len(trace_nacc_in_125.t)):
#  		 	results_innacc.write(str(trace_nacc_in_125.t[i]/ms)+","+str(trace_nacc_in_125[9].v[i]/mV)+"\n")
# =============================================================================
print(" - - -     ventral file print     - - - ")


filename_rasterplot_aca="00111_rasterplot_aca_111_"+filename+".dat"
with open(filename_rasterplot_aca,"w") as filename_rasterplot_aca:
    for i_c in range(len(spikes_ACA_Pyramid_111.t/ms)):
        filename_rasterplot_aca.write(str(spikes_ACA_Pyramid_111.t[i_c]/ms)+","+str(spikes_ACA_Pyramid_111.i[i_c])+"\n")

filename_rasterplot_aca_in="00112_rasterplot_in_112_"+filename+".dat"
with open(filename_rasterplot_aca_in,"w") as filename_rasterplot_aca_in:
    for i_c in range(len(spikes_ACA_in_112.t/ms)):
        filename_rasterplot_aca_in.write(str(spikes_ACA_in_112.t[i_c]/ms)+","+str(spikes_ACA_in_112.i[i_c])+"\n")

filename_rasterplot_VTA_DA="00141_rasterplot_VTA_DA_141_"+filename+".dat"
with open(filename_rasterplot_VTA_DA,"w") as filename_rasterplot_VTA_DA:
    for i_c in range(len(spikes_VTA_DA_131.t/ms)):
        filename_rasterplot_VTA_DA.write(str(spikes_VTA_DA_131.t[i_c]/ms)+","+str(spikes_VTA_DA_131.i[i_c])+"\n")


filename_rasterplot_msnd1c="00121_rasterplot_msnd1_core_121_"+filename+".dat"
with open(filename_rasterplot_msnd1c,"w") as results_rasterplot__msnd1c:
    for i_c in range(len(spikes_msnd1_core_121.t/ms)):
        results_rasterplot__msnd1c.write(str(spikes_msnd1_core_121.t[i_c]/ms)+","+str(spikes_msnd1_core_121.i[i_c])+"\n")

    
filename_rasterplot_msnd2c="00122_rasterplot_msnd2_core_122_"+filename+".dat"
with open(filename_rasterplot_msnd2c,"w") as results_rasterplot__msnd2c:
    for i_c in range(len(spikes_msnd2_core_122.t/ms)):
        results_rasterplot__msnd2c.write(str(spikes_msnd2_core_122.t[i_c]/ms)+","+str(spikes_msnd2_core_122.i[i_c])+"\n")

    
filename_rasterplot_msnd1s="00123_rasterplot_msnd1_shell_123_"+filename+".dat"
with open(filename_rasterplot_msnd1s,"w") as results_rasterplot__msnd1s:
    for i_c in range(len(spikes_msnd1_shell_123.t/ms)):
        results_rasterplot__msnd1s.write(str(spikes_msnd1_shell_123.t[i_c]/ms)+","+str(spikes_msnd1_shell_123.i[i_c])+"\n")

    
filename_rasterplot_msnd2s="00124_rasterplot_msnd2_shell_124_"+filename+".dat"
with open(filename_rasterplot_msnd2s,"w") as results_rasterplot__msnd2s:
    for i_c in range(len(spikes_msnd2_shell_124.t/ms)):
        results_rasterplot__msnd2s.write(str(spikes_msnd2_shell_124.t[i_c]/ms)+","+str(spikes_msnd2_shell_124.i[i_c])+"\n")        

filename_rasterplot_nacc_in="00125_rasterplot_nacc_in_125_"+filename+".dat"
with open(filename_rasterplot_nacc_in,"w") as filename_rasterplot_nacc_in:
    for i_c in range(len(spikes_nacc_in_125.t/ms)):
        filename_rasterplot_nacc_in.write(str(spikes_nacc_in_125.t[i_c]/ms)+","+str(spikes_nacc_in_125.i[i_c])+"\n")        



filename_rasterplot_vGPe="00141_rasterplot_vGPe_141_"+filename+".dat"
with open(filename_rasterplot_vGPe,"w") as filename_rasterplot_vGPe:
    for i_c in range(len(spikes_Ventral_GPe_141.t/ms)):
        filename_rasterplot_vGPe.write(str(spikes_Ventral_GPe_141.t[i_c]/ms)+","+str(spikes_Ventral_GPe_141.i[i_c])+"\n") 


filename_rasterplot_vGPi="00142_rasterplot_vGPi_142_"+filename+".dat"
with open(filename_rasterplot_vGPi,"w") as filename_rasterplot_vGPi:
    for i_c in range(len(spikes_Ventral_GPi_142.t/ms)):
        filename_rasterplot_vGPi.write(str(spikes_Ventral_GPi_142.t[i_c]/ms)+","+str(spikes_Ventral_GPi_142.i[i_c])+"\n") 

filename_rasterplot_vSTN="00143_rasterplot_vSTN_143_"+filename+".dat"
with open(filename_rasterplot_vSTN,"w") as filename_rasterplot_vSTN:
    for i_c in range(len(spikes_Ventral_STN_143.t/ms)):
        filename_rasterplot_vSTN.write(str(spikes_Ventral_STN_143.t[i_c]/ms)+","+str(spikes_Ventral_STN_143.i[i_c])+"\n") 

filename_rasterplot_vTHL="00151_rasterplot_vTHL_151_"+filename+".dat"
with open(filename_rasterplot_vTHL,"w") as filename_rasterplot_vTHL:
    for i_c in range(len(spikes_Ventral_THL_151.t/ms)):
        filename_rasterplot_vTHL.write(str(spikes_Ventral_THL_151.t[i_c]/ms)+","+str(spikes_Ventral_THL_151.i[i_c])+"\n")
# 
# ##############################################################################
# ##############################################################################
# ##############################################################################
print(" - - -     dorsal file print     - - - ")



filename_rasterplot_pfc="00211_rasterplot_pfc_211_"+filename+".dat"
with open(filename_rasterplot_pfc,"w") as filename_rasterplot_pfc:
    for i_c in range(len(spikes_PFC_Pyramid_211.t/ms)):
        filename_rasterplot_pfc.write(str(spikes_PFC_Pyramid_211.t[i_c]/ms)+","+str(spikes_PFC_Pyramid_211.i[i_c])+"\n")

filename_rasterplot_pfc_in="00212_rasterplot_in_212_"+filename+".dat"
with open(filename_rasterplot_pfc_in,"w") as filename_rasterplot_pfc_in:
    for i_c in range(len(spikes_PFC_in_212.t/ms)):
        filename_rasterplot_pfc_in.write(str(spikes_PFC_in_212.t[i_c]/ms)+","+str(spikes_PFC_in_212.i[i_c])+"\n")

filename_rasterplot_SNc_DA="00241_rasterplot_SNc_DA_241_"+filename+".dat"
with open(filename_rasterplot_SNc_DA,"w") as filename_rasterplot_SNc_DA:
    for i_c in range(len(spikes_SNc_DA_231.t/ms)):
        filename_rasterplot_SNc_DA.write(str(spikes_SNc_DA_231.t[i_c]/ms)+","+str(spikes_SNc_DA_231.i[i_c])+"\n")


filename_rasterplot_msnd1caudate="00221_rasterplot_msnd1_caudate_221_"+filename+".dat"
with open(filename_rasterplot_msnd1caudate,"w") as filename_rasterplot_msnd1caudate:
    for i_c in range(len(spikes_msnd1_caudate_221.t/ms)):
        filename_rasterplot_msnd1caudate.write(str(spikes_msnd1_caudate_221.t[i_c]/ms)+","+str(spikes_msnd1_caudate_221.i[i_c])+"\n")

    
filename_rasterplot_msnd2caudate="00222_rasterplot_msnd2_caudate_222_"+filename+".dat"
with open(filename_rasterplot_msnd2caudate,"w") as filename_rasterplot_msnd2caudate:
    for i_c in range(len(spikes_msnd2_caudate_222.t/ms)):
        filename_rasterplot_msnd2caudate.write(str(spikes_msnd2_caudate_222.t[i_c]/ms)+","+str(spikes_msnd2_caudate_222.i[i_c])+"\n")

    
filename_rasterplot_caudate_in="00223_rasterplot_caudate_in_223_"+filename+".dat"
with open(filename_rasterplot_caudate_in,"w") as filename_rasterplot_caudate_in:
    for i_c in range(len(spikes_caudate_in_223.t/ms)):
        filename_rasterplot_caudate_in.write(str(spikes_caudate_in_223.t[i_c]/ms)+","+str(spikes_caudate_in_223.i[i_c])+"\n")        



filename_rasterplot_dGPe="00241_rasterplot_dGPe_241_"+filename+".dat"
with open(filename_rasterplot_dGPe,"w") as filename_rasterplot_dGPe:
    for i_c in range(len(spikes_Dorsal_GPe_241.t/ms)):
        filename_rasterplot_dGPe.write(str(spikes_Dorsal_GPe_241.t[i_c]/ms)+","+str(spikes_Dorsal_GPe_241.i[i_c])+"\n") 


filename_rasterplot_dGPi="00242_rasterplot_dGPi_242_"+filename+".dat"
with open(filename_rasterplot_dGPi,"w") as filename_rasterplot_dGPi:
    for i_c in range(len(spikes_Dorsal_GPi_242.t/ms)):
        filename_rasterplot_dGPi.write(str(spikes_Dorsal_GPi_242.t[i_c]/ms)+","+str(spikes_Dorsal_GPi_242.i[i_c])+"\n") 

filename_rasterplot_dSTN="00243_rasterplot_dSTN_243_"+filename+".dat"
with open(filename_rasterplot_dSTN,"w") as filename_rasterplot_dSTN:
    for i_c in range(len(spikes_Dorsal_STN_243.t/ms)):
        filename_rasterplot_dSTN.write(str(spikes_Dorsal_STN_243.t[i_c]/ms)+","+str(spikes_Dorsal_STN_243.i[i_c])+"\n") 

filename_rasterplot_dTHL="00151_rasterplot_dTHL_151_"+filename+".dat"
with open(filename_rasterplot_dTHL,"w") as filename_rasterplot_dTHL:
    for i_c in range(len(spikes_Dorsal_THL_251.t/ms)):
        filename_rasterplot_dTHL.write(str(spikes_Dorsal_THL_251.t[i_c]/ms)+","+str(spikes_Dorsal_THL_251.i[i_c])+"\n")



# ##############################################################################
# ##############################################################################
# ##############################################################################
# ##############################################################################
# ##############################################################################
# 
# 
# from scipy.signal import square, sawtooth, correlate
# from numpy import pi, random
# import numpy
# import matplotlib.pyplot as pilot
# 
# 
# C=fr_cortex/fr_cortex[fr_cortex.argmax()]
#                   
#                   
# V=fr_VTA_DA_131/fr_VTA_DA_131[fr_VTA_DA_131.argmax()]
#                   
# T=fr_Ventral_THL_151/fr_Ventral_THL_151[fr_Ventral_THL_151.argmax()]
#                                   
#                                   
# N=fr_nacc/fr_nacc[fr_nacc.argmax()]
# 
# print("arg max kullanimi")
# print(fr_cortex[fr_cortex.argmax()])
# 
# 
# 
# 
# #print(len(C))
# #print(len(N))
# 
# 
# # calculate cross correlation of the two signals
# xcorrCN = correlate(C, N)
# xcorrVN = correlate(V, N)
# xcorrTN = correlate(T, N)
# 
# # The peak of the cross-correlation gives the shift between the two signals
# # The xcorr array goes from -nsamples to nsamples
# # construct time array
# #t = numpy.linspace(0.0, len(A), 0.1, endpoint=False)
# #dt = numpy.linspace(-t[-1], t[-1], 1)
# 
# #print(xcorr.argmax())
# #recovered_time_shift = dt[xcorr.argmax()]
# 
# #print(recovered_time_shift)
# 
# 
# # force the phase shift to be in [-pi:pi]
# #recovered_phase_shift = 2*pi*(((0.5 + recovered_time_shift/period) % 1.0) - 0.5)
# recovered_phase_shift = 1-(abs(xcorrCN.argmax()-len(C)))/len(C)  #rahmi elibol ekledi....
# #relative_error = (recovered_phase_shift - phase_shift)/(2*pi)
# 
# #print ("Original phase shift: %.2f pi" % (phase_shift/pi))
# print ("Recovered phase shift: %.2f " % (recovered_phase_shift))
# #print ("Relative error: %.4f" % (relative_error))
# 
# 
# x1 = numpy.arange(0, len(C), 1)
# x2 = numpy.arange(0, len(xcorrCN), 1)
# 
# pilot.Figure()
# pilot.subplot(411)
# pilot.step(x1,C)
# pilot.ylabel('Korteks')
# pilot.subplot(412)
# pilot.step(x1,V)
# pilot.ylabel('VTA')
# pilot.subplot(413)
# pilot.step(x1,T)
# pilot.ylabel('Talamus')
# pilot.subplot(414)
# pilot.step(x1,N)
# pilot.ylabel('NAc')
# pilot.show()
# 
# 
# pilot.Figure()
# pilot.subplot(311)
# pilot.step(x2,xcorrCN)
# pilot.ylabel('Kortex & NAc')
# pilot.subplot(312)
# pilot.step(x2,xcorrVN)
# pilot.ylabel('VTA & NAc')
# pilot.subplot(313)
# pilot.step(x2,xcorrTN)
# pilot.ylabel('Ventral_THL_151 & NAc')
# pilot.show()
# 
# pilot.Figure()
# pilot.step(x1,fr_D1)
# pilot.step(x1,fr_D2)
# pilot.ylabel('D1 and D2 fr')
# pilot.savefig('example2.pdf', dpi=300)
# pilot.show()
# #==============================================================================
# 
# import numpy as np
# Fsamp=100
# 
# pilot.Figure()
# pilot.subplot(311)
# pilot.step(x1,N)
# pilot.subplot(312)
# fr_nacc_av=fr_nacc-np.average(fr_nacc)
# pilot.psd(fr_nacc_av,NFFT=4096,Fs=Fsamp,label='LFP_av',color='k')
# pilot.subplot(313)
# pilot.specgram(fr_nacc_av,NFFT=1024, Fs=Fsamp)
# pilot.show()
# 
# 
#==============================================================================






#########################################
#############---------LFP-------#########
#########################################
        
        
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



print('Elektrot 1 (Core): ('+str(xE11)+','+str(yE11)+','+str(zE11)+')')
print('Elektrot 2 (Shell): ('+str(xE21)+','+str(yE21)+','+str(zE21)+')')
print('Elektrot 3 (Caudate): ('+str(xE3)+','+str(yE3)+','+str(zE3)+')')








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
    if i_lfp<50:  #NAcc IN grubu dahil
        d_nacc_in_E11=((xE11-xN_nacc_in[i_lfp])*(xE11-xN_nacc_in[i_lfp])+(yE11-yN_nacc_in[i_lfp])*(yE11-yN_nacc_in[i_lfp])+(zE11-zN_nacc_in[i_lfp])*(zE11-zN_nacc_in[i_lfp]))
        d_nacc_in_E21=((xE21-xN_nacc_in[i_lfp])*(xE21-xN_nacc_in[i_lfp])+(yE21-yN_nacc_in[i_lfp])*(yE21-yN_nacc_in[i_lfp])+(zE21-zN_nacc_in[i_lfp])*(zE21-zN_nacc_in[i_lfp]))
        LFP_Is_vE1=LFP_Is_vE1+trace_I_s_msnd1_core[i_lfp].Is/d_msnd1_core_E11+trace_I_s_msnd2_core[i_lfp].Is/d_msnd2_core_E11+trace_I_s_msnd1_shell[i_lfp].Is/d_msnd1_shell_E11+trace_I_s_msnd2_shell[i_lfp].Is/d_msnd2_shell_E11+trace_I_s_nacc_in[i_lfp].Is/d_nacc_in_E11
        LFP_Is_vE2=LFP_Is_vE2+trace_I_s_msnd1_core[i_lfp].Is/d_msnd1_core_E21+trace_I_s_msnd2_core[i_lfp].Is/d_msnd2_core_E21+trace_I_s_msnd1_shell[i_lfp].Is/d_msnd1_shell_E21+trace_I_s_msnd2_shell[i_lfp].Is/d_msnd2_shell_E21+trace_I_s_nacc_in[i_lfp].Is/d_nacc_in_E21
    else:         #NAcc IN grubu yok
        LFP_Is_vE1=LFP_Is_vE1+trace_I_s_msnd1_core.Is[i_lfp]/d_msnd1_core_E11+trace_I_s_msnd2_core.Is[i_lfp]/d_msnd2_core_E11+trace_I_s_msnd1_shell.Is[i_lfp]/d_msnd1_shell_E11+trace_I_s_msnd2_shell.Is[i_lfp]/d_msnd2_shell_E11
        LFP_Is_vE2=LFP_Is_vE2+trace_I_s_msnd1_core.Is[i_lfp]/d_msnd1_core_E21+trace_I_s_msnd2_core.Is[i_lfp]/d_msnd2_core_E21+trace_I_s_msnd1_shell.Is[i_lfp]/d_msnd1_shell_E21+trace_I_s_msnd2_shell.Is[i_lfp]/d_msnd2_shell_E21
    


for i_lfp in range(number_of_neurons_in_msnd1_caudate_221):
    d_msnd1_caudate_E3=((xE3-xN_msnd1_caudate[i_lfp])*(xE3-xN_msnd1_caudate[i_lfp])+(yE3-yN_msnd1_caudate[i_lfp])*(yE3-yN_msnd1_caudate[i_lfp])+(zE3-zN_msnd1_caudate[i_lfp])*(zE3-zN_msnd1_caudate[i_lfp]))
    d_msnd2_caudate_E3=((xE3-xN_msnd2_caudate[i_lfp])*(xE3-xN_msnd2_caudate[i_lfp])+(yE3-yN_msnd2_caudate[i_lfp])*(yE3-yN_msnd2_caudate[i_lfp])+(zE3-zN_msnd2_caudate[i_lfp])*(zE3-zN_msnd2_caudate[i_lfp]))
    if i_lfp<50:  #NAcc IN grubu dahil
        d_caudate_in_E3=((xE3-xN_caudate_in[i_lfp])*(xE3-xN_caudate_in[i_lfp])+(yE3-yN_caudate_in[i_lfp])*(yE3-yN_caudate_in[i_lfp])+(zE3-zN_caudate_in[i_lfp])*(zE3-zN_caudate_in[i_lfp]))
        LFP_Is_vE3=LFP_Is_vE3+trace_I_s_msnd1_caudate[i_lfp].Is/d_msnd1_caudate_E3+trace_I_s_msnd2_caudate[i_lfp].Is/d_msnd2_caudate_E3+trace_I_s_caudate_in[i_lfp].Is/d_caudate_in_E3
    else:
        LFP_Is_vE3=LFP_Is_vE3+trace_I_s_msnd1_caudate[i_lfp].Is/d_msnd1_caudate_E3+trace_I_s_msnd2_caudate[i_lfp].Is/d_msnd2_caudate_E3
        
        

#########################################
#############------LFP (end)----#########
#########################################
#########################################
############------LFP Graphs----#########
#########################################







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



filename_LFP_Is="LFP_Is_"+filename+".dat"
with open(filename_LFP_Is,"w") as results_LFP_Is:
    for i in range(len(LFP_Is_vE1)):
 		 	results_LFP_Is.write(str(LFP_Is_vE1[i]/amp)+","+str(LFP_Is_vE2[i]/amp)+","+str(LFP_Is_vE3[i]/amp)+"\n")

show()



#########################################
#########--- LFP Graphs-- (end) --#######
#########################################








print("----------------------------------------------")
print(" - - -     finish     - - - ")
print("----------------------------------------------")
