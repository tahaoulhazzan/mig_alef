# +
import CoolProp.CoolProp as cp
import numpy as np
from matplotlib import pyplot as plt

def p_inlet(t, p_i_fcv, aprr): # calcule le pressure inlet au niveau du dispenser
    return p_i_fcv+(aprr/60)*t*1e6
# p_i_fcv : pression dans le réservoir fuel cell vehicle
# aprr : average pressure rate ramp d'après la norme SAE J 2601/2

def redvalve(p_i, p_o, T_in, kp, rho_o): # calcule le débit massique à partir de Delta pression
    del_p = abs(p_i-p_o)/1e5
    rho_i = cp.PropsSI('Dmass', 'T', T_in, 'P', p_i, 'H2')
    kp = kp
    vdot = (2*del_p/(kp*rho_i))**0.5
    mdot = (vdot/3600)*rho_o
    return mdot

nb_reservoirs = 5 #nombre de réservoirs à charger successivement
t=0
dt = 0.1
aprr = 7.5
T_ambient = 25 + 273.15
p_fcv_ini = 20e5
T_ini = T_ambient
V_fcv = 0.35 #350litres par réservoir
kp_valve = 0.035
cascade = [[300e5,50,2.435], [300e5,50,2.435], [300e5,50,2.435], [300e5,50,2.435],
           [450e5,50,1.758], [450e5,50,1.758], [450e5,50,1.758], [450e5,50,1.758]]
time_array = np.array([])
mdot_array = np.array([])
pin_array = np.array([]) #pression in (à la sortie des tanks du cascade storage)
cooling_array = np.array([]) #puissance de refroidissement du H2
cascade_track = [ [ np.array([]), np.array([]), np.array([]) ] ] #tracking des paramètres du cascade storage en décharge
fcv_track =  [ ] #tracking des paramètres des réservoirs du fcv en recharge

for reservoir in range(nb_reservoirs)  :
    #initialisation des paramètres du réservoir du fcv
    dm_dt = 0
    du_dt_fcv = 0
    p_fcv = p_fcv_ini
    u_fcv = cp.PropsSI('U', 'P', p_fcv_ini, 'T', T_ini, 'H2')
    m_fcv = V_fcv*cp.PropsSI('D', 'P', p_fcv_ini, 'T', T_ini, 'H2')
    fcv_track += [ [ t, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]) ] ]
    #initialisation des paramètres du tank du cascade storage system
    stage = 0
    p_tank = cascade[stage][0]
    m_tank = cascade[stage][1]
    T_tank = T_ambient
    u_tank = cp.PropsSI('U', 'P', p_tank, 'T', T_tank, 'H2')
    du_dt_tank = 0
    #début recharge du réservoir
    while p_fcv<350e5:
        u_fcv += du_dt_fcv*dt
        u_tank += du_dt_tank*dt
        m_fcv += dm_dt*dt
        m_tank -= dm_dt*dt
        rho_fcv = m_fcv/V_fcv
        rho_tank = m_tank/cascade[stage][2]
        p_aprr = p_inlet(t-fcv_track[reservoir][0], p_fcv_ini, aprr)
        p_fcv = cp.PropsSI('P', 'U', u_fcv, 'Dmass', rho_fcv, 'H2')
        h_tank = cp.PropsSI('H', 'P', p_tank, 'T', T_tank, 'H2' )
        rho_m = cp.PropsSI('D', 'P', p_fcv, 'H', h_tank, 'H2') #insentalpic
        T_i = - 30 + 273.15
        dm_dt = min(redvalve(p_aprr, p_fcv, T_i, kp_valve, rho_fcv), redvalve(p_tank, p_aprr, T_tank, kp_valve, rho_m))
        hin = cp.PropsSI('H', 'P', p_aprr, 'T', T_i, 'H2')
        p_tank = cp.PropsSI('P', 'U', u_tank, 'Dmass', rho_tank, 'H2')
        cooling = dm_dt*(h_tank-hin)/0.9
        du_dt_fcv = dm_dt*(hin-u_fcv)/m_fcv
        du_dt_tank = dm_dt*(u_tank-h_tank)/m_tank
        T_fcv = cp.PropsSI('T', 'U', u_fcv, 'Dmass', rho_fcv, 'H2')
        T_tank = cp.PropsSI('T', 'U', u_tank, 'Dmass', rho_tank, 'H2')
        time_array = np.append(t, time_array) 
        fcv_track[reservoir][1] = np.append(m_fcv, fcv_track[reservoir][1]) # masse H2 dans le fcv
        fcv_track[reservoir][2] = np.append(T_fcv, fcv_track[reservoir][2]) # température H2 dans le fcv
        fcv_track[reservoir][3] = np.append(p_fcv, fcv_track[reservoir][3]) # pression H2 dans le fcv
        fcv_track[reservoir][4] = np.append(dm_dt, fcv_track[reservoir][4]) # débit H2 dans le fcv
        fcv_track[reservoir][5] = np.append(p_tank, fcv_track[reservoir][5]) # pin H2 dans le fcv
        fcv_track[reservoir][6] = np.append(cooling, fcv_track[reservoir][6]) # pin H2 dans le fcv
        mdot_array = np.append(dm_dt, mdot_array)
        pin_array = np.append(p_tank, pin_array)
        cooling_array = np.append(cooling, cooling_array)
        cascade_track[stage][0] = np.append(p_tank, cascade_track[stage][0])
        cascade_track[stage][1] = np.append(m_tank, cascade_track[stage][1])
        cascade_track[stage][2] = np.append(T_tank, cascade_track[stage][2])
        t += dt
        #condition de switch au tank suivant du cascade storage system
        if p_tank-p_aprr < 1e4 :
            #condition H2 insuffisant
            if stage >= len(cascade)-1 :
                print(f'stock insuffisant pour le bus {(reservoir)//5+1}, réservoir {reservoir+1} chargé à {p_fcv/1e5} bar après {t/60} minutes')
                break
            else :    
                cascade[stage] = [p_tank, m_tank, cascade[stage][2]]
                cascade_track += [ [ np.array([]), np.array([]), np.array([]) ] ]
                stage += 1
                p_tank = cascade[stage][0]
                m_tank = cascade[stage][1]
                T_tank = T_ambient
                u_tank = cp.PropsSI('U', 'P', p_tank, 'T', T_tank, 'H2')
                du_dt_tank = 0
    cascade[stage] = [p_tank, m_tank, cascade[stage][2]]

#fontions pour visualiser
def p_in_plot():
    plt.plot(time_array/60, pin_array/1e5)
    plt.xlabel('Temps ($minutes$)')
    plt.ylabel('Pression IN ($Bar$)')

def cooling_plot():
    plt.plot(time_array/60, cooling_array/1e3)
    plt.xlabel('Temps ($minutes$)')
    plt.ylabel('Puissance de refroidissement ($kiloWatts$)')

def T_fcv_plot(reservoir):
    n = len( fcv_track[reservoir-1][2] )
    l_t=np.linspace( fcv_track[reservoir-1][0], fcv_track[reservoir][0], num = int(n) )
    plt.plot(l_t/60, fcv_track[reservoir-1][2][::-1]-273.15)
    plt.xlabel('Temps ($minutes$)')
    plt.ylabel('Température de $H_2$ dans le réservoir ($Celsius$)')

def m_fcv_plot(reservoir):
    n = len( fcv_track[reservoir-1][1] )
    l_t=np.linspace( fcv_track[reservoir-1][0], fcv_track[reservoir][0], num = int(n) )
    plt.plot(l_t, fcv_track[reservoir-1][1][::-1])

def p_fcv_plot(reservoir):
    n = len( fcv_track[reservoir-1][3] )
    l_t=np.linspace( fcv_track[reservoir-1][0], fcv_track[reservoir][0], num = int(n) )
    plt.plot(l_t, fcv_track[reservoir-1][3][::-1])

def mdot_plot(reservoir):
    n = len( fcv_track[reservoir-1][4] )
    l_t=np.linspace( fcv_track[reservoir-1][0], fcv_track[reservoir][0], num = int(n) )
    plt.plot(l_t, fcv_track[reservoir-1][4][::-1])

def p_tank_plot(tank):
    plt.plot(cascade_track[tank-1][0][::-1])

def m_tank_plot(tank):
    plt.plot(cascade_track[tank-1][1][::-1])

def T_tank_plot(tank):
    plt.plot(cascade_track[tank-1][2][::-1])
# -


