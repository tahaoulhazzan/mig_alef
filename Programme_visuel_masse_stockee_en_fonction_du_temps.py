
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
import numpy as np
import math as m


'''
ON CREE LES LISTES DES DIFFERENTS COMPOSANTS
--------------------------------------------
'''


T=293  #K

n_stockage=3       #nombre d'éléments livrés par camion
Pstockage=350e5   #en Pa
list_mass_stock=[350]*n_stockage    #liste des masses
v_in = list_mass_stock[0]/cp.PropsSI('Dmass' ,'T', 293,'P', Pstockage, 'H2')


list_mass_MP=[50]*4
list_pressure_MP=[450e5]*4
V_MP = list_mass_MP[0]/cp.PropsSI('Dmass' ,'T', 293,'P', list_pressure_MP, 'H2')


list_mass_LP=[50]*4
list_pressure_LP=[300e5]*4
V_LP = list_mass_LP[0]/cp.PropsSI('Dmass' ,'T', T,'P', list_pressure_LP[0], 'H2')



nb_passages_bus= 7    #par jour
nb_passages_trains=  0     #par an

reservoirs_bus =  35    #en kg
reservoir_train =  215    #en kg
V_bus = 0
#le compresseur ne fontionne que pour les middle pressure storages
#Je vais donc faire à la main pour les low pressure storage




# REMPLISSAGE DES LP TANKS

'''
La modélisation du remplissage de Low Pressure storage
-------------------------------------------------------
'''



# étudions désormais une fonction qui suit la quantité d'hydrogène dans le LP tank 

#on fait donc des programmes de modélisation de la recharge
def remplissage_LP(list_mass_stock, list_mass_LP, i, debit_normo =500/3600, volume_stockage = 14.8, T = 293) :
    
    print("lp")
    if i >= len(stock_tab)-1 : 
        return list_mass_stock, list_mass_LP, i, stock_tab, LP_tab
    dt = 1
    debit_massique = debit_normo*cp.PropsSI('Dmass', 'T', T, 'P', 1e5,'H2')  ####on peut modifier ici les caractéristiques du compresseur, nous avons pris un débit constant, ce qui est assez grossier
    for k in range (len(list_mass_LP)):  #on parcourt le stock LowPressure

        for n in range (len(list_mass_stock)): #on parcourt le stock
                #si on a une delta P favorable, on fait un remplissage naturel
            pressure = cp.PropsSI('P' ,'T', T,'Dmass', list_mass_stock[n]/v_in, 'H2')
            time = 0
            while pressure - cp.PropsSI('P' ,'T', T,'Dmass', list_mass_LP[k]/V_LP, 'H2') > 0 and list_mass_LP[k] < 50 and list_mass_stock[n] > 60:   
                dm=0.01    #débit massique de 36kg/h (cohérent)
                list_mass_LP[k], list_mass_stock[n] =  list_mass_LP[k] + dm, list_mass_stock[n] - dm
                # Dmass= (list_mass_LP[k] + list_mass_stock[n])/(V_LP + v_in) Mauvaise idée, ça remplit trop
                # list_mass_LP[k], list_mass_stock[n] = V_LP*cp.PropsSI('Dmass' ,'T', T,'P', cp.PropsSI('P' ,'T', T,'Dmass', Dmass, 'H2'), 'H2'), v_in*cp.PropsSI('Dmass' ,'T', T,'P', cp.PropsSI('P' ,'T', T,'Dmass',Dmass, 'H2'), 'H2')
                if abs ( round(time) - time ) < 0.01:
                    i += 1
                    if i >= len(stock_tab) - 1:
                        return list_mass_stock, list_mass_LP, i, stock_tab,LP_tab
                    stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP
                    
                pressure = cp.PropsSI('P' ,'T', T,'Dmass', list_mass_stock[n]/v_in, 'H2')
            i = round(i)
            stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP

            if i < len(temps):
                stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP
                #sinon on fait avec le compresseur
        for n in range (len(list_mass_stock)):
            while list_mass_LP[k] < 50 and list_mass_stock[n] > 60 : 
                list_mass_stock[n] -= debit_massique*dt  # derivee de la masse vaut -debit, m[i+1]=m[i]-debit*dt
                list_mass_LP[k] += debit_massique*dt
                i += dt
                if i >= len(stock_tab) - 1:
                    return list_mass_stock, list_mass_LP, i, stock_tab, LP_tab
                stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP
    
    return list_mass_stock, list_mass_LP, i, stock_tab, LP_tab



''''
La modélisation du remplissage de Middle Pressure storage
-------------------------------------------------------
'''



def remplissage_MP(list_mass_stock, list_mass_MP, i, debit_normo =500/3600,volume_stockage = 14.8,T = 293):   #on définit ce débit à partir de ceux présents dans le commerce https://www.atlascopco.com/fr-fr/compressors/products/gas-compressors/h2y-high-pressure-hydrogen-compressor
    if i >= len(stock_tab):
        return list_mass_stock, list_mass_MP, i, stock_tab, MP_tab
    print("mp")
    
    
    dt = 1
    debit_massique = debit_normo*cp.PropsSI('Dmass', 'T', T, 'P', 1013e2,'H2') ####on peut modifier ici les caractéristiques du compresseur, nous avons pris un débit constant, ce qui est assez grossier
    for k in range (len(list_mass_MP)):  #on parcourt le stck LP

        for n in range (len(list_mass_stock)): #le procesus dure pendant temps secondes
            pressure = cp.PropsSI('P' ,'T', T,'Dmass', list_mass_stock[n]/v_in, 'H2')
            time = 0
            while pressure - cp.PropsSI('P' ,'T', T,'Dmass', list_mass_LP[k]/V_LP, 'H2') > 0 and list_mass_LP[k] < 50 and list_mass_stock[n] > 60:   
                dm=0.01    #débit massique de 36kg/h (cohérent)
                list_mass_LP[k], list_mass_stock[n] =  list_mass_LP[k] + dm, list_mass_stock[n] - dm
                # Dmass= (list_mass_LP[k] + list_mass_stock[n])/(V_LP + v_in) Mauvaise idée, ça remplit trop
                # list_mass_LP[k], list_mass_stock[n] = V_LP*cp.PropsSI('Dmass' ,'T', T,'P', cp.PropsSI('P' ,'T', T,'Dmass', Dmass, 'H2'), 'H2'), v_in*cp.PropsSI('Dmass' ,'T', T,'P', cp.PropsSI('P' ,'T', T,'Dmass',Dmass, 'H2'), 'H2')
                if abs ( round(time) - time ) < 0.01:
                    i += 1
                    if i >= len(stock_tab) - 1:
                        return list_mass_stock, list_mass_MP, i, stock_tab, MP_tab
                    stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP
                pressure = cp.PropsSI('P' ,'T', T,'Dmass', list_mass_stock[n]/v_in, 'H2')
            i = round(i)
            if i >= len(stock_tab) - 1:
                return list_mass_stock, list_mass_MP, i, stock_tab, MP_tab
            stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP

            
                #sinon on fait avec le compresseur
            if i < len(stock_tab):
                stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP
        for n in range (len(list_mass_stock)):
            while list_mass_MP[k] < 50 and list_mass_stock[n] > 60 : 
                list_mass_stock[n] -= debit_massique*dt  
                list_mass_MP[k] += debit_massique*dt
                i += dt
                if i >= len(stock_tab) -1:
                    return list_mass_stock, list_mass_MP, i, stock_tab, MP_tab
                stock_tab[i], MP_tab[i], LP_tab[i] = list_mass_stock, list_mass_MP, list_mass_LP

    return list_mass_stock, list_mass_MP, i, stock_tab, MP_tab











'''
Modélisation du remplissage d'un bus
------------------------------------
'''

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


def plein_bus(nb_reservoirs, list_mass_LP, list_mass_MP, i, T_ambient= 25 + 273.15 ):     #fonction simulant le remplissage d'un bus
    t=0
    dt = 0.1
    aprr = 5      
    p_fcv_ini = 20e5
    T_ini = T_ambient
    V_fcv = 0.35 #350litres par réservoir 
    kp_valve = 0.035    ###d'où ça vient
    cascade = [ [cp.PropsSI('P','Dmass', list_mass_LP[0]/2.435, 'T', T_ini, 'H2'), list_mass_LP[0], 2.435] , [cp.PropsSI('P','Dmass', list_mass_LP[1]/2.435, 'T', T_ini, 'H2') , list_mass_LP[1],2.435] , [cp.PropsSI('P','Dmass', list_mass_LP[2]/2.435, 'T', T_ini, 'H2'), list_mass_LP[2],2.435] , [cp.PropsSI('P','Dmass', list_mass_LP[3]/2.435, 'T', T_ini, 'H2'), list_mass_LP[3],2.435],
                [cp.PropsSI('P','Dmass', list_mass_MP[0]/1.758, 'T', T_ini, 'H2'), list_mass_MP[0], 1.758] , [cp.PropsSI('P','Dmass', list_mass_MP[1]/1.758, 'T', T_ini, 'H2'), list_mass_MP[1], 1.758] , [cp.PropsSI('P','Dmass', list_mass_MP[2]/1.758, 'T', T_ini, 'H2'), list_mass_MP[2], 1.758] , [cp.PropsSI('P','Dmass', list_mass_MP[3]/1.758, 'T', T_ini, 'H2'), list_mass_MP[3], 1.758] ]
    print(cascade)
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
        fcv_track += [ [ t, np.array([]), np.array([]), np.array([]), np.array([]) ] ]
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
            fcv_track[reservoir][3] = np.append(T_fcv, fcv_track[reservoir][3]) # pression H2 dans le fcv
            fcv_track[reservoir][4] = np.append(dm_dt, fcv_track[reservoir][4]) # débit H2 dans le fcv
            mdot_array = np.append(dm_dt, mdot_array)
            pin_array = np.append(p_tank, pin_array)
            cooling_array = np.append(cooling, cooling_array)
            cascade_track[stage][0] = np.append(p_tank, cascade_track[stage][0])
            cascade_track[stage][1] = np.append(m_tank, cascade_track[stage][1])
            cascade_track[stage][2] = np.append(T_tank, cascade_track[stage][2])
            t += dt
            
            if abs(t - round(t))<0.001 :
                i+=1
                if stage in range(4) :
                    list_mass_LP[stage] = m_tank
                else :
                    list_mass_MP[stage%4] = m_tank      #???
                if i >= len(stock_tab)-1 : 
                    return i
                LP_tab[i] = list_mass_LP
                MP_tab[i] = list_mass_MP
                stock_tab[i] = list_mass_stock
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

    return i




##### ICI IL FAUT RENTRER LES CONDITIONS INITIALES
Minit_LP1 = cp.PropsSI('Dmass' ,'T', T,'P', 300e5, 'H2') * V_LP
list_mass_LP = [Minit_LP1]*4


temps=np.arange(3600*15)    #on découpe la journée  de 15h en secondes
stock_tab = np.empty((len(temps), 3))   #tableau des données dans les stockages en f(t)
stock_tab[0] = list_mass_stock
LP_tab = np.empty((len(temps), 4)) #idem pour LP
LP_tab[0] = list_mass_LP
MP_tab = np.empty((len(temps), 4)) #idem pour MP
MP_tab[0] = list_mass_MP

3600*1.75

i=0
cycle = 0
while i < len(temps)-2:
    if i % 3600*2 == 0:   ####toutes les 2h
        cycle+=1
        print('charge bus :', cycle)
        print('i = ', i)
        i = plein_bus(5, list_mass_LP, list_mass_MP, i)
        print (i)
        print('sortie')
    
    #elif i % 3600*2 == 0 or i == 0: #Tous les combien de temps on recharge?
        #cycle+=1
        if cycle % 2 == 0 :
            print('recharge cascade :', cycle)
            list_mass_stock,list_mass_LP, i, stock_tab, LP_tab= remplissage_LP(list_mass_stock, list_mass_LP, i)
            list_mass_stock, list_mass_MP, i, stock_tab, MP_tab = remplissage_MP(list_mass_stock, list_mass_MP, i)
            

    else:
        i += 1
        LP_tab[i] = list_mass_LP
        MP_tab[i] = list_mass_MP
        stock_tab[i] = list_mass_stock




# press = np.vectorize(lambda x,v : cp.PropsSI('P' ,'T', T,'Dmass', x/v, 'H2'))
# LP_press_tab = press(LP_tab, V_LP)
# MP_press_tab = press(MP_tab, V_MP)

temps=temps/3600#on remet en heures
#Création des subplots pour visualiser la quantité
fig, axs = plt.subplots(2, 4, figsize=(12, 6))


stock_tab[-1] = stock_tab[-2]
LP_tab[-1] = LP_tab[-2]
MP_tab[-1] = MP_tab[-2]
# Sous-plot 1
axs[0, 0].plot(temps, stock_tab[:, 0])
axs[0, 0].set_title('stock1')
axs[0, 0].set_xlabel('time (h)')
axs[0, 0].set_ylabel('mass (kg)')


# Sous-plot 2
axs[0, 1].plot(temps, stock_tab[:, 1])
axs[0, 1].set_title('stock2')

# Sous-plot 3
axs[0, 2].plot(temps, LP_tab[:, 0])
axs[0, 2].set_title('LP1')

# Sous-plot 4
axs[0, 3].plot(temps, LP_tab[:, 1])
axs[0, 3].set_title('LP2')

# Sous-plot 5
axs[1, 0].plot(temps, LP_tab[:, 2])
axs[1, 0].set_title('LP_3')

# Sous-plot 6
axs[1, 1].plot(temps, MP_tab[:, 0])
axs[1, 1].set_title('MP1')

# Sous-plot 7
axs[1, 2].plot(temps, MP_tab[:, 1])
axs[1, 2].set_title('MP_2')

# Supprimer le huitième sous-plot (il n'y a pas assez de sous-plots pour 2x4)
fig.delaxes(axs[1, 3])

# Ajuster automatiquement les espacements pour éviter le chevauchement
plt.tight_layout()

# Afficher le graphique
plt.show()





