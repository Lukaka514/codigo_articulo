import numpy as np
import matplotlib.pyplot as plt 

fs = 44000
error_suavizado = []
error_t = 1/(2*fs)
error_f = []
error_cruce = 1/(3000 * np.sqrt(12))
datos = np.loadtxt("/home/lukaka/Descargas/Escritorio/topicos/dt4.dat", comments=';')
Nt = 50
g_o = 9.81
pendiente_t = 85.8
mejor_error = 1000
mejor_tiempo = []
mejor_frecuencia = []

Nint = []
t = datos[:, 0]
x = datos[:, 1]
f_h = 3000
v = 343
gravedad = []
def cruces_k(x):
    s = np.sign(x)
    s = s[s != 0]
    return int(np.sum(np.diff(s) != 0))
for j in range(Nt):
    tiempo = []
    frecuencia = []
    error_f = []
    error_zc = []
    N = j+2
    intervalo = np.array_split(datos,N, axis = 0)
    for i, df in enumerate(intervalo):
        ti = df [0, 0]
        tf = df [-1, 0]
        DelT = tf - ti
        x = df[:,1]
        n = cruces_k(x)
        f = (n/(DelT * 2))
        t_centro =((ti + tf)/2)
        frecuencia.append(f)
        tiempo.append(t_centro)
        error_zc = (n * error_t)
        error_tc = 1/f
        error_suavizado = (pendiente_t * (error_tc / 2))
        error_f.append(error_suavizado)
    # print(f"Sección {i}:{t_centro} s,  {f} hz")
    error_fabs = (np.array(error_f) + error_zc)*2000
    t_mean = np.mean(tiempo)
    f_mean = np.mean(frecuencia) 
    tiempo_np = np.array(tiempo, dtype=float)
    frecuencia_np = np.array(frecuencia, dtype=float)
    pendiente = np.sum((tiempo_np - t_mean)*(frecuencia_np - f_mean)) / np.sum((tiempo_np-t_mean)**2)
    intercept = f_mean - pendiente * t_mean
    error_pendiente = np.sum(( - t_mean)*(frecuencia_np - f_mean)) / np.sum((tiempo_np-t_mean)**2)
    f_p = pendiente * tiempo_np + intercept
    g = -pendiente * (v/f_h)
    g_error2 = error_pendiente 
    g_error = abs(g_o - g)
    Y = f_h - ((g_o * f_h)/(v))


    print(n)
    if g_error < mejor_error:      
         mejor_N = N
         mejor_tiempo = tiempo
         mejor_frecuencia = frecuencia
         mejor_g = g 
         mejor_pendiente = pendiente
         mejor_intercept = intercept
         mejor_f_p = f_p
         mejor_error = g_error
         mejor_incertudumbre = error_zc
         mejor_error_f = error_fabs
    gravedad.append(g)
    Nint.append(N)
print(mejor_g )
print(mejor_N)  




fig, axes = plt.subplots(1,1 ,figsize=(10,4))
#axes.scatter(Nint, gravedad, color='steelblue', s=30, label='gravedad vs intervalos', zorder=3)
#axes.grid(True)
#axes.set_xlabel("Número de intervalos")
#axes.set_ylabel("Gravedad calculada")
axes.scatter(mejor_tiempo, mejor_frecuencia, color='steelblue', s=30, label='frecuencia vs tiempo', zorder=3)
axes.plot(mejor_tiempo, mejor_f_p, linewidth=2, label=f'mejor_frecuencia = {mejor_pendiente:.2f}_mejor_tiempo + {mejor_intercept:.2f}')
axes.text(0.8, 0.95, f"g = {mejor_g:.2f}m/s²", fontsize=12,
             transform=axes.transAxes,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
axes.set_xlabel("Tiempo s")
axes.set_ylabel("Frecuencia hz")
#axes.errorbar(mejor_tiempo, mejor_frecuencia, yerr=mejor_error_f ,fmt= 'none', ecolor = 'gray', capsize=4)
plt.show()
