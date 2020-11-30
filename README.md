## Universidad de Costa Rica
### Escuela de Ingeniería Eléctrica
#### IE0405 - Modelos Probabilísticos de Señales y Sistemas

Segundo semestre del 2020

---

* Estudiante: **Jose Miguel Pizarro Viales**
* Carné: **B86079**
* Grupo: **1**

---
# `P4` - *Modulación digital IQ*
### 4.1. - Modulación QPSK

En esta sección se modificaron los métodos de modulación brindados para que funcionen en la modulación QPSK. Para lograr esto solo fue necesario adaptar los métodos `modulador` y `demodulador`. 

En primer lugar `modulador` se adapto para que ahora divida los bits totales de la imagen en 2 grupos distintos, dependiendo de los indices de cada bit son pares o impares.

```python
    bits_Q = [] #Bits con fase Q
    bits_I = [] #bits con fase I
    
    #Se dividen los bits de 2 en 2
    for i, bit in enumerate(bits):
        if i%2 == 0:
            bits_I.append(bit)
        else:
            bits_Q.append(bit)

    bits_Q = np.array(bits_Q) 
    bits_I = np.array(bits_I)
```

De esta forma se aplica el algoritmo original, solo que ahora hay una señal portadora diferente para cada grupos de bits, esto según el esquema de modulación QPSK. De esta forma se tiene el siguiente algoritmo para la modulación:

```python
    #Numero de bits por cada fase
    N_Q = len(bits_Q)
    N_I = len(bits_I)

    # 2. Construyendo un periodo de las señales portadoras c(t)_I y c(t)_Q
    Tc = 1 / fc  # periodo [s]
    t_periodo = np.linspace(0, Tc, mpp)
    #NUevas señales portadoras
    portadora_Q = np.sin(2*np.pi*fc*t_periodo)
    portadora_I = np.cos(2*np.pi*fc*t_periodo)  

    # 3. Inicializar las señales moduladas s(t)_I y s(t)_Q
    t_simulacion_Q = np.linspace(0, N_Q*Tc, N_Q*mpp)
    t_simulacion_I = np.linspace(0, N_I*Tc, N_I*mpp)
    senal_Tx_I = np.zeros(t_simulacion_I.shape)
    senal_Tx_Q = np.zeros(t_simulacion_Q.shape)
    moduladora_I = np.zeros(t_simulacion_I.shape)  
    moduladora_Q = np.zeros(t_simulacion_Q.shape)

    # 4. Asignar las formas de onda según los bits (QPSK)
    
    #Bits de fase Q
    for i, bit in enumerate(bits_Q):
        if bit == 1:
            senal_Tx_Q[i*mpp : (i+1)*mpp] = portadora_Q
            moduladora_Q[i*mpp : (i+1)*mpp] = 1
        else:
            senal_Tx_Q[i*mpp : (i+1)*mpp] = portadora_Q * -1
            moduladora_Q[i*mpp : (i+1)*mpp] = 0
    
    #Bits de fase I
    for i, bit in enumerate(bits_I):
        if bit == 1:
            senal_Tx_I[i*mpp : (i+1)*mpp] = portadora_I
            moduladora_I[i*mpp : (i+1)*mpp] = 1
        else:
            senal_Tx_I[i*mpp : (i+1)*mpp] = portadora_I * -1
            moduladora_I[i*mpp : (i+1)*mpp] = 0
  
    # Nueva señal modulada
    senal_Tx = senal_Tx_I + senal_Tx_Q
    
    
    # 5. Calcular la potencia promedio de la señal modulada
    Pm_Q = (1 / (N_Q*Tc)) * np.trapz(pow(senal_Tx_Q, 2), t_simulacion_Q)
    Pm_I = (1 / (N_I*Tc)) * np.trapz(pow(senal_Tx_I, 2), t_simulacion_I)
    
    Pm = Pm_I + Pm_Q
```
De esta forma, se modifico a `modulador` para que regrese la señal modulada `Senal_Tx`, la potencia promedio de la señal modulada `Pm`, las señales portadoras `Portadora_I` y `Portadora_Q` y las señales moduladoras `Moduladora_I` y `Moduladora_Q`.


En cuando a `Demodulador`, solo fue necesario aplicar el algoritmo original a ambas fases Q e I al aplicar la demoludación por detección de energía con la señal portadora correspondiente a cada fase. De esta forma se hallaron los bits transportados en cada una de las fases:

```python
    # Cantidad de muestras en senal_Rx
    M = len(senal_Rx)

    # Cantidad de bits en transmisión
    N =  int(M / (mpp))

    # Vector para bits obtenidos por la demodulación
    bits_Rx_Q  = np.zeros(N)
    bits_Rx_I  = np.zeros(N)
    bits_Rx = np.zeros(2*N)

    # Vector para la señal demodulada
    senal_demodulada_Q = np.zeros(M)
    senal_demodulada_I = np.zeros(M)

    # Energía de un período de la portadora
    Es_Q = np.sum(portadora_Q**2)
    Es_I = np.sum(portadora_I**2)

    # Demodulación
    for i in range(N):
        #Demodulación de la fase I
        #Producto de dos funciones
        producto_I = senal_Rx[i*mpp : (i+1)*mpp] * portadora_I
        senal_demodulada_I[i*mpp : (i+1)*mpp] = producto_I
        Ep_I = np.sum(producto_I) 

        if Ep_I > 0:
            bits_Rx_I[i] = 1
        else:
            bits_Rx_I[i] = 0
        
        #Demodulación de la fase Q
        # Producto interno de dos funciones
        producto_Q = senal_Rx[i*mpp : (i+1)*mpp] * portadora_Q
        senal_demodulada_Q[i*mpp : (i+1)*mpp] = producto_Q
        Ep_Q = np.sum(producto_Q) 

        # Criterio de decisión por detección de energía
        if Ep_Q > 0:
            bits_Rx_Q[i] = 1
        else:
            bits_Rx_Q[i] = 0
```

De esta forma solo resta reordenar los bits recuperados según el criterio de división empleado en el modulador. De esta forma se tiene lo siguiente:

```python
    #Recuperación de los bits
    #Contadores
    N_I = 0
    N_Q = 0
    
    #Algoritmo de recuperación
    for i in range(0, len(bits_Rx)):
        if i%2 == 0: #Debido a que los bits pares se enviaron con fase I
            bits_Rx[i] = bits_Rx_I[N_I]
            N_I = N_I+1
        else: #Debido a que los bits impares se enviaron con fase Q
            bits_Rx[i] = bits_Rx_Q[N_Q]
            N_Q = N_Q+1
```

De esta forma, se modifico a `demodulador` para que regrese los bits recuperados `bits_Rx` y las señales demoduladas de cada fase `senal_demodulada_I` y `senal_demodulada_Q`.

Al hacer esto se obtuvieron los siguientes resultados:

<img align='center' src='https://github.com/pizarrin737/Proyecto4/blob/main/Imagen_Recuperada.png?raw=true' width='350'/>

<img align='center' src='https://github.com/pizarrin737/Proyecto4/blob/main/Modulacion.png?raw=true' width='350'/>


### 4.2. - Estacionaridad y ergodicidad

En esta sección se aplicaron pruebas de estacionaridad y ergodicidad a la señal modulada `Senal_Tx`. Para ello, en primer lugar se realizaron pruebas de ergodicidad en el promedio a lo larga de un periodo completo de la onda. Para ello se tomaron en cuenta todas las posibles combinaciones de `Senal_Tx`. Estas están dadas por los simbolos correspondientes a las variables aleatorias dadas por los bits del canal I e Q (1 si el bit es 1 o -1 si el bit es 0). De esta forma se implementaron las 4 cmbinaciones posibles y se calculo el promedio de las mismas.

```python
#Posibles combinaciones de las portadoras
comb1 = portadora_I + portadora_Q
comb2 = -portadora_I - portadora_Q
comb3 = portadora_I - portadora_Q
comb4 = -portadora_I + portadora_Q

combinaciones = np.array([comb1, comb2, comb3, comb4])

#Promedio de las combinaciones
Prom_combinaciones = [np.mean(combinaciones[:,i]) for i in range(len(comb1))]
```

El mismo procedimiento se repitió para `Senal_Tx`, por lo que se obtienen los siguientes resultados:


<img align='center' src='https://github.com/pizarrin737/Proyecto4/blob/main/Promedio.png?raw=true' width='350'/>

De esta forma se compruba que para este periodo, `Senal_Tx` es igual a la combinacion 3 y además, su promedio estadisitico es igual constante a lo largo del tiempo con un valor de 0. Por este motivo se puede afirmar que `Senal_Tx` constituye un proceso aleatrio estacionario de primer orden. Además, por la forma de lo onda que describe `Senal_Tx` en el gráfico, tambien se puede afirmar que su promedio temporal es 0 dado que es una función periodica. De esta forma se tine que el promedio temporal de `Senal_Tx` es igual a su promedio estadístico. Por lo tanto `Senal_Tx` constituye un proceso aleatorio ergódico.

Por otra parte, tambien se analizó la autocorrelación de `Senal_Tx`. Para ello se calculo la correlación para las 4 combinacioes de la siguiente manera:

```python
#Posibles combinaciones de las portadoras
plt.figure(figsize=(14, 7))

# T valores de desplazamiento tau
desplazamiento = np.arange(len(comb1))

# Inicialización de matriz de valores de correlación para las 4 combinaciones
corr_comb = np.empty((4, len(desplazamiento)))


# Cálculo de correlación en las combinaciones para cada valor de tau
for n in range(4):
	for i, tau in enumerate(desplazamiento):
		corr_comb[n, i] = np.correlate(combinaciones[n,:], np.roll(combinaciones[n,:], tau))/len(comb1)
```

El mismo proceso se aplicó a `Senal_Tx` por lo que se obtuvieron los siguientes resultados:

<img align='center' src='https://github.com/pizarrin737/Proyecto4/blob/main/Autocorrelacion.png?raw=true' width='350'/>

Notese que todos los valores de correlación de las combinaciones y de `Senal_Tx` son iguales en cada instaate de tiempo a lo largo de un periodo. De esta forma se demustra que `Senal_Tx` tambien constituye un proceso alatoria con estacionaridad de segundo orden. 

Entonces se concluye que `Senal_Tx` constituye un proceso aletorio estacionario en sentido amplio. Esta conclusión refuerza el hecho de que `Senal_Tx` constituye un proceso aletorio ergódico.


### 4.3. - Densidad espectral de potencia

Dado que en el inciso anterior se demostró que el proceso alatroio descrito por `Senal_Tx` es ampliamente estacionario, se teien que la densidad espectral de potencia está descrita por la siguiente relación:

$$
S_{XX}(\omega) = \lim_{T\rightarrow \infty } \frac{E\left [ \left | X_T(\omega) \right |^2 \right ]}{2T}
$$

De esta forma se empleo el siguiente código para recrear dicha operación:

```python
from scipy import fft

#Fourier
senal_f = fft(senal_Tx)

#Muestras por seañal
Nm = len(senal_Tx)

#Numero de simbolos (198 x 89 x 8 x 3) )
Ns = Nm//mpp

#Tiempo del simbolo
Tc = 1/fc

#Periodo de muestreo
Tm = Tc/mpp

#Tiempo de la simulación
T = Ns*Tc

#Espacio de frecuencias
f = np.linspace(0.0, 1.0/(2.0*Tm),Nm//2)

#Grafico
plt.plot(f,2.0/Nm*np.power(np.abs(senal_f[0:Nm//2]),2))
plt.xlim(0,20000)
plt.grid()
plt.show()
```
Entonces, se obtuviron los siguientes resultados:

<img align='center' src='https://github.com/pizarrin737/Proyecto4/blob/main/Den_Espc_Potencia.png?raw=true' width='350'/>
