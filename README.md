# Malla_nanofibras_fortran

Proyecto en fortran para realizar operaciones con mallas de fibras con mayor rapidez.
Principalmente interseccion de fibras (muy costoso), simplificacion de malla (de mallas curvas a mallas tipo resorte) y equilibrbio mecanico.

# Formatos de archivos

Hay dos formatos, Malla Completa y Malla Simplificada.
En el primer caso es la malla como se deposita por el algoritmo de deposición virtual.
El segundo caso es luego de subdividir las fibras largas en las intersecciones y quedarse con fibras más cortes entre dos puntos de unión.

DETALLE: La numeración de fibras y nodos en los archivos es de base 0, porque el programa original estaba codificado en Python.
La numeración internamente en este programa de Fortran90 es de base 1, así que hay que tener cuidado con eso.

# Unidades 

Se trabaja en un sistema de unidades de: micron, kg, seg.
En este sistema, las fuerzas quedan en microNewton, y las tensiones en MPa.

# Archivo de configuración  
##(ConfigurationFile.txt) 

Sirve para seleccionar las opciones de la corrida.
Se usan etiquetas con asteriscos.

### * Numero de acciones
Es la primera etiqueta que se busca y se lee.

Líneas: 

1. [Num_acciones (integer)]. Integer indicando cuántas etiquetas más se van a buscar en el archivo.
2. [List_acciones (integer(Num_acciones))]. Lista numerada de cuáles de todas las acciones que se van a buscar, se aplicarán (es para agilizar la modificación del archivo para realizar diferentes operaciones sin tener que cambiar todo el texto).

EJemplo: 

<pre>
* Numero de acciones 
3
1	3	4 
</pre>

En este caso se van a realizar 3 acciones, respectivamente en las etiquetas 1 3 y 4. Puede estar desordenado e incluso haber repeticiones, aunque no se recomienda.
NOTA: claramente las demás etiquetas necesitan tener numeración.

### * Parametros Constitutivos 

Acá se dan los parámetros constitutivos. No hace falta que esté esto si se van a calcular intersecciones o simplificar una malla. Pero para calcular el equilibrio bajo deformación, si esto no está, va a saltar algún error con seguridad.

Líneas: 
1. [nParamCon]. Número de parámetros constitutivos 
2. [ParamCon]. Array con parámetros constitutivos, el primer valor es un selector de ley constitutiva.

Ejemplo: 

<pre>
* Parametros Constitutivos 
5
3   2.9d3   2.9d0   0.1d0   1.05
</pre>

### * Intersectar
Se calculan las intersecciones de la malla y se agregan los nodos correspondientes, partiendo los segmentos intersectados en dos y rearreglando la conectividad de las fibras para que la malla siga teniendo sentido.
Es una operación sobre el objeto "Malla Completa", tanto input como output (el output tiene mismo nombre pero se agrega "_i" antes del ".txt".

Lineas:

1. "Intersectar"
2. [Opcion_archivo (integer)] [Archivo (character(len=120))]. Opcion archivo indica si se da una malla (1) o un archivo con una lista de mallas (2). Archivo es el nombre del archivo a abrir, cualquiera sea la opción. Si es una lista de mallas, en cada línea se pone un nombre de archivo de malla presente en la carpeta, y en la primera línea de todas va el número de mallas total.
3. [Num_pasadas (integer)]. Numero máximo de pasadas de interseccion que se van a hacer. Esto es porque en cada pasada el algoritmo solo puede capturar una interseccion por segmento. Entonces para poder tener varias por segmento se hacen varias pasadas. Si se advierte que no se encuentra ninguna intersección en alguna pasada, se completa la búsqueda con éxito.
4. [Periodicidad (integer)]. Opcion indicando si se toman intersecciones entre la ultima y primera capas (1) o no (0).
	
Ejemplo: 

<pre>
* 2 
Intersectar 
1	Malla_completa.txt 
10
1
</pre>

En este caso se tiene que la etiqueta Nro. 2 es "Intersectar".
Se da un archivo de malla, de nombre "Malla\_completa.txt". El output va a ser un archivo de malla completa de nombre "Malla\_completa\_i.txt".
Se harán como máximo 10 pasadas buscando intersecciones y se tendrá en cuenta intersecciones entre la primera y última capa.

### * Simplificar
A partir de una malla de tipo "Malla completa" intersectada se calcula una malla de tipo "Malla Simple".
Los nombres de los archivos a leer se modifican para agregar "\_i" antes de la extensión ".txt". Hay que tener cuidado al poner el nombre en ConfigurationFile.txt o en la lista de archivos, ya que ahí va sin el "\_i".
Los archivos modificados se escriben con un "\_s" previo a la extensión ".txt", por lo que quedan finalizando en "\_i\_s.txt".
Simplificar la malla implica deshacerse de las fibras largas compuestas por concatenaciones de segmentos lineales pequeños, y conformar una malla de fibras cortas entre dos puntos de unión (intersecciones entre fibras largas).
Además se agregan parámetros constitutivos para la tensión ingenieril de las fibras.

Líneas 

1. "Simplificar" 
2. [Opcion_archivo (integer)] [Archivo (character(len=120))].
3. [nParamCon]. Número de parámetros constitutivos 
4. [ParamCon]. Array con parámetros constitutivos, el primer valor es un selector de ley constitutiva.

#### Leyes constitutivas programadas 

Según el valor de selector ( ParamCon(1) ) se toma una ley constitutiva u otra, por lo que cambian también el significado del resto de los parámetros constitutivos, así como su cantidad necesaria.

1. Ley bilineal elástica. Parámetros: k1, k2.
	+ k1 es la constante elástica en fuerza de fibras rectas a la tracción.
	+ k2 es la constante elástica en fuerza de fibras enruladas.
2. Ley bilineal elástica. Parámetros: Et, Eb.
	+ Et es el módulo elástico de tensión ingenieril de fibras rectas a la tracción.
	+ Eb es el módulo elástico de tensión ingenieril de fibras enruladas.

Ejemplo: 

<pre>
* 4
Simplificar 
1	Malla_completa.txt 
3
2	2.9d9	2.9d6
</pre>

En este caso se tiene que la etiqueta Nro. 4 es "Simplificar".
Se da un archivo de malla, de nombre "Malla\_completa.txt", por lo que se va a leer al archivo "Malla\_completa\_i.txt".
El output va a ser un archivo de malla completa de nombre "Malla\_completa\_i\_s.txt".
Además se dan 3 parámetros constitutivos en la línea siguiente, de los cuales el primero indica que se usa la ley constitutiva número 2, con los dos parámetros necesarios siguiendo.


### * Equilibrar 
A partir de una malla simplificada se calcula el equilibrio **elástico** dada una deformación (tensor F macroscópico, que se aplica a los bordes).
Si hay deformación plástica en las fibras se mantiene constante, no hay un avance temporal en esta instrucción.
También se puede dar un archivo con muchos F para calcular una curva constitutiva.

Líneas 

1. "Equilibrar" 
2. [Opcion_archivo (integer)] [Archivo (character(len=120))].
3. [npasos]. Número de pasos de vibración.
4. [vec_veces]. Número de veces que se vibra en cada paso.
5. [vec_drmags]. Magnitud del desplazamiento de cada nodo en cada paso.
6. [fza_ref] [fza_tol]. Fuerza de referencia y fuerza tolerancia (para calcular el equilibrio). La fuerza de referencia se utiliza para linealizar el equilibrio y hallar los desplazamientos nodales. La fuerza de tolerancia indica el valor umbral para el cual se puede considerar que un nodo esta en equilibrio. Estos valores se deben obtener acorde a los parametros constitutivos de las nanofibras.
7. [Opcion_F(integer)]. Opción si se da un solo F (1) o si se da un archivo con un F por línea para calcular una curva constitutiva (2).
8. [Fmacro] o [archivo_F]. Si es [F_macro] se dan los 4 valores en orden: F11, F21, F12, F22 (orden Fortran, por columnas). Si es [archivo_F] se debe tener un archivo con el número de Fs en la primera línea y un F por línea. En el último caso, los archivos de salida tienen una numeración agregada al final así se guardan todas las mallas equilibradas sin sobreescribir archivos.

Ejemplo: 

<pre>
* 1
Equilibrar 
1	Malla_completa.txt 
4
10    10    10    10
10.d0    1.d0    0.1d0    0.01d0
100.d0	0.01d0
1
1.1d0    0.d0    0.d0    1.d0
</pre>

En este caso, la instrucción de la etiqueta número 1 es "Equilibrar".
Se lee una malla simple en el archivo "Malla\_completa\_i\_s.txt" ("\_s" indica que no es malla completa sino que se ha simplificado en algúna instrucción previa).
Luego se indica que para equilibrar se van a llevar a cabo 4 pasos de vibración. Cada paso tiene 10 iteraciones y las magnitudes de desplazamientos en cada paso son 10, 1, 0.1, 0.01.
Luego, el 1 indica que el tensor de deformaciones se da en la linea siguiente, el cual se pone en orden típico de Fortran: F11, F21, F12, F22.


### * Traccion 

Esta instruccion es para llevar a cabo una simulacion de un ensayo de traccion con F11 y F22 determinados en cada paso. Se realiza un esquema explícito temporal. En cada paso de tiempo se calcula el equilibrio de la malla y la evolución de los parámetros que puedan variar con el tiempo (plasticidad, rotura de fibras, etc.).

Líneas:

1. Traccion 
2. [Opcion_archivo (integer)] [Archivo (character(len=120))]. Si Opcion_archivo = 3, en este caso, se empieza desde una malla previamente deformada, con un Fmacro dado en la malla.
3. [npasos]. Número de pasos de vibración (para cada paso temporal se calcula un equilibrio mediante método vibracional). 
4. [vec_veces]. Número de veces que se vibra en cada paso.
5. [vec_drmags]. Magnitud del desplazamiento de cada nodo en cada paso.
6. [fza_ref] [fza_tol]. Fuerza de referencia y fuerza tolerancia (ver en sección "Equilibrar").
7. [dT] [dotF11] [dotF22] [F11fin]. Paso temporal, tasa de deformación axial, tasa de deformación transversal y deformación axial final.
8. [Archivo_curva]. Nombre de archivo de curva constitutiva para output de Fmacro y Tmacro.
9. [Opcion_guardar] [dF_guardar]. Indica si se guardan mallas intermedas (si se escriben los archivos, 1=si) y con que intervalos en F11 se hace.

Ejemplo: 

<pre>
* 1
Traccion 
1	Malla_completa.txt 
4
10    10    10    10
10.d0    1.d0    0.1d0    0.01d0
100.d0	0.01d0
0.1d0    0.01d0    -0.01d0    1.2d0
1	0.01d0
</pre>

### * Uniaxial 

Esta instruccion es para llevar a cabo una simulacion de un ensayo de traccion uniaxial. Se realiza un esquema explícito temporal. En cada paso de tiempo se calcula el equilibrio de la malla y la evolución de los parámetros que puedan variar con el tiempo (plasticidad, rotura de fibras, etc.).

Líneas:

1. Traccion 
2. [Opcion_archivo (integer)] [Archivo (character(len=120))]. Si Opcion_archivo = 3, en este caso, se empieza desde una malla previamente deformada, con un Fmacro dado en la malla.
3. [npasos]. Número de pasos de vibración (para cada paso temporal se calcula un equilibrio mediante método vibracional). 
4. [vec_veces]. Número de veces que se vibra en cada paso.
5. [vec_drmags]. Magnitud del desplazamiento de cada nodo en cada paso.
6. [fza_ref] [fza_tol]. Fuerza de referencia y fuerza tolerancia (ver en sección "Equilibrar").
7. [dT] [dotF11] [t22] [F11fin]. Paso temporal, tasa de deformación axial, tracción transversal (=0.d0 para tracción uniaxial, otro valor si se desea) y deformación axial final.
8. [Archivo_curva]. Nombre de archivo de curva constitutiva para output de Fmacro y Tmacro.
9. [Opcion_guardar] [dF_guardar]. Indica si se guardan mallas intermedas (si se escriben los archivos, 1=si) y con que intervalos en F11 se hace.

Ejemplo: 

<pre>
* 1
Traccion 
1	Malla_completa.txt 
4
10    10    10    10
10.d0    1.d0    0.1d0    0.01d0
100.d0	0.01d0
0.1d0    0.01d0    -0.01d0    1.2d0
1	0.01d0
</pre>
