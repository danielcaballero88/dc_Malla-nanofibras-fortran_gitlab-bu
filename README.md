# Malla_nanofibras_fortran

Proyecto en fortran para realizar operaciones con mallas de fibras con mayor rapidez.
Principalmente interseccion de fibras (muy costoso), simplificacion de malla (de mallas curvas a mallas tipo resorte) y equilibrbio mecanico.

# Formatos de archivos

Hay dos formatos, Malla Completa y Malla Simplificada.
En el primer caso es la malla como se deposita por el algoritmo de deposición virtual.
El segundo caso es luego de subdividir las fibras largas en las intersecciones y quedarse con fibras más cortes entre dos puntos de unión.

DETALLE: La numeración de fibras y nodos en los archivos es de base 0, porque el programa original estaba codificado en Python.
La numeración internamente en este programa de Fortran90 es de base 1, así que hay que tener cuidado con eso.

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
1 3 4 
</pre>

En este caso se van a realizar 3 acciones, respectivamente en las etiquetas 1 3 y 4. Puede estar desordenado e incluso haber repeticiones, aunque no se recomienda.
NOTA: claramente las demás etiquetas necesitan tener numeración.

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
Los nombres de los archivos a leer se modifican para agregar "\_i" antes de la extensión ".txt". Hay que tener cuidado al poner el nombre en ConfigurationFile.txt o en la lista de archivos.
Los archivos modificados se escriben con un "\_s" previo a la extensión ".txt", por lo que quedan finalizando en "\_i\_s.txt".

Líneas 

1. "Simplificar" 
2. [Opcion_archivo (integer)] [Archivo (character(len=120))].

Ejemplo: 

<pre>
* 4
Simplificar 
1	Malla_completa.txt 
</pre>

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
6. [Opcion_F(integer)]. Opción si se da un solo F (1) o si se da un archivo con un F por línea para calcular una curva constitutiva (2).
7. [Fmacro] o [archivo_F]. Si es [F_macro] se dan los 4 valores en orden: F11, F21, F12, F22 (orden Fortran, por columnas). Si es [archivo_F] se debe tener un archivo con el número de Fs en la primera línea y un F por línea. En el último caso, los archivos de salida tienen una numeración agregada al final así se guardan todas las mallas equilibradas sin sobreescribir archivos.

Ejemplo: 

<pre>
* 1
Equilibrar 
1	Malla_completa.txt 
4
10    10    10    10
10.d0    1.d0    0.1d0    0.01d0
1
1.1d0    0.d0    0.d0    1.d0
</pre>