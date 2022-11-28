# Implementación en paralelo de los algoritmos que especifican las pruebas SAC y BIC

## Autor

| **Nombre y Apellidos** |            **Correo**            | **Grupo** |
| :-----------------------: | :------------------------------: | :-------: |
|    Eziel C. Ramos Piñón   |     ezielramos498@gmail.com      |   C-512   |

## Tutor

|    Evaristo José Madarro Capó  |

## Implementación y Ejecución

### Implementación

La parte computacional del proyecto está implementada completamente en c++. 

#### Agregar nuevo cifrador de flujo

Para agregar un nuevo cifrador a evaluar, añadir la función que define a este cifrador en el archivo main.cpp. 
Todas las variables globales del cifrador deben estar contenidas dentro de esta funcion para evitar [race-condition]((https://en.wikipedia.org/wiki/Race_condition)).

### Ejecución

Para ejecutar el proyecto es necesario escribir las siguientes líneas desde una terminal abierta en esta misma dirección:

```
cd src
g++ main.cpp -o main
.\main.exe
```