#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <sys/time.h>

using namespace std;

//Definicion del objeto donde se tiene: existencia o no del objeto, vectores de posicion, fuerza, aceleracion y velocidad y masa
struct object {
        bool * exists;
        double ** p;
        double ** f;
        double ** a;
        double ** v;
        double * mass;
};

// Se comprueba el numero de argumentos de entrada
int checkNumberOfParams (int argc, char * argv[], int numberOfParams){
        if (numberOfParams!= 5){
                string auxArray[6] = {"?"};
                for (int i = 0; i < 6; i++){
                        auxArray[i] = "?";     
                }
                for (int i = 0; i < argc; i++){
                        auxArray[i] = (string)argv[i];
                }
                cerr << "sim-soa invoked with " << numberOfParams  << " parameters.\n";
                cerr << "Arguments:\n";
                cerr << "  num_objects: " << auxArray[1] << "\n" << "  num_iterations: " << auxArray[2] << "\n";
                cerr << "  random_seed: " << auxArray[3] << "\n" << "  size_enclosure: " << auxArray[4] << "\n";
                cerr << "  time_step: " << auxArray[5] << "\n";
                return -1;      
        }
        return 0;
}

//Se comprueba que los parametros de entrada son validos
int checkParams(int num_objects, int num_iterations, int random_seed, double size_enclosure, double time_step, int numberOfParams){
        if ((num_objects < 0) | (num_iterations < 0) | (random_seed <= 0) | (size_enclosure <= 0) | (time_step <= 0)){
                if (num_objects <= 0){
                        cerr << "Error: Invalid number of objects\n";
                }
                if (num_iterations <= 0){
                        cerr << "Error: Invalid number of iterations\n";
                }
                if (random_seed <= 0){
                        cerr << "Error: Invalid random seed\n";
                }
                if (size_enclosure <= 0){
                        cerr << "Error: Invalid box size\n";
                }
                if (time_step <= 0){
                        cerr << "Error: Invalid time step\n";
                }
                cerr << "sim-soa invoked with " << numberOfParams  << " parameters.\n";
                cerr << "Arguments:\n";
                cerr << "  num_objects: " << num_objects << "\n" << "  num_iterations: " << num_iterations << "\n";
                cerr << "  random_seed: " << random_seed << "\n" << "  size_enclosure: " << size_enclosure << "\n";
                cerr << "  time_step: " << time_step << "\n";
                return -1;    
        }
        return 0;
}

//Se inicializan los parametros de los objetos
void initParam(object objects, int num_objects, int random_seed, double size_enclosure ){
        mt19937_64 generator(random_seed);
        uniform_real_distribution <double> dis_uniform(0, size_enclosure);
        normal_distribution <double> dis_normal(1E21, 1E15);
        for (int i = 0; i < num_objects; i++){
                objects.p[i][0] = dis_uniform(generator); 
                objects.p[i][1] = dis_uniform(generator); 
                objects.p[i][2] = dis_uniform(generator); 
                objects.v[i][0] = 0; 
                objects.v[i][1] = 0; 
                objects.v[i][2] = 0; 
                objects.a[i][0] = 0; 
                objects.a[i][1] = 0; 
                objects.a[i][2] = 0; 
                objects.f[i][0] = 0; 
                objects.f[i][1] = 0; 
                objects.f[i][2] = 0; 
                objects.mass[i] = dis_normal(generator);
                objects.exists[i] = true;
        }

}

//Se muestran por la salida estandar los parametros del programa
void initParamExit(int num_objects, int num_iterations, int random_seed, double size_enclosure, double time_step){
                cout << "Creating simulation:\n";
                cout << "  num_objects: " << num_objects << "\n" << "  num_iterations: " << num_iterations << "\n";
                cout << "  random_seed: " << random_seed << "\n" << "  size_enclosure: " << size_enclosure << "\n";
                cout << "  time_step: " << time_step << "\n";
}

//Metodo para generar los documentos de salida con las configuraciones iniciales y finales
void generateDocuments(string path, object objects, double size_enclosure, double time_step, int num_objects){
        ofstream file;
        file.open(path, ofstream::out | ofstream::trunc);
        file  << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects <<"\n";
        for(int i = 0; i < num_objects; i++){
                file << objects.p[i][0] << " " << objects.p[i][1] <<  " " << objects.p[i][2] 
                <<  " " <<  objects.v[i][0] <<  " " << objects.v[i][1] << " " << objects.v[i][2]
                <<  " " << objects.mass[i] <<  "\n";
        }
        file.close();       
}

//Se controlan las colisiones de los objetos
void controlColisions(int num_objects, object objects){
        for (int i = 0; i < num_objects-1; i++){
                if (objects.exists[i]) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (objects.exists[j]){
                                        if (objects.p[i][0] == objects.p[j][0] && objects.p[i][1] == objects.p[j][1] 
                                        && objects.p[i][2] == objects.p[j][2]){
                                        objects.mass[i] +=  objects.mass[j];
                                        objects.v[i][0] +=  objects.v[j][0];
                                        objects.v[i][1] +=  objects.v[j][1];
                                        objects.v[i][2] +=  objects.v[j][2];                                                               
                                        objects.exists[j] = false;
                                        }
                                }
                        }
                }
        }
}

//Al inicio de cada iteraccion se reinicializan las fuerzas de cada objeto a cero
void resetForces(int num_objects, object  objects){
                for (int i = 0; i < num_objects; i++){
                        if (objects.exists[i] == true) {
                                objects.f[i][0] = 0;
                                objects.f[i][1] = 0;
                                objects.f[i][2] = 0;
                        }                      
                }
}

//Calculo de las fuerzas, se reliza simultaneamente el calculo de las fuerzas ij y ji
void calculateForces(int num_objects, object  objects, double gConst){
        //Calculo de fuerzas
        for (int i = 0; i < num_objects - 1; i++){
                if (objects.exists[i]) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (objects.exists[j]){                                              
                                        double auxVectorX = objects.p[j][0] - objects.p[i][0];
                                        double auxVectorY = objects.p[j][1] - objects.p[i][1];
                                        double auxVectorZ = objects.p[j][2] - objects.p[i][2];
                                        double botPart = (auxVectorX * auxVectorX) + (auxVectorY * auxVectorY) + (auxVectorZ * auxVectorZ); 
                                        botPart = sqrt(botPart);
                                        botPart = botPart * botPart * botPart;
                                        double forceX = (gConst * objects.mass[i] *  objects.mass[j] * auxVectorX)/botPart; 
                                        double forceY = (gConst * objects.mass[i] *  objects.mass[j] * auxVectorY)/botPart; 
                                        double forceZ = (gConst * objects.mass[i] *  objects.mass[j] * auxVectorZ)/botPart; 
                                        objects.f[i][0]  += forceX; 
                                        objects.f[i][1]  += forceY;  
                                        objects.f[i][2]  += forceZ;   
                                        objects.f[j][0]  += -forceX; 
                                        objects.f[j][1]  += -forceY; 
                                        objects.f[j][2]  += -forceZ;                      
                                }
                        }
                }
        }
}

// Calculo de velocidades, aceleraciones y posiciones.
// Se comprueba ademas que el objeto se encuentra dentro del cubo, de lo contrario se le reposiciona
void calculateParams(int num_objects, object  objects, double size_enclosure, double time_step){
        for (int i = 0; i < num_objects; i++){
                if (objects.exists[i]){
                        //Calculo aceleracion, velocidad y posicion de cada objeto  incluyendo colisiones con el contenedor
                        objects.a[i][0] = objects.f[i][0] / objects.mass[i];
                        objects.a[i][1] = objects.f[i][1] / objects.mass[i];  
                        objects.a[i][2] = objects.f[i][2] / objects.mass[i];    
                        objects.v[i][0] += time_step * objects.a[i][0];
                        objects.v[i][1] += time_step * objects.a[i][1]; 
                        objects.v[i][2] += time_step * objects.a[i][2];   
                        objects.p[i][0] += time_step * objects.v[i][0];
                        objects.p[i][1] += time_step * objects.v[i][1];
                        objects.p[i][2] += time_step * objects.v[i][2];
                        if(objects.p[i][0] <= 0){
                                objects.p[i][0] = 0;
                                objects.v[i][0] = objects.v[i][0] * -1;
                        }
                        else if(objects.p[i][0] >= size_enclosure){
                                objects.p[i][0] = size_enclosure;
                                objects.v[i][0] = objects.v[i][0] * -1;
                        }
                        if(objects.p[i][1] <= 0){
                                objects.p[i][1] = 0;
                                objects.v[i][1] = objects.v[i][1] * -1;
                        }
                        else if(objects.p[i][1] >= size_enclosure){
                                objects.p[i][1] = size_enclosure;
                                objects.v[i][1] = objects.v[i][1] * -1;
                        }
                        if(objects.p[i][2] <= 0){
                                objects.p[i][2] = 0;
                                objects.v[i][2] = objects.v[i][2] * -1;
                        }
                        else if(objects.p[i][2] >= size_enclosure){
                                objects.p[i][2] = size_enclosure;
                                objects.v[i][2] = objects.v[i][2] * -1;
                        }
                }
        }
}

//Bucle principal del programa realiza las iteracciones
void iterate(int num_objects, object  objects, double size_enclosure, double time_step, int num_iterations){
        static double gConst = 6.674 / 10E11;
        for (int iteration = 0; iteration < num_iterations; iteration++){
                // Se reinicializan las fuerzas a 0
                resetForces(num_objects, objects);
                // Se calculan las fuerzas sobre cada objeto
                calculateForces(num_objects, objects, gConst);
                // Se calculan aceleracion velocidad y posicion sobre cada objeto, ademas de reposicionamientos si es necesario
                calculateParams(num_objects, objects, size_enclosure, time_step);
                // Se controlan las colisiones de los objetos
                controlColisions(num_objects, objects);
        }      
}

int main (int argc, char * argv[]){
        // Comprobacion del numero de entradas
        int numberOfParams = argc - 1;
        if (checkNumberOfParams(argc, argv, numberOfParams) == -1){
                return -1;
        }
        // Comprobacion de los datos
        int num_objects = stoi(argv[1]);
        int num_iterations = stoi(argv[2]);
        int random_seed = stoi(argv[3]);
        double size_enclosure = stod(argv[4]);
        double time_step = stod(argv[5]);
        if (checkParams(num_objects, num_iterations, random_seed, size_enclosure, time_step, numberOfParams) == -1){
                return -2;
        }
        // Se generan las estructuras de datos
        object objects;
        objects.exists = new bool[num_objects];
        objects.p = new double *[num_objects];
        objects.f = new double *[num_objects];
        objects.a = new double *[num_objects];
        objects.v = new double *[num_objects];
        for (int i = 0; i < num_objects; i++){
                objects.p[i] = new double[3];
                objects.f[i] = new double[3];
                objects.a[i] = new double[3];
                objects.v[i] = new double[3];
        }
        objects.mass = new double[num_objects];

        // Se generan las condiciones iniciales
        initParam(objects, num_objects, random_seed, size_enclosure);
        
        // Se imprimen los parametros 
        initParamExit(num_objects, num_iterations, random_seed, size_enclosure, time_step);

        // Se generan el documento de la configuracion inicial
        generateDocuments("./init_config.txt", objects, size_enclosure, time_step, num_objects);

        //Comprobacion de colisiones inicial
        controlColisions(num_objects, objects);

        //Bucle principal
        iterate(num_objects, objects, size_enclosure, time_step, num_iterations);

        // Se generan el documento de la configuracion final
        generateDocuments("./final_config.txt", objects, size_enclosure, time_step, num_objects);

        return 0;
}

