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
        bool exists;
        double px, py, pz;
        double fx, fy, fz;
        double ax, ay, az;
        double vx, vy, vz;
        double mass;
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
                cerr << "sim-aos invoked with " << numberOfParams  << " parameters.\n";
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
                cerr << "sim-aos invoked with " << numberOfParams  << " parameters.\n";
                cerr << "Arguments:\n";
                cerr << "  num_objects: " << num_objects << "\n" << "  num_iterations: " << num_iterations << "\n";
                cerr << "  random_seed: " << random_seed << "\n" << "  size_enclosure: " << size_enclosure << "\n";
                cerr << "  time_step: " << time_step << "\n";
                return -1;    
        }
        return 0;
}

//Se inicializan los parametros de los objetos
void initParam(object * objects, int num_objects, int random_seed, double size_enclosure ){
        mt19937_64 generator(random_seed);
        uniform_real_distribution <double> dis_uniform(0, size_enclosure);
        normal_distribution <double> dis_normal(1E21, 1E15);
        for (int i = 0; i < num_objects; i++){
                objects[i].exists = true;
                objects[i].px = dis_uniform(generator); 
                objects[i].py = dis_uniform(generator); 
                objects[i].pz = dis_uniform(generator); 
                objects[i].mass = dis_normal(generator);
                objects[i].fx = 0; 
                objects[i].fy = 0; 
                objects[i].fz = 0;
                objects[i].ax = 0; 
                objects[i].ay = 0; 
                objects[i].az = 0;  
                objects[i].vx = 0; 
                objects[i].vy = 0; 
                objects[i].vz = 0; 
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
void generateDocuments(string path, object * objects, double size_enclosure, double time_step, int num_objects){
        ofstream file;
        file.open(path, ofstream::out | ofstream::trunc);
        file  << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects <<"\n";
        for(int i = 0; i < num_objects; i++){
                file << objects[i].px << " " << objects[i].py <<  " " << objects[i].pz 
                <<  " " <<  objects[i].vx <<  " " << objects[i].vy << " " << objects[i].vz
                <<  " " << objects[i].mass <<  "\n";
        }
        file.close();       
}

//Se controlan las colisiones de los objetos
void controlColisions(int num_objects, object * objects){
        for (int i = 0; i < num_objects-1; i++){
                if (objects[i].exists) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (objects[j].exists){
                                        if (objects[i].px == objects[j].px && objects[i].py == objects[j].py
                                        && objects[i].pz == objects[j].pz){
                                                objects[j].exists = false;
                                                objects[i].vx +=  objects[j].vx;
                                                objects[i].vy +=  objects[j].vy; 
                                                objects[i].vz +=  objects[j].vz;  
                                                objects[i].mass += objects[j].mass;                                                              
                                        }
                                }
                        }
                }
        }
}

//Al inicio de cada iteraccion se reinicializan las fuerzas de cada objeto a cero
void resetForces(int num_objects, object * objects){
                for (int i = 0; i < num_objects; i++){
                        if (objects[i].exists) {
                                objects[i].fx = 0;
                                objects[i].fy = 0;
                                objects[i].fz = 0; 
                        }                      
                }
}

//Calculo de las fuerzas, se reliza simultaneamente el calculo de las fuerzas ij y ji
void calculateForces(int num_objects, object * objects, double gConst){
        //Calculo de fuerzas
        for (int i = 0; i < num_objects - 1; i++){
                if (objects[i].exists) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (objects[j].exists){                                              
                                        double auxVectorX = objects[j].px - objects[i].px;
                                        double auxVectorY = objects[j].py - objects[i].py;
                                        double auxVectorZ = objects[j].pz - objects[i].pz;
                                        double botPart = (auxVectorX * auxVectorX) + (auxVectorY * auxVectorY) + (auxVectorZ * auxVectorZ); 
                                        botPart = sqrt(botPart);
                                        botPart = botPart * botPart * botPart;
                                        double forceX = (gConst * objects[i].mass * objects[j].mass * auxVectorX)/botPart;
                                        double forceY = (gConst * objects[i].mass * objects[j].mass * auxVectorY)/botPart; 
                                        double forceZ = (gConst * objects[i].mass * objects[j].mass * auxVectorZ)/botPart;  
                                        objects[i].fx  += forceX;  
                                        objects[i].fy  += forceY;                                           
                                        objects[i].fz  += forceZ;
                                        objects[j].fx  += -forceX;
                                        objects[j].fy  += -forceY;  
                                        objects[j].fz  += -forceZ;              
                                }
                        }
                }
        }
}

// Calculo de velocidades, aceleraciones y posiciones.
// Se comprueba ademas que el objeto se encuentra dentro del cubo, de lo contrario se le reposiciona
void calculateParams(int num_objects, object * objects, double size_enclosure, double time_step){
        // Calculo de velocidades, aceleraciones y posiciones
        for (int i = 0; i < num_objects; i++){
                if (objects[i].exists){
                        //Calculo aceleracion, velocidad y posicion de cada objeto  incluyendo colisiones con el contenedor
                        objects[i].ax = objects[i].fx/objects[i].mass;
                        objects[i].ay = objects[i].fy/objects[i].mass; 
                        objects[i].az = objects[i].fz/objects[i].mass; 
                        objects[i].vx += time_step * objects[i].ax; 
                        objects[i].vy += time_step * objects[i].ay;  
                        objects[i].vz += time_step * objects[i].az; 
                        objects[i].px += time_step * objects[i].vx;
                        objects[i].py += time_step * objects[i].vy;
                        objects[i].pz += time_step * objects[i].vz;
                        if(objects[i].px <= 0){
                                objects[i].px = 0;
                                objects[i].vx = objects[i].vx * -1;
                        }
                        else if(objects[i].px >= size_enclosure){
                                objects[i].px = size_enclosure;
                                objects[i].vx = objects[i].vx * -1;
                        }
                        if(objects[i].py <= 0){
                                objects[i].py = 0;
                                objects[i].vy = objects[i].vy * -1;
                        }
                        else if(objects[i].py >= size_enclosure){
                                objects[i].py = size_enclosure;
                                objects[i].vy = objects[i].vy * -1;
                        }
                        if(objects[i].pz <= 0){
                                objects[i].pz = 0;
                                objects[i].vz = objects[i].vz * -1;
                        }
                        else if(objects[i].pz >= size_enclosure){
                                objects[i].pz = size_enclosure;
                                objects[i].vz = objects[i].vz * -1;
                        }
                }
        }
}

//Bucle principal del programa realiza las iteracciones
void iterate(int num_objects, object * objects, double size_enclosure, double time_step, int num_iterations){
        static double gConst = 6.674 / 1E11;
        for (int iteration = 0; iteration < num_iterations; iteration++){
                // Se reinicializan las fuerzas a 0
                resetForces(num_objects, objects);
                // Se calculan las fuerzas sobre cada objeto
                calculateForces(num_objects, objects, gConst);
                // Se calculan aceleracion velocidad y posicion sobre cada objeto
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
        object * objects = new object [num_objects];

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

