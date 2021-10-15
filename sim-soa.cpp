#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <sys/time.h>

// Un struct con arrays p[500][3], v[500],[3].... o con [500][14]
using namespace std;
//Definicion del objeto
struct object {
        bool ** exists;
        double *** f;
        double *** p;
        double *** v;
        double *** a;
        double ** mass;
};


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

int checkParams(int num_objects, int num_iterations, int random_seed, double size_enclosure, double time_step, int numberOfParams){
        if (num_objects < 0 | num_iterations < 0 | random_seed <= 0 | size_enclosure <= 0 | time_step <= 0){
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

void generateDocuments(string path, object objects, double size_enclosure, double time_step, int num_objects){
        ofstream file;
        file.open(path, ofstream::out | ofstream::trunc);
        file  << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects <<"\n";
        for(int i = 0; i < num_objects; i++){
                file << *objects.p[i][0] << " " << *objects.p[i][1] <<  " " << *objects.p[i][2] 
                <<  " " <<  *objects.v[i][0] <<  " " << *objects.v[i][1] << " " << *objects.v[i][2]
                <<  " " << *objects.mass[i] <<  "\n";
        }
        file.close();       
}

void controlColisions(int num_objects, object objects){
        for (int i = 0; i < num_objects-1; i++){
                if (*objects.exists[i]) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (*objects.exists[j]){
                                        if (*objects.p[i][0] == *objects.p[j][0] && *objects.p[i][1] == *objects.p[j][1] 
                                        && *objects.p[i][2] == *objects.p[j][2]){
                                                *objects.mass[i] += *objects.mass[j];
                                                for (int k = 0; k < 3 ; k++){
                                                        *objects.v[i][k] +=  *objects.v[j][k];                                                               
                                                }
                                                *objects.exists[j] = false;
                                        }
                                }
                        }
                }
        }
}


void resetForces(int num_objects, object  objects){
                for (int i = 0; i < num_objects; i++){
                        if (*objects.exists[i] == true) {
                                for (int k = 0; k < 3; k++){
                                        *objects.f[i][k] = 0;
                                }  
                        }                      
                }
}

void calculateForces(int num_objects, object  objects, double gConst){
        //Calculo de fuerzas
        for (int i = 0; i < num_objects - 1; i++){
                if (*objects.exists[i]) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (*objects.exists[j] & (i != j)){                                              
                                        double auxVector[3] = {0};
                                        double botPart = 0;
                                        for (int k = 0; k < 3; k++){
                                                auxVector[k] = *objects.p[j][k] - *objects.p[i][k];
                                                botPart += auxVector[k] * auxVector[k]; 
                                        }
                                        botPart = sqrtf(botPart);
                                        botPart = botPart * botPart * botPart;
                                        for (int k = 0; k < 3 ; k++){
                                                double force = (gConst * *objects.mass[i] * *objects.mass[j] * auxVector[k])/botPart; 
                                                *objects.f[i][k]  += force;  
                                                *objects.f[j][k]  += -force;     
                                        }
                                }
                        }
                }
        }
}

void calculateParams(int num_objects, object  objects, double size_enclosure, double time_step){
        // Calculo de velocidades, aceleraciones y posiciones
        for (int i = 0; i < num_objects; i++){
                if (*objects.exists[i]){
                        //Calculo aceleracion, velocidad y posicion de cada objeto  incluyendo colisiones con el contenedor
                        for (int k = 0; k < 3 ; k++){
                                *objects.a[i][k] = *objects.f[i][k]/ *objects.mass[i];  
                                *objects.v[i][k] += *objects.a[i][k] * time_step;  
                                *objects.p[i][k] += *objects.v[i][k] * time_step;
                                if(*objects.p[i][k] <= 0){
                                        *objects.p[i][k] = 1;
                                        *objects.v[i][k] = *objects.v[i][k] * -1;
                                }
                                else if(*objects.p[i][k] >= size_enclosure){
                                        *objects.p[i][k] = size_enclosure;
                                        *objects.v[i][k] = *objects.v[i][k] * -1;
                                }
                        }
                }
        }
}

void iterate(int num_objects, object  objects, double size_enclosure, double time_step, int num_iterations){
        static double gConst = 6.674 / 10E11;
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
        // Calculo del tiempo de ejecucion
        struct timeval start;
        gettimeofday(&start, NULL);
        long int startms = start.tv_sec * 1000 + start.tv_usec / 1000;
        // Comprobacion del numero de entradas
        int numberOfParams = argc - 1;
        if (checkNumberOfParams(argc, argv, numberOfParams) == -1){
                return -1;
        }
        // Comprobacion de los datos
        int num_objects = stoi(argv[1]);
        int num_iterations = stoi(argv[2]);
        int random_seed = stoi(argv[3]);
        double size_enclosure = stof(argv[4]);
        double time_step = stof(argv[5]);
        if (checkParams(num_objects, num_iterations, random_seed, size_enclosure, time_step, numberOfParams) == -1){
                return -2;
        }
        // Se generan las condiciones iniciales
        mt19937_64 generator(random_seed);
        uniform_real_distribution <double> dis_uniform(0.0, size_enclosure);
        normal_distribution <double> dis_normal(10E21, 10E15);
        object objects;
        objects.exists = (bool**) malloc (num_objects * sizeof(bool));
        objects.p = (double***) malloc (num_objects * 3 * sizeof(double));
        objects.f = (double***) malloc (num_objects * 3 * sizeof(double));
        objects.a = (double***) malloc (num_objects * 3 * sizeof(double));
        objects.v = (double***) malloc (num_objects * 3 * sizeof(double));
        objects.mass = (double**) malloc (num_objects * sizeof(double));
        cout << "funciona1\n";
        cout << "funciona2" << *objects.p[0];
        for (int i = 0; i < num_objects; i++){
                for (int k = 0 ; k < 3; k++){
                        *objects.p[i][k] = dis_uniform(generator); 
                        *objects.v[i][k] = 0; 
                        *objects.a[i][k] = 0; 
                        *objects.f[i][k] = 0; 
                }
                *objects.mass[i] = dis_normal(generator);
                *objects.exists[i] = true;
        }
        cout << "funciona3\n";
        // Se generan el documento de la configuracion inicial
        generateDocuments("./init_config.txt", objects, size_enclosure, time_step, num_objects);

        //Comprobacion de colisiones inicial
        controlColisions(num_objects, objects);

        //Bucle principal
        iterate(num_objects, objects, size_enclosure, time_step, num_iterations);

        // Se generan el documento de la configuracion final
        generateDocuments("./final_config.txt", objects, size_enclosure, time_step, num_objects);
        
        // Calculo del tiempo de ejecucion
        struct timeval end;
        gettimeofday(&end, NULL);
        long int endms = end.tv_sec * 1000 + end.tv_usec / 1000;
        long int timeRunning = endms - startms;
        cout << "Execution time; "<<  timeRunning << " milisecond/s" << endl;      
        return 0;
};

