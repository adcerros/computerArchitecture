#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <sys/time.h>

using namespace std;
//Definicion del objeto
struct object {
        bool exists;
        float p[3], v[3], a[3], f[3];
        double mass;
};

void generateDocuments(string path, object * objects, float size_enclosure, float time_step, int num_objects){
        ofstream file;
        file.open(path, ofstream::out | ofstream::trunc);
        file  << size_enclosure << " " << time_step << " " << num_objects <<"\n";
        for(int i = 0; i < num_objects; i++){
                if (objects[i].exists == false){
                        file << "The object " << i << " has merged" << "\n";
                }
                else{
                        file << objects[i].p[0] << " " << objects[i].p[1] <<  " " << objects[i].p[2] 
                        <<  " " <<  objects[i].v[0] <<  " " << objects[i].v[1] << " " << objects[i].v[2]
                        <<  " " << objects[i].mass <<  "\n";
                }
        }
        file.close();       
}

void controlColisions(int num_objects, object * objects){
        for (int i = 0; i < num_objects-1; i++){
                if (objects[i].exists == true) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (objects[j].exists == true){
                                        if (objects[i].p[0] == objects[j].p[0] && objects[i].p[1] && objects[j].p[1] 
                                        && objects[i].p[2] && objects[j].p[2]){
                                                objects[i].mass += objects[j].mass;
                                                for (int k = 0; k < 3 ; k++){
                                                        objects[i].v[k] +=  objects[j].v[k];                                                               
                                                }
                                                objects[j].exists = false;
                                        }
                                }
                        }
                }
        }
}


void resetForces(int num_objects, object * objects, int DIMENSIONS){
                for (int i = 0; i < num_objects; i++){
                        for (int k = 0; k < DIMENSIONS; k++){
                                objects[i].f[k] = 0;
                        }                        
        }
}

void calculateForces(int num_objects, object * objects, int DIMENSIONS){
        float gConst = 6.674 * (1/pow(10,11));
        //Calculo de fuerzas
        for (int i = 0; i < num_objects - 1; i++){
                if (objects[i].exists == true) {
                        for (int j = i + 1; j < num_objects; j++){
                                if (objects[j].exists == true & i != j){                                              
                                        float auxVector[3] = {0};
                                        for (int k = 0; k < DIMENSIONS; k++){
                                                auxVector[k] = objects[j].p[k] - objects[i].p[k];
                                        }  
                                        float botPart = pow(sqrt(pow(auxVector[0],2) + pow(auxVector[1],2) + pow(auxVector[2],2)),3);  
                                        for (int k = 0; k < 3 ; k++){
                                                float force = (gConst * objects[i].mass * objects[j].mass * auxVector[k])/botPart; 
                                                objects[i].f[k]  += force;  
                                                objects[j].f[k]  += force * -1;     
                                        }
                                }
                        }
                }
        }
}

void calculateParams(int num_objects, object * objects, int DIMENSIONS, float size_enclosure, float time_step){
        // Calculo de velocidades, aceleraciones y posiciones
        for (int i = 0; i < num_objects; i++){
                //Calculo aceleracion, velocidad y posicion de cada objeto  incluyendo colisiones con el contenedor
                for (int k = 0; k < DIMENSIONS ; k++){
                        objects[i].a[k] = objects[i].f[k]/objects[i].mass;  
                        objects[i].v[k] += objects[i].a[k] * time_step;  
                        objects[i].p[k] += objects[i].v[k] * time_step;
                        if(objects[i].p[k] <= 0){
                                objects[i].p[k] = 1;
                                objects[i].v[k] = objects[i].v[k] * -1;
                        }
                        else if(objects[i].p[k] >= size_enclosure){
                                objects[i].p[k] = size_enclosure;
                                objects[i].v[k] = objects[i].v[k] * -1;
                        }
                }
        }
}

int main (int argc, char * argv[]){
        // Calculo del tiempo de ejecucion
        struct timeval start;
        gettimeofday(&start, NULL);
        long int startms = start.tv_sec * 1000 + start.tv_usec / 1000;
        
        // Comprobacion del numero de entradas
        int numberOfParams = argc - 1;
        if (numberOfParams!= 5){
                cerr << "sim-soa invoked with " + to_string(numberOfParams) + " parameters.\n";
                return -1;      
        }
       
        // Comprobacion de los datos
        int num_objects = stoi(argv[1]);
        int num_iterations = stoi(argv[2]);
        int random_seed = stoi(argv[3]);
        float size_enclosure = stof(argv[4]);
        float time_step = stof(argv[5]);
        int DIMENSIONS = 3;
        if (num_objects < 0 | num_iterations < 0 | random_seed <= 0 | size_enclosure <= 0 | time_step <= 0){
                cerr << "sim-soa invoked with " + to_string(numberOfParams) + " parameters.\n";
                return -2;    
        }
        
        // Se generan las condiciones iniciales
        mt19937_64 generator(random_seed);
        uniform_real_distribution<float> dis_uniform(0.0, size_enclosure);
        normal_distribution<float> dis_normal(pow(10,21), pow(10,15));
        object objects [num_objects];
        for (int i = 0; i < num_objects; i++){
                for (int k = 0 ; k < DIMENSIONS; k++){
                        objects[i].p[k] = dis_uniform(generator); 
                        objects[i].v[k] = 0; 
                        objects[i].a[k] = 0; 
                        objects[i].f[k] = 0; 
                }
                objects[i].mass = dis_normal(generator);
                objects[i].exists = true;
        }
        
        // Se generan el documento de la configuracion inicial
        generateDocuments("./init_config.txt", objects, size_enclosure, time_step, num_objects);

        //Comprobacion de colisiones inicial
        controlColisions(num_objects, objects);

        //Bucle principal
        for (int iteration = 0; iteration < num_iterations; iteration++){
                // Se reinicializan las fuerzas a 0
                resetForces(num_objects, objects, DIMENSIONS);
                // Se calculan las fuerzas sobre cada objeto
                calculateForces(num_objects, objects, DIMENSIONS);
                // Se calculan aceleracion velocidad y posicion sobre cada objeto
                calculateParams(num_objects, objects, DIMENSIONS, size_enclosure, time_step);
                // Se controlan las colisiones de los objetos
                controlColisions(num_objects, objects);
        }
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

