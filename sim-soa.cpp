#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace std;
//Definicion del objeto
struct object {
        bool exists;
        float positionVector[3];
        float velocityVector[3];
        float acelerationVector[3];
        float mass;
};

void generateDocuments(string path, object * objects, float size_enclosure, float time_step, int num_objects){
        ofstream file;
        file.open(path, ofstream::out | ofstream::trunc);
        file  << size_enclosure << " " << time_step << " " << num_objects <<"\n";
        for(int i = 0; i < num_objects; i++){
                object current = objects[i];
                file << current.positionVector[0] << " " << current.positionVector[1] <<  " " << current.positionVector[2] 
                <<  " " <<  0 <<  " " << 0 <<  " " << 0
                <<  " " << current.mass <<  "\n";
        }
        file.close();       
};



int main (int argc, char * argv[]){
        
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
        if (num_objects < 0 | num_iterations < 0 | random_seed <= 0 | size_enclosure <= 0 | time_step <= 0){
                cerr << "sim-soa invoked with " + to_string(numberOfParams) + " parameters.\n";
                return -2;    
        }
        // Se inicializa el generador de numeros aleatorios
        mt19937_64 generator(random_seed);
        uniform_real_distribution<float> dis_uniform(0.0, size_enclosure);
        normal_distribution<float> dis_normal(pow(10,21), pow(10,15));

        
        // Se generan las condiciones iniciales
        object objects [num_objects];
        for (int i = 0; i < num_objects; i++){
                object current;
                current.positionVector[0] = dis_uniform(generator);
                current.positionVector[1] = dis_uniform(generator);
                current.positionVector[2] = dis_uniform(generator);
                current.mass = dis_normal(generator);
                current.exists = 1;
                objects[i] = current;
        }
        // Se generan el documento de la configuracion inicial
        generateDocuments("./init_config.txt", objects, size_enclosure, time_step, num_objects);
        
        //Bucle principal
        float gConst = 6.674 * (1/pow(10,11));
        for (int iteration = 0; iteration < num_iterations; iteration++){
                cout << objects[0].positionVector[0] << " " << objects[0].positionVector[1] <<" "<< objects[0].positionVector[2] << objects[0].mass<< "\n";
                //Comprobacion de colisiones
                for (int i = 0; i < num_objects-1; i++){
                        object current = objects[i];
                        if (current.exists == 1) {
                                for (int j = i+1; j< num_objects; j++){
                                        object next = objects[j];
                                        if (next.exists == 1){
                                                if (current.positionVector[0] == next.positionVector[0] && current.positionVector[1] == next.positionVector[1] && current.positionVector[2] == next.positionVector[2]){
                                                        current.mass = current.mass + next.mass;
                                                        current.velocityVector[0] = current.velocityVector[0] + next.velocityVector[0];
                                                        current.velocityVector[1] = current.velocityVector[1] + next.velocityVector[1];
                                                        current.velocityVector[2] = current.velocityVector[2] + next.velocityVector[2];
                                                        next.exists = 0;
                                                }
                                        }
                                }
                        }
                }
                //Calculo de fuerzas, velocidades y aceleraciones
                for (int i = 0; i < num_objects; i++){
                        object current = objects[i];
                        if (current.exists == 1) {
                                float currentForces[3];
                                for (int j = 0; j< num_objects; j++){
                                        object next = objects[j];
                                        if (next.exists == 1 && i != j){
                                                // Se calculan las distancias del objeto i respecto al j
                                                float auxVector[3];
                                                for (int k = 0; k < 3; k++){   
                                                        auxVector[k] = current.positionVector[k] - next.positionVector[k];
                                                }
                                                float botPart = pow(sqrt(pow(auxVector[0],2) + pow(auxVector[1],2) + pow(auxVector[2],2)),3);
                                                for (int k = 0; k < 3; k++){                                                        
                                                        currentForces[k] += (gConst * current.mass * next.mass * auxVector[k])/botPart;
                                                        cout << "object: " << i << " " << "iteration: " << iteration << " " << "coord: " << k << " " << currentForces[k] << "\n";
                                                }
                                        }
                                }
                                //Calculo aceleracion, velocidad y posicion de cada objeto
                                for (int j = 0; j< 3; j++){                                                        
                                        current.acelerationVector[j] = currentForces[j]/current.mass;
                                        current.velocityVector[j] = current.velocityVector[j] + (current.acelerationVector[j] * time_step);
                                        current.positionVector[j] = current.positionVector[j] + (current.velocityVector[j] * time_step);
                                        //colisiones con paredes
                                }
                        }
                }
        }             
        return 0;
};

