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
        float px, py, pz;
        float vx, vy, vz;
        float ax, ay, az;
        float fx, fy, fz;
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
                        file << objects[i].px << " " << objects[i].py <<  " " << objects[i].pz 
                        <<  " " <<  0 <<  " " << 0 <<  " " << 0
                        <<  " " << objects[i].mass <<  "\n";
                }
        }
        file.close();       
};



int main (int argc, char * argv[]){
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
                objects[i].px = dis_uniform(generator); 
                objects[i].py = dis_uniform(generator); 
                objects[i].pz = dis_uniform(generator); 
                objects[i].vx = 0; 
                objects[i].vy = 0; 
                objects[i].vz = 0; 
                objects[i].ax = 0; 
                objects[i].ay = 0; 
                objects[i].az = 0; 
                objects[i].fx = 0; 
                objects[i].fy = 0; 
                objects[i].fz = 0; 
                objects[i].mass = dis_normal(generator);
                objects[i].exists = true;
        }
        // Se generan el documento de la configuracion inicial
        generateDocuments("./init_config.txt", objects, size_enclosure, time_step, num_objects);

        //Bucle principal
        float gConst = 6.674 * (1/pow(10,11));
        for (int iteration = 0; iteration < num_iterations; iteration++){
                //Comprobacion de colisiones
                for (int i = 0; i < num_objects-1; i++){
                        if (objects[i].exists == true) {
                                for (int j = i+1; j< num_objects; j++){
                                        if (objects[j].exists == true){
                                                if (objects[i].px == objects[j].px && objects[i].py && objects[j].py && objects[i].pz && objects[j].pz){
                                                        objects[i].mass = objects[i].mass + objects[j].mass;
                                                        objects[i].vx = objects[i].vx + objects[j].vx;
                                                        objects[i].vy = objects[i].vy + objects[j].vy;
                                                        objects[i].vz = objects[i].vz + objects[j].vz;
                                                        objects[j].exists = false;
                                                }
                                        }
                                }
                        }
                }
                //Calculo de fuerzas
                for (int i = 0; i < num_objects; i++){
                        if (objects[i].exists == true) {
                                float xforce = 0, yforce = 0, zforce = 0;  
                                for (int j = 0; j < num_objects; j++){
                                        if (objects[j].exists == true & i != j){                                              
                                                float x = 0, y = 0, z = 0;  
                                                x = objects[j].px - objects[i].px;
                                                y = objects[j].py - objects[i].py;
                                                z = objects[j].pz - objects[i].pz;
                                                float botPart = pow(sqrt(pow(x,2) + pow(y,2) + pow(z,2)),3);  
                                                xforce += (gConst * objects[i].mass * objects[j].mass * x)/botPart;
                                                yforce += (gConst * objects[i].mass * objects[j].mass * y)/botPart;
                                                zforce += (gConst * objects[i].mass * objects[j].mass * z)/botPart;
                                        }
                                }
                        objects[i].fx = xforce;
                        objects[i].fy = yforce;
                        objects[i].fz = zforce;
                        }
                }
                // Calculo de velocidades, aceleraciones y posiciones
                for (int i = 0; i < num_objects; i++){
                        //Calculo aceleracion, velocidad y posicion de cada objeto  
                        objects[i].ax = objects[i].fx/objects[i].mass;
                        objects[i].ay = objects[i].fy/objects[i].mass;
                        objects[i].az = objects[i].fz/objects[i].mass;
                        objects[i].vx += objects[i].ax * time_step;
                        objects[i].vy += objects[i].ay * time_step;
                        objects[i].vz += objects[i].az * time_step;
                        objects[i].px += objects[i].vx * time_step;
                        objects[i].py += objects[i].vy * time_step;
                        objects[i].pz += objects[i].vz * time_step;
                        if(objects[i].px <= 0){
                                objects[i].px = 1;
                                objects[i].vx = objects[i].vx * -1;
                        }
                        else if(objects[i].px >= size_enclosure){
                                objects[i].px = size_enclosure;
                                objects[i].vx = objects[i].vx * -1;
                        }
                        if(objects[i].py <= 0){
                                objects[i].py = 1;
                                objects[i].vy = objects[i].vy * -1;
                        }
                        else if(objects[i].py >= size_enclosure){
                                objects[i].py = size_enclosure;
                                objects[i].vy = objects[i].vy * -1;
                        }
                        if(objects[i].pz <= 0){
                                objects[i].pz = 1;
                                objects[i].vz = objects[i].vz * -1;
                        }
                        else if(objects[i].pz >= size_enclosure){
                                objects[i].pz = size_enclosure;
                                objects[i].vz = objects[i].vz * -1;
                        }
                }
        }
        generateDocuments("./final_config.txt", objects, size_enclosure, time_step, num_objects);
        struct timeval end;
        gettimeofday(&end, NULL);
        long int endms = end.tv_sec * 1000 + end.tv_usec / 1000;
        long int timeRunning = endms - startms;
        cout << "Execution time; "<<  timeRunning << " milisecond/s" << endl;      
        return 0;
};

