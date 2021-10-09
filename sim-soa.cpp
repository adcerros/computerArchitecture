#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <iomanip>

using namespace std;
//Definicion del objeto
struct object {
        float x, y, z, vx, vy, vz, mass;
};

int main (int argc, char * argv[]){
        
        // Comprobacion del numero de entradas
        int numberOfParams = argc - 1;
        if (numberOfParams!= 5){
                cerr << "sim-soa invoked with " + to_string(numberOfParams) + " parameters.\n";
                return -1;      
        };
        // Comprobacion de los datos
        int num_objects = stoi(argv[1]);
        int num_interations = stoi(argv[2]);
        int random_seed = stoi(argv[3]);
        float size_enclosure = stof(argv[4]);
        float time_step = stof(argv[5]);
        if (num_objects < 0 | num_interations < 0 | random_seed <= 0 | size_enclosure <= 0 | time_step <= 0){
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
                current.x = dis_uniform(generator);
                current.y = dis_uniform(generator);
                current.z = dis_uniform(generator);
                current.mass = dis_normal(generator);
                current.vx = dis_uniform(generator);
                current.vy = dis_uniform(generator);
                current.vz = dis_uniform(generator);
                objects[i] = current;
        }
        // Se generan los documentos
        ofstream file;
        file.open("./init_config.txt", ofstream::out | ofstream::trunc);
        file  << size_enclosure << " " << time_step << " " << num_objects <<"\n";
        for(int i = 0; i < num_objects; i++){
                object current = objects[i];
                file << current.x << " " << current.y <<  " " << current.z 
                <<  " " << current.vx <<  " " << current.vy <<  " " << current.vz 
                <<  " " << current.mass <<  "\n";
        };
        file.close();
        return 0;
};

