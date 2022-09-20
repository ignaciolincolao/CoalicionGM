// coalition.cpp: define el punto de entrada de la aplicación.
//
#include <iostream> //general
#include <nlohmann/json.hpp> //json
#include <string> //string
#include <cstring>
#include <fstream> //ifstream
#include <algorithm>
#include <vector>
#include <chrono>
#include "functions.h"
using namespace std;
using json = nlohmann::json;
double alpha = 0.99;//revisar
double delta = 0.9; //0.95 es el limit
double delta2 = -1;
string votes = "votes.json";
int main(int argc, char* argv[])
{
    if (argc > 1)
    {
        alpha = stoi(argv[1]);
        delta= stoi(argv[2]);
        delta2= stoi(argv[3]);
        votes = stoi(argv[4]);
    }
    //cargar file de votaciones
    ifstream file(votes);
    json data = json::parse(file);

    //se crea y abre el file de salida
    ofstream results;
    results.open("results.json");

    //numero de parlamentario
    int n = data["rollcalls"][0]["votes"].size();
    //creacion de la matriz de distancia
    double** distance_matrix = (double**)malloc(n * sizeof(double*));
    for (size_t i = 0; i < n; i++)
    {
        distance_matrix[i] = (double*)malloc(n * sizeof(double));
    }
    //Creacion matriz de posiciones
    double** position_matrix = (double**)malloc(n * sizeof(double*));
    for (size_t i = 0; i < n; i++)
    {
        position_matrix[i] = (double*)malloc(2 * sizeof(double));
        position_matrix[i][0] = data["rollcalls"][0]["votes"][i]["x"];
        //cout << fixed << position_matrix[i][0] << setprecision(9) << ",";
        position_matrix[i][1] = data["rollcalls"][0]["votes"][i]["y"];
        //break;
    }
    //rellenado de la matriz de distancia
    // Matriz rellenando por completo
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            distance_matrix[i][j] = eucledian_distance(position_matrix[i][0], position_matrix[i][1], position_matrix[j][0], position_matrix[j][1]);
        }
    }

    //inicializacion de quorum
    int quorum = trunc(n / 2) + 1;

    
    //Calculo Tiempo Inicial
    auto initial_time = chrono::high_resolution_clock::now();
    /////////////////////////////////////
    ////// Ordenamiento algoritmo B
    /////////////////////////////////////
    int** congressmen = (int**)malloc(n * sizeof(int*));
    for (size_t i = 0; i < n; i++)
    {
        congressmen[i] = (int*)malloc(quorum * sizeof(int));
    }

    double* initial_fitness = (double*)malloc(n * sizeof(double));
    int* initial_fitness_index = (int*)malloc(n * sizeof(int));

    for (int j = 0; j < n; j++) {
        minDistEdge(congressmen[j], distance_matrix[j], n, quorum);
        sort_bubble(congressmen[j], quorum);
        initial_fitness[j] = evaluate_solution(congressmen[j], distance_matrix, quorum);
    }

    /// Ordena los resultados
    sort_bubble_index(initial_fitness_index, initial_fitness, n, n);
    //sacar el mejor
    int* coalition = (int*)malloc(quorum * sizeof(int));
    double fitness_minimum_winning_coalition;

    memcpy(coalition, congressmen[initial_fitness_index[0]],sizeof(int)*quorum);
    fitness_minimum_winning_coalition = initial_fitness[initial_fitness_index[0]];
    /////////////////////////////////////////////////////////////
    /////// Fin poblacion inicial
    ////////////////////////////////////////////////////////////
    /*for (size_t i = 0; i < quorum; i++)
    {
        cout << coalition[i] << ",";
    }
    cout << endl << fitness << endl;*/
    // Calculo Centroide de Coalicion
    bool improve_posibility = true;
    //inicio de punteros
    double* centroid = (double*)malloc(2 * sizeof(double));
    bool* grid_minimum_winning_coalition = (bool*)malloc(n * sizeof(bool));
    int* not_in_minimum_winning_coalition = (int*)malloc((n - quorum) * sizeof(int));
    int* possible_winning_coalition = (int*)malloc(quorum * sizeof(int));
    //inicio vectores
    struct Distance_vector* distance_vector_minimum_winning_coalition=(Distance_vector*)malloc(sizeof(struct Distance_vector)*(n-quorum));
    vector<Possible_improvement> improvement_vector;
    //inicio variables
    int counter;
    double sum = 0;
    double new_fitness;
    double fitness_copy;
    //double* matDisHull;
    // Cambiar 214 por tamaño Quorum
    struct Point* Pts = (Point*)malloc(sizeof(struct Point) *quorum);
    //Point Pts[quorum];
    
    calcularCentroide(centroid, position_matrix, coalition, quorum);

    for (size_t i = 0; i < quorum; i++)
    {
        Pts[i].x = position_matrix[coalition[i]][0];
        Pts[i].y = position_matrix[coalition[i]][1];
        Pts[i].position = coalition[i];
        Pts[i].index = i;
    }
    //size = sizeof(Pts) / sizeof(Pts[0]);
    //auto hull = convexHull(Pts, size);
    auto hull = convexHull(Pts, quorum);
    for (size_t i = 0; i < n; i++)
    {
        grid_minimum_winning_coalition[i] = 0;
    }
    for (size_t i = 0; i < quorum; i++)
    {
        grid_minimum_winning_coalition[coalition[i]] = 1;
    }
    counter = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (!grid_minimum_winning_coalition[i]) {
            not_in_minimum_winning_coalition[counter] = i;
            counter++;
        }
    }
    //matDisHull = nullptr;
    //matDisHull = (double*)malloc(hull.size() * sizeof(double));
    vector<DisHull> vector_distance_hull;
    vector_distance_hull = distance_of_hull_to_point(hull, centroid);
    distance_of_points_to_coalition(distance_vector_minimum_winning_coalition, not_in_minimum_winning_coalition, coalition, centroid, position_matrix, n, quorum);
    /*
    --1-Calcular cual de los elementos del convex hull estan mas lejos del centroid de la coalition
    --2-Para cada punto que no forme la coalition calcular la sumatoria de las distancias a todos los puntos que si la forman
    --3-Ordenar los que esten mas cerca a los que esten mas lejos
    4-Tomar 1 punto que este mas cerca y probar intercambiando ese punto con el punto mas lejano del centroid y ver si improvement
    5-Si no improvement con ninguno, tomar el segundo punto mas lejano del centroid y repetir proceso
    */

    //double limit = (vector_distance_hull[0].distance) * (alpha);

    double limit;
    fitness_copy = fitness_minimum_winning_coalition;
    int counter_posibility_of_improvement;
    bool improvement = false;
    int number_of_points=0;
    while (improve_posibility)
    {
        delta2 = delta2 + 1;
        limit = vector_distance_hull[0].distance * alpha * pow(delta, delta2);
        for (size_t i = 0; i < (n-quorum); i++) if (distance_vector_minimum_winning_coalition[i].centroid_distance < limit) number_of_points++;
        counter_posibility_of_improvement = 0;
        for (size_t i = 0; i < hull.size(); i++)
        {
            improvement_vector.clear();
            for (size_t j = 0; j < (number_of_points); j++)
            {
                memcpy(possible_winning_coalition, coalition, sizeof(int) * quorum);
                possible_winning_coalition[hull[vector_distance_hull[i].hull_index].index] = distance_vector_minimum_winning_coalition[j].position;
                sort_bubble(possible_winning_coalition, quorum);
                new_fitness = evaluate_solution(possible_winning_coalition, distance_matrix, quorum);
                if (new_fitness < fitness_minimum_winning_coalition)
                {
                    cout << endl;
                    cout << "Cantidad de puntos tomados:" << number_of_points << endl;
                    cout << "Fitness Anterior:" << fixed << fitness_minimum_winning_coalition << setprecision(9) << endl;
                    cout << "Fitness Recalculado:" << fixed << new_fitness << setprecision(9) << endl;
                    cout << "Coalicion Anterior:";
                    for (size_t i = 0; i < quorum; i++)
                    {
                        cout << coalition[i] << ",";
                    }
                    cout << endl;
                    cout << "Coalicion Recalculado:";
                    for (size_t i = 0; i < quorum; i++)
                    {
                        cout << possible_winning_coalition[i] << ",";
                    }
                    cout << endl;
                    improvement_vector.push_back(Possible_improvement());
                    improvement_vector[counter_posibility_of_improvement].fitness = new_fitness;
                    improvement_vector[counter_posibility_of_improvement].index = j;
                    fitness_minimum_winning_coalition = new_fitness;
                    improvement = true;
                }
            }
            sort(improvement_vector.begin(), improvement_vector.end(), &vecMejora_Sort);
            if (improvement)
            {
                coalition[hull[vector_distance_hull[i].hull_index].index] = distance_vector_minimum_winning_coalition[improvement_vector[0].index].position;
                sort_bubble(coalition, quorum);
                fitness_minimum_winning_coalition = evaluate_solution(coalition, distance_matrix, quorum);
                improvement = false;
                break;
            }
        }
        if (fitness_minimum_winning_coalition == fitness_copy)
        {
            improve_posibility = false;
        }
        else {
            fitness_copy = fitness_minimum_winning_coalition;
            //nuevo convex hull
            calcularCentroide(centroid, position_matrix, coalition, quorum);
            Pts = nullptr;
            Pts = (Point*)malloc(sizeof(struct Point) * quorum);
            for (size_t i = 0; i < quorum; i++)
            {
                Pts[i].x = position_matrix[coalition[i]][0];
                Pts[i].y = position_matrix[coalition[i]][1];
                Pts[i].position = coalition[i];
                Pts[i].index = i;
            }
            //size = sizeof(Pts) / sizeof(Pts[0]);
            //auto hull = convexHull(Pts, size);
            auto hull = convexHull(Pts, quorum);
            for (size_t i = 0; i < n; i++)
            {
                grid_minimum_winning_coalition[i] = 0;
            }
            for (size_t i = 0; i < quorum; i++)
            {
                grid_minimum_winning_coalition[coalition[i]] = 1;
            }
            counter = 0;
            for (size_t i = 0; i < n; i++)
            {
                if (!grid_minimum_winning_coalition[i]) {
                    not_in_minimum_winning_coalition[counter] = i;
                    counter++;
                }
            }
            vector_distance_hull.clear();
            vector_distance_hull = distance_of_hull_to_point(hull, centroid);
            distance_vector_minimum_winning_coalition = nullptr;
            distance_vector_minimum_winning_coalition = (Distance_vector*)malloc(sizeof(struct Distance_vector) * (n - quorum));
            distance_of_points_to_coalition(distance_vector_minimum_winning_coalition, not_in_minimum_winning_coalition, coalition, centroid, position_matrix, n, quorum);
        }
    }
    auto final_time = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(final_time - initial_time).count();
    time_taken *= 1e-9;

    cout << "Algoritmo terminado - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    cout << "Tiempo total:" << fixed << time_taken << setprecision(9) << endl;
    cout << "Fitness Final:" << fixed << fitness_minimum_winning_coalition << setprecision(9) << endl;
    cout << "Coalicion:";
    for (size_t i = 0; i < quorum; i++)
    {
        cout << coalition[i] << ",";
    }
    cout << endl;
}
