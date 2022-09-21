#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include "functions.h"

using namespace std;
using json = nlohmann::json;

// Default params
double alpha = 0.99;
double delta = 0.9; // 0.95 its a limit
double delta2 = -1;
string votes = "votes.json";

int main(int argc, char *argv[])
{
    // Params
    if (argc > 1)
    {
        alpha = stoi(argv[1]);
        delta = stoi(argv[2]);
        delta2 = stoi(argv[3]);
        votes = stoi(argv[4]);
    }
    // Load voting file
    ifstream file(votes);
    json data = json::parse(file);

    // Create output file
    ofstream results;
    results.open("results.json");

    // Number of congressmen
    int n = data["rollcalls"][0]["votes"].size();
    // Creating the distance matrix
    double **distance_matrix = (double **)malloc(n * sizeof(double *));
    for (size_t i = 0; i < n; i++)
    {
        distance_matrix[i] = (double *)malloc(n * sizeof(double));
    }
    // Creation and filling of position matrix
    double **position_matrix = (double **)malloc(n * sizeof(double *));
    for (size_t i = 0; i < n; i++)
    {
        position_matrix[i] = (double *)malloc(2 * sizeof(double));
        position_matrix[i][0] = data["rollcalls"][0]["votes"][i]["x"];
        position_matrix[i][1] = data["rollcalls"][0]["votes"][i]["y"];
    }
    // Filling the distance matrix
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            distance_matrix[i][j] = eucledian_distance(position_matrix[i][0], position_matrix[i][1], position_matrix[j][0], position_matrix[j][1]);
        }
    }

    // Quorum initialization
    int quorum = trunc(n / 2) + 1;

    // Calculate start time
    auto initial_time = chrono::high_resolution_clock::now();
    // Congressmen matrix Initialization
    int **congressmen = (int **)malloc(n * sizeof(int *));
    for (size_t i = 0; i < n; i++)
    {
        congressmen[i] = (int *)malloc(quorum * sizeof(int));
    }
    // Fitness initialization
    double *initial_fitness = (double *)malloc(n * sizeof(double));
    int *initial_fitness_index = (int *)malloc(n * sizeof(int));
    // Creation of initial solutions
    for (int j = 0; j < n; j++)
    {
        minDistEdge(congressmen[j], distance_matrix[j], n, quorum);
        sort_bubble(congressmen[j], quorum);
        initial_fitness[j] = evaluate_solution(congressmen[j], distance_matrix, quorum);
    }
    // Sort the solutions
    sort_bubble_index(initial_fitness_index, initial_fitness, n, n);
    // Create the variables for the best solution
    int *coalition = (int *)malloc(quorum * sizeof(int));
    double fitness_minimum_winning_coalition;
    // bring out the best
    memcpy(coalition, congressmen[initial_fitness_index[0]], sizeof(int) * quorum);
    fitness_minimum_winning_coalition = initial_fitness[initial_fitness_index[0]];
    // Pointer initialization
    double *centroid = (double *)malloc(2 * sizeof(double));
    bool *grid_minimum_winning_coalition = (bool *)malloc(n * sizeof(bool));
    int *not_in_minimum_winning_coalition = (int *)malloc((n - quorum) * sizeof(int));
    int *possible_winning_coalition = (int *)malloc(quorum * sizeof(int));
    // Vector Initialization
    struct Distance_vector *distance_vector_minimum_winning_coalition = (Distance_vector *)malloc(sizeof(struct Distance_vector) * (n - quorum));
    vector<Possible_improvement> improvement_vector;
    // Variable initialization
    bool possibility_of_improvement = true;
    int counter;
    double sum = 0;
    double new_fitness;
    double fitness_copy;
    double limit;
    fitness_copy = fitness_minimum_winning_coalition;
    int counter_posibility_of_improvement;
    bool improvement = false;
    int number_of_points = 0;
    // Point Structure Initialization
    struct Point *Pts = (Point *)malloc(sizeof(struct Point) * quorum);
    // Calculate the centroid of the best solution
    calculate_centroid(centroid, position_matrix, coalition, quorum);
    // save the points of the best solution
    for (size_t i = 0; i < quorum; i++)
    {
        Pts[i].x = position_matrix[coalition[i]][0];
        Pts[i].y = position_matrix[coalition[i]][1];
        Pts[i].position = coalition[i];
        Pts[i].index = i;
    }
    // Calculate the convex hull of the best solution
    auto hull = convexHull(Pts, quorum);
    // Calculate the grid of the voting
    for (size_t i = 0; i < n; i++)
    {
        grid_minimum_winning_coalition[i] = 0;
    }
    for (size_t i = 0; i < quorum; i++)
    {
        grid_minimum_winning_coalition[coalition[i]] = 1;
    }
    // Calculate number of points not in the best solution
    counter = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (!grid_minimum_winning_coalition[i])
        {
            not_in_minimum_winning_coalition[counter] = i;
            counter++;
        }
    }
    // Calculate the distance vector of the points not in the best solution towards the convex hull
    vector<DisHull> vector_distance_hull;
    vector_distance_hull = distance_of_hull_to_point(hull, centroid);
    distance_of_points_to_coalition(distance_vector_minimum_winning_coalition, not_in_minimum_winning_coalition, coalition, centroid, position_matrix, n, quorum);
    while (possibility_of_improvement)
    {
        delta2 = delta2 + 1;
        // Calculate the limit
        limit = vector_distance_hull[0].distance * alpha * pow(delta, delta2);
        // Calculate the number of points that are within the limit
        for (size_t i = 0; i < (n - quorum); i++)
            if (distance_vector_minimum_winning_coalition[i].centroid_distance < limit)
                number_of_points++;
        // Reset the improvement counter
        counter_posibility_of_improvement = 0;

        for (size_t i = 0; i < hull.size(); i++)
        {
            // Reset the improvement vector
            improvement_vector.clear();
            for (size_t j = 0; j < (number_of_points); j++)
            {
                // Obtain the possible winning coalition
                memcpy(possible_winning_coalition, coalition, sizeof(int) * quorum);
                possible_winning_coalition[hull[vector_distance_hull[i].hull_index].index] = distance_vector_minimum_winning_coalition[j].position;
                // sort the possible winning coalition
                sort_bubble(possible_winning_coalition, quorum);
                // Calculate the fitness of the possible winning coalition
                new_fitness = evaluate_solution(possible_winning_coalition, distance_matrix, quorum);
                if (new_fitness < fitness_minimum_winning_coalition)
                {
                    // prints
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
                    // Save the improvement
                    improvement_vector.push_back(Possible_improvement());
                    // Save the fitness of the possible winning coalition
                    improvement_vector[counter_posibility_of_improvement].fitness = new_fitness;
                    // Save the position of the possible winning coalition
                    improvement_vector[counter_posibility_of_improvement].index = j;
                    fitness_minimum_winning_coalition = new_fitness;
                    improvement = true;
                }
            }
            // sort the improvement vector
            sort(improvement_vector.begin(), improvement_vector.end(), &vecMejora_Sort);
            if (improvement)
            {
                // A congressman from the convex hull is changed for one that is within the vector of possible improvement
                coalition[hull[vector_distance_hull[i].hull_index].index] = distance_vector_minimum_winning_coalition[improvement_vector[0].index].position;
                // sort the coalition
                sort_bubble(coalition, quorum);
                // Calculate the new fitness
                fitness_minimum_winning_coalition = evaluate_solution(coalition, distance_matrix, quorum);
                // We get out of the loop
                improvement = false;
                break;
            }
        }
        // if the fitness does not change, we get out of the loop
        if (fitness_minimum_winning_coalition == fitness_copy)
        {
            possibility_of_improvement = false;
        }
        else
        {
            // save the fitness
            fitness_copy = fitness_minimum_winning_coalition;
            // Calculate the centroid
            calculate_centroid(centroid, position_matrix, coalition, quorum);
            // Initialize the points again
            Pts = nullptr;
            Pts = (Point *)malloc(sizeof(struct Point) * quorum);
            for (size_t i = 0; i < quorum; i++)
            {
                Pts[i].x = position_matrix[coalition[i]][0];
                Pts[i].y = position_matrix[coalition[i]][1];
                Pts[i].position = coalition[i];
                Pts[i].index = i;
            }
            // Calculate the convex hull
            auto hull = convexHull(Pts, quorum);
            // Calculate the grid of the voting
            for (size_t i = 0; i < n; i++)
            {
                grid_minimum_winning_coalition[i] = 0;
            }
            for (size_t i = 0; i < quorum; i++)
            {
                grid_minimum_winning_coalition[coalition[i]] = 1;
            }
            // Calculate number of points not in the best solution
            counter = 0;
            for (size_t i = 0; i < n; i++)
            {
                if (!grid_minimum_winning_coalition[i])
                {
                    not_in_minimum_winning_coalition[counter] = i;
                    counter++;
                }
            }
            // Clear the distance vector and calculate again
            vector_distance_hull.clear();
            vector_distance_hull = distance_of_hull_to_point(hull, centroid);
            distance_vector_minimum_winning_coalition = nullptr;
            distance_vector_minimum_winning_coalition = (Distance_vector *)malloc(sizeof(struct Distance_vector) * (n - quorum));
            distance_of_points_to_coalition(distance_vector_minimum_winning_coalition, not_in_minimum_winning_coalition, coalition, centroid, position_matrix, n, quorum);
        }
    }
    // Calculate the time of the algorithm
    auto final_time = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(final_time - initial_time).count();
    time_taken *= 1e-9;
    // Print the results
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
