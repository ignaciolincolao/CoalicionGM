// coalicion.cpp: define el punto de entrada de la aplicación.
//
#include <iostream> //general
#include <nlohmann/json.hpp> //json
#include <string> //string
#include <cstring>
#include <fstream> //ifstream
#include <algorithm>
#include <vector>
using namespace std;
using json = nlohmann::json;
float minDist(float d[], bool genSet[], int n)
{
    float min = FLT_MAX, min_index;
    for (int v = 0; v < n; v++)
        if (genSet[v] == false && d[v] < min)
            min = d[v], min_index = v;
    return min_index;
}

float eval_sol(int* pos, float** mat, int largo) {
    float suma = 0;
    for (size_t i = 0; i <= (largo - 2); i++)
    {

        for (size_t j = i + 1; j <= (largo - 1); j++)
        {
            suma = suma + mat[pos[i]][pos[j]];
        }
    }
    return suma;
}

//funcion para ordenar un arreglo de menor a mayor
void sort_bubble(int* array, int largo)
{
    int temp = 0;
    for (size_t i = 0; i < largo; i++)
    {
        bool already_sorted = true;
        for (size_t j = 0; j < largo - i - 1; j++)
        {
            if (array[j] > array[j + 1])
            {
                temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
                already_sorted = false;
            }
        }
        if (already_sorted)
            break;
    }
}

void sort_BubbleIndex(int* arrayIndex, float* arrayDist, int n, int quorum) {
    int temp = 0;
    //int orderIndex[n];
    int* orderIndex = new int[n];
    for (size_t i = 0; i < n; i++)
        orderIndex[i] = i;

    for (size_t i = 0; i < n; i++)
    {
        bool already_sorted = true;
        for (size_t j = 0; j < n - i - 1; j++)
        {
            if (arrayDist[orderIndex[j]] > arrayDist[orderIndex[j + 1]])
            {
                temp = orderIndex[j];
                orderIndex[j] = orderIndex[j + 1];
                orderIndex[j + 1] = temp;
                already_sorted = false;
            }
        }
        if (already_sorted)
            break;
    }
    memcpy(arrayIndex, orderIndex, sizeof(int) * quorum);
}

void minDistEdge(int* arrayIndex, float* arrayDist, int n, int quorum)
{
    //bool genSet[n]; 
    bool* genSet = new bool[n];
    for (size_t i = 0; i < n; i++) {
        genSet[i] = false;
    }
    for (size_t i = 0; i < quorum; i++) {
        int u = minDist(arrayDist, genSet, n);
        genSet[u] = true;
        arrayIndex[i] = u;
    }
}

float dis_euc(float x1, float y1, float x2, float y2)
{
    float calculo = pow(pow((x2 - x1), 2) + pow((y2 - y1), 2), 1 / (float)2);
    return calculo;
}
struct Point
{
    float x, y;
    int pos;
    int indice;
};
struct VectDis {
    float Distancia;
    int pos;
};
struct pMejora {
    float fitness;
    int indice;
};
struct DisHull {
    float Distancia;
    int indiceHull;
};
float orientation(Point p, Point q, Point r)
{
    float val = (q.y - p.y) * (r.x - q.x) -
        (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // collinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// Prints convex hull of a set of n points.
vector<Point> convexHull(Point points[], int n)
{
    // There must be at least 3 points
    // Initialize Result
    vector<Point> hull;
    if (n < 3) return hull;
    // Find the leftmost point
    int l = 0;
    for (int i = 1; i < n; i++)
        if (points[i].x < points[l].x)
            l = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int p = l, q;
    do
    {
        // Add current point to result
        hull.push_back(points[p]);

        // Search for a point 'q' such that orientation(p, q,
        // x) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        q = (p + 1) % n;
        for (int i = 0; i < n; i++)
        {
            // If i is more counterclockwise than current q, then
            // update q
            if (orientation(points[p], points[i], points[q]) == 2)
                q = i;
        }

        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;

    } while (p != l);  // While we don't come to first point
    // Print Result
    return hull;
}
bool VectDist_Sort(VectDis const& lvd, VectDis const& rvd) {
    return lvd.Distancia < rvd.Distancia;
}

bool vecMejora_Sort(pMejora const& lvd, pMejora const& rvd) {
    return lvd.fitness < rvd.fitness;
}
bool vecDisHull_Sort(DisHull const& lvd, DisHull const& rvd) {
    return lvd.Distancia > rvd.Distancia;
}
void calcularCentroide(float* centroide,float** matPos,int* coalicion,int quorum) {
    centroide[0] = 0;
    centroide[1] = 0;
    for (size_t i = 0; i < quorum; i++)
    {
        centroide[0] = centroide[0] + matPos[coalicion[i]][0];
        centroide[1] = centroide[1] + matPos[coalicion[i]][1];
    }
    centroide[0] = centroide[0] / quorum;
    centroide[1] = centroide[1] / quorum;
}
// ConvexHull obtenido de https://www.geeksforgeeks.org/convex-hull-set-1-jarviss-algorithm-or-wrapping/
int main()
{
    //cargar archivo de votaciones
    ifstream archivo("votacion.json");
    json data = json::parse(archivo);

    //se crea y abre el archivo de salida
    ofstream resultados;
    resultados.open("resultados.json");

    //numero de parlamentario
    int n = data["rollcalls"][0]["votes"].size();

    //creacion de la matriz de distancia
    float** matDis = (float**)malloc(n * sizeof(float*));
    for (size_t i = 0; i < n; i++)
    {
        matDis[i] = (float*)malloc(n * sizeof(float));
    }
    //Creacion matriz de posiciones
    float** matPos = (float**)malloc(n * sizeof(float*));
    for (size_t i = 0; i < n; i++)
    {
        matPos[i] = (float*)malloc(2 * sizeof(float));
        matPos[i][0] = data["rollcalls"][0]["votes"][i]["x"];
        matPos[i][1] = data["rollcalls"][0]["votes"][i]["y"];
    }
    //rellenado de la matriz de distancia
    // Matriz rellenando por completo
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            matDis[i][j] = dis_euc(matPos[i][0], matPos[i][1], matPos[j][0], matPos[j][1]);
        }
    }
    //inicializacion de quorum
    int quorum = trunc(n / 2) + 1;

    

    /////////////////////////////////////
    ////// Ordenamiento algoritmo B
    /////////////////////////////////////
    int** congresistas = (int**)malloc(n * sizeof(int*));
    for (size_t i = 0; i < n; i++)
    {
        congresistas[i] = (int*)malloc(quorum * sizeof(int));
    }

    float* fitnessInit = (float*)malloc(n * sizeof(float));
    int* fitnessInitIndex = (int*)malloc(n * sizeof(int));

    for (int j = 0; j < n; j++) {
        minDistEdge(congresistas[j], matDis[j], n, quorum);
        sort_bubble(congresistas[j], quorum);
        fitnessInit[j] = eval_sol(congresistas[j], matDis, quorum);
    }

    /// Ordena los resultados
    sort_BubbleIndex(fitnessInitIndex, fitnessInit, n, n);
    //sacar el mejor
    int* coalicion = (int*)malloc(quorum * sizeof(int));
    float fitnessCGM;

    memcpy(coalicion, congresistas[fitnessInitIndex[0]],sizeof(int)*quorum);
    fitnessCGM = fitnessInit[fitnessInitIndex[0]];
    /////////////////////////////////////////////////////////////
    /////// Fin poblacion inicial
    ////////////////////////////////////////////////////////////
    /*for (size_t i = 0; i < quorum; i++)
    {
        cout << coalicion[i] << ",";
    }
    cout << endl << fitness << endl;*/
    // Calculo Centroide de Coalicion
    bool posibilidadMejora = true;
    //inicio de punteros
    float* centroide = (float*)malloc(2 * sizeof(float));
    float* vecDis = (float*)malloc(n * sizeof(float));
    bool* mallaCGM = (bool*)malloc(n * sizeof(bool));
    int* notCGM = (int*)malloc((n - quorum) * sizeof(int));
    int* CGP = (int*)malloc(quorum * sizeof(int));
    //inicio vectores
    vector<VectDis> vectDisCGM;
    vector<pMejora> vecMejora;
    //inicio variables
    int size;
    int cont;
    float maxDis;
    int maxDisPos;
    float sum = 0;
    int cantidad = 50;
    float fitnessNuevo;
    float copiaFit;
    //float* matDisHull;
    // Cambiar 214 por tamaño Quorum
    Point Pts[214];
    
    calcularCentroide(centroide, matPos, coalicion, quorum);
    for (size_t i = 0; i < quorum; i++)
    {
        cout << coalicion[i] << ",";
    }
    cout << endl << fitnessCGM << endl;

    for (size_t i = 0; i < quorum; i++)
    {
        Pts[i].x = matPos[coalicion[i]][0];
        Pts[i].y = matPos[coalicion[i]][1];
        Pts[i].pos = coalicion[i];
        Pts[i].indice = i;
    }
    size = sizeof(Pts) / sizeof(Pts[0]);
    auto hull = convexHull(Pts, size);
    for (size_t i = 0; i < n; i++)
    {
        mallaCGM[i] = 0;
    }
    for (size_t i = 0; i < quorum; i++)
    {
        mallaCGM[coalicion[i]] = 1;
    }
    cont = 0;
    for (size_t i = 0; i < n; i++)
    {
        if (!mallaCGM[i]) {
            notCGM[cont] = i;
            cont++;
        }
    }
    //matDisHull = nullptr;
    //matDisHull = (float*)malloc(hull.size() * sizeof(float));
    vector<DisHull> vecDisHull;
    for (size_t i = 0; i < hull.size(); i++) {
        //matDisHull[i] = dis_euc(hull[i].x, hull[i].y, centroide[0], centroide[1]);
        vecDisHull.push_back(DisHull());
        vecDisHull[i].Distancia = dis_euc(hull[i].x, hull[i].y, centroide[0], centroide[1]);
        vecDisHull[i].indiceHull = i;
    }
    sort(vecDisHull.begin(), vecDisHull.end(), &vecDisHull_Sort);
    vectDisCGM.clear();
    for (int i = 0; i < (n - quorum); i++) {
        for (size_t j = 0; j < quorum; j++)
        {
            sum = sum + dis_euc(matPos[notCGM[i]][0], matPos[notCGM[i]][1], matPos[coalicion[j]][0], matPos[coalicion[j]][1]);
            //sum = sum + matDis[notCGM[i]][coalicion[j]];
        }
        vectDisCGM.push_back(VectDis());
        vectDisCGM[i].Distancia = sum;
        vectDisCGM[i].pos = notCGM[i];
        sum = 0;
    }

    sort(vectDisCGM.begin(), vectDisCGM.end(), &VectDist_Sort);
    /*
    --1-Calcular cual de los elementos del convex hull estan mas lejos del centroide de la coalicion
    --2-Para cada punto que no forme la coalicion calcular la sumatoria de las distancias a todos los puntos que si la forman
    --3-Ordenar los que esten mas cerca a los que esten mas lejos
    4-Tomar 1 punto que este mas cerca y probar intercambiando ese punto con el punto mas lejano del centroide y ver si mejora
    5-Si no mejora con ninguno, tomar el segundo punto mas lejano del centroide y repetir proceso
    */
    copiaFit = fitnessCGM;
    int contPM;
    bool mejora = false;
    while (posibilidadMejora)
    {
        contPM = 0;
        for (size_t i = 0; i < hull.size(); i++)
        {
            vecMejora.clear();
            for (size_t j = 0; j < (vectDisCGM.size() - cantidad); j++)
            {
                memcpy(CGP, coalicion, sizeof(int) * quorum);
                CGP[hull[vecDisHull[i].indiceHull].indice] = vectDisCGM[j].pos;
                sort_bubble(CGP, quorum);
                fitnessNuevo = eval_sol(CGP, matDis, quorum);
                if (fitnessNuevo < fitnessCGM)
                {
                    vecMejora.push_back(pMejora());
                    vecMejora[contPM].fitness = fitnessNuevo;
                    vecMejora[contPM].indice = j;
                    cout << endl << fitnessNuevo << endl;
                    fitnessCGM = fitnessNuevo;
                    mejora = true;
                }
            }
            sort(vecMejora.begin(), vecMejora.end(), &vecMejora_Sort);
            if (mejora)
            {
                coalicion[hull[vecDisHull[i].indiceHull].indice] = vectDisCGM[vecMejora[0].indice].pos;
                sort_bubble(coalicion, quorum);
                fitnessCGM = eval_sol(coalicion, matDis, quorum);
                mejora = false;
                break;
            }
        }
        if (fitnessCGM == copiaFit)
        {
            posibilidadMejora = false;
        }
        else {
            copiaFit = fitnessCGM;
            //nuevo convex hull
            calcularCentroide(centroide, matPos, coalicion, quorum);
            for (size_t i = 0; i < quorum; i++)
            {
                Pts[i].x = matPos[coalicion[i]][0];
                Pts[i].y = matPos[coalicion[i]][1];
                Pts[i].pos = coalicion[i];
                Pts[i].indice = i;
            }
            size = sizeof(Pts) / sizeof(Pts[0]);
            auto hull = convexHull(Pts, size);
            for (size_t i = 0; i < n; i++)
            {
                mallaCGM[i] = 0;
            }
            for (size_t i = 0; i < quorum; i++)
            {
                mallaCGM[coalicion[i]] = 1;
            }
            cont = 0;
            for (size_t i = 0; i < n; i++)
            {
                if (!mallaCGM[i]) {
                    notCGM[cont] = i;
                    cont++;
                }
            }
            vecDisHull.clear();
            for (size_t i = 0; i < hull.size(); i++) {
                //matDisHull[i] = dis_euc(hull[i].x, hull[i].y, centroide[0], centroide[1]);
                vecDisHull.push_back(DisHull());
                vecDisHull[i].Distancia = dis_euc(hull[i].x, hull[i].y, centroide[0], centroide[1]);
                vecDisHull[i].indiceHull = i;
            }
            sort(vecDisHull.begin(), vecDisHull.end(), &vecDisHull_Sort);
            vectDisCGM.clear();
            for (int i = 0; i < (n - quorum); i++) {
                for (size_t j = 0; j < quorum; j++)
                {
                    //sum = sum + matDis[notCGM[i]][coalicion[j]];
                    sum = sum + dis_euc(matPos[notCGM[i]][0], matPos[notCGM[i]][1], matPos[coalicion[j]][0], matPos[coalicion[j]][1]);
                }
                vectDisCGM.push_back(VectDis());
                vectDisCGM[i].Distancia = sum;
                vectDisCGM[i].pos = notCGM[i];
                sum = 0;
            }
            sort(vectDisCGM.begin(), vectDisCGM.end(), &VectDist_Sort);
        }
    }
    cout << "Algoritmo terminado" << endl;
    //cout << "Fitness Final:" << fitnessCGM<<endl;
    cout << "Coalicion:";
    for (size_t i = 0; i < quorum; i++)
    {
        cout << coalicion[i] << ",";
    }
    cout << endl;
    cout << "Hull:" << endl;
    for (auto punto : hull) {
        cout << punto.pos << ",";
    }
}
