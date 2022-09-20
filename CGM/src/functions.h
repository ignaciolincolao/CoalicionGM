#include <vector>
using namespace std;
double minimum_distance(double d[], bool congressman_set[], int n)
{
    double min = FLT_MAX, min_index;
    for (int v = 0; v < n; v++)
        if (congressman_set[v] == false && d[v] < min)
            min = d[v], min_index = v;
    return min_index;
}

double evaluate_solution(int* position, double** mat, int length) {
    double sum = 0;
    for (size_t i = 0; i <= (length - 2); i++)
    {

        for (size_t j = i + 1; j <= (length - 1); j++)
        {
            sum = sum + mat[position[i]][position[j]];
        }
    }
    return sum;
}

//funcion para ordenar un arreglo de menor a mayor
void sort_bubble(int* array, int length)
{
    int temp = 0;
    for (size_t i = 0; i < length; i++)
    {
        bool already_sorted = true;
        for (size_t j = 0; j < length - i - 1; j++)
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

void sort_bubble_index(int* index_array, double* arrayDist, int n, int quorum) {
    int temp = 0;
    //int order_index[n];
    int* order_index = new int[n];
    for (size_t i = 0; i < n; i++)
        order_index[i] = i;

    for (size_t i = 0; i < n; i++)
    {
        bool already_sorted = true;
        for (size_t j = 0; j < n - i - 1; j++)
        {
            if (arrayDist[order_index[j]] > arrayDist[order_index[j + 1]])
            {
                temp = order_index[j];
                order_index[j] = order_index[j + 1];
                order_index[j + 1] = temp;
                already_sorted = false;
            }
        }
        if (already_sorted)
            break;
    }
    memcpy(index_array, order_index, sizeof(int) * quorum);
}

void minDistEdge(int* index_array, double* arrayDist, int n, int quorum)
{
    bool* congressman_set = new bool[n];
    for (size_t i = 0; i < n; i++) {
        congressman_set[i] = false;
    }
    for (size_t i = 0; i < quorum; i++) {
        int u = minimum_distance(arrayDist, congressman_set, n);
        congressman_set[u] = true;
        index_array[i] = u;
    }
}

double eucledian_distance(double x1, double y1, double x2, double y2)
{
    double calculation = pow(pow((x2 - x1), 2) + pow((y2 - y1), 2), 1 / (double)2);
    return calculation;
}
struct Point
{
    double x, y;
    int position;
    int index;
};
struct Distance_vector {
    double distance;
    double centroid_distance;
    int position;
};
struct Possible_improvement {
    double fitness;
    int index;
};
struct DisHull {
    double distance;
    int hull_index;
};
double orientation(Point p, Point q, Point r)
{
    double value = (q.y - p.y) * (r.x - q.x) -
        (q.x - p.x) * (r.y - q.y);

    if (value == 0) return 0;  // collinear
    return (value > 0) ? 1 : 2; // clock or counterclock wise
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
bool VectDist_Sort(Distance_vector const& lvd, Distance_vector const& rvd) {
    return lvd.distance < rvd.distance;
}

bool vecMejora_Sort(Possible_improvement const& lvd, Possible_improvement const& rvd) {
    return lvd.fitness < rvd.fitness;
}
bool vector_distance_hull_sort(DisHull const& lvd, DisHull const& rvd) {
    return lvd.distance > rvd.distance;
}
void calcularCentroide(double* centroide, double** matPos, int* coalicion, int quorum) {
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

//Sort
void swap_double(Distance_vector* a, Distance_vector* b) {
    Distance_vector temp = *a;
    *a = *b;
    *b = temp;
}

void selection_sort_distance_vector(Distance_vector array[], int size) {
    for (int step = 0; step < size - 1; step++) {
        int index_of_minimum = step;
        for (int i = step + 1; i < size; i++) if (array[i].distance < array[index_of_minimum].distance) index_of_minimum = i;
        swap_double(&array[index_of_minimum], &array[step]);
    }
}

vector<DisHull> distance_of_hull_to_point(vector<Point> hull,double* centroid) {
    vector<DisHull> vector_distance_hull;
    for (size_t i = 0; i < hull.size(); i++) {
        vector_distance_hull.push_back(DisHull());
        vector_distance_hull[i].distance = eucledian_distance(hull[i].x, hull[i].y, centroid[0], centroid[1]);
        vector_distance_hull[i].hull_index = i;
    }
    sort(vector_distance_hull.begin(), vector_distance_hull.end(), &vector_distance_hull_sort);
    return vector_distance_hull;
}
void distance_of_points_to_coalition(struct Distance_vector* distance_vector_minimum_winning_coalition,int* not_in_minimum_winning_coalition, int* coalition, double* centroid, double** position_matrix,int n, int quorum) {
    double sum = 0;
    for (int i = 0; i < (n - quorum); i++) {
        for (size_t j = 0; j < quorum; j++)
        {
            sum = sum + eucledian_distance(position_matrix[not_in_minimum_winning_coalition[i]][0], position_matrix[not_in_minimum_winning_coalition[i]][1], position_matrix[coalition[j]][0], position_matrix[coalition[j]][1]);
        }
        distance_vector_minimum_winning_coalition[i].distance = sum;
        distance_vector_minimum_winning_coalition[i].position = not_in_minimum_winning_coalition[i];
        distance_vector_minimum_winning_coalition[i].centroid_distance = eucledian_distance(position_matrix[not_in_minimum_winning_coalition[i]][0], position_matrix[not_in_minimum_winning_coalition[i]][1], centroid[0], centroid[1]);
        sum = 0;
    }
    selection_sort_distance_vector(distance_vector_minimum_winning_coalition, (n - quorum));
}