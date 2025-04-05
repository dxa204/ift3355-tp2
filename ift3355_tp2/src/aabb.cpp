#include "aabb.h"
#include <limits>
#include <vector> 

// @@@@@@ VOTRE CODE ICI
// Implémenter l'intersection d'un rayon avec un AABB dans l'intervalle décrit.

// Intersection d'un rayon avec une boîte englobante alignée sur les axes (AABB)
// Cette fonction retourne vrai s'il y a une intersection dans l'intervalle [t_min, t_max]
bool AABB::intersect(Ray ray, double t_min, double t_max) {
    for (int i = 0; i < 3; ++i) {
        // Calculer la direction inverse pour éviter la division
        double invD = 1.0 / ray.direction[i];
        
        // Calculer les distances d'intersection avec les deux plans de l'AABB sur cet axe
        double t0 = (min[i] - ray.origin[i]) * invD;
        double t1 = (max[i] - ray.origin[i]) * invD;

        // Assurer que t0 est le minimum et t1 est le maximum
        if (invD < 0.0) std::swap(t0, t1);

        // Mettre à jour t_min et t_max pour affiner l'intervalle d'intersection
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        // Si les intervalles ne se chevauchent pas, retourner faux car il n'y a pas d'intersection
        if (t_max <= t_min) return false;
    }
    return true;
}

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction qui permet de trouver les 8 coins de notre AABB.
std::vector<double3> retrieve_corners(AABB aabb) {
    // Récupérer les points min et max
    double3 min = aabb.min;
    double3 max = aabb.max;

    // Calculer les 8 coins en combinant les valeurs min et max pour chaque axe
    std::vector<double3> corners = {
        {min.x, min.y, min.z}, {min.x, min.y, max.z}, {min.x, max.y, min.z}, {min.x, max.y, max.z},
        {max.x, min.y, min.z}, {max.x, min.y, max.z}, {max.x, max.y, min.z}, {max.x, max.y, max.z}
    };

    return corners;
}

// @@@@@@ VOTRE CODE ICI
// Implémenter la fonction afin de créer un AABB qui englobe tous les points.

AABB construct_aabb(const std::vector<double3>& points) {
    // Initialiser les points min et max avec des valeurs extrêmes
    double3 min_point = {std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity()};
    double3 max_point = {-std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()};

    // Itérer sur chaque point pour mettre à jour les limites min et max
    for (const auto& point : points) {
        min_point.x = std::min(min_point.x, point.x);
        min_point.y = std::min(min_point.y, point.y);
        min_point.z = std::min(min_point.z, point.z);

        max_point.x = std::max(max_point.x, point.x);
        max_point.y = std::max(max_point.y, point.y);
        max_point.z = std::max(max_point.z, point.z);
    }

    // Construire et retourner l'AABB avec les points min et max calculés
    return AABB{min_point, max_point};
}

AABB combine(AABB a, AABB b) {
    return AABB{min(a.min, b.min), max(a.max, b.max)};
}

bool compare(AABB a, AABB b, int axis) {
    return a.min[axis] < b.min[axis];
}
