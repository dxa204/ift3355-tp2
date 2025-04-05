#pragma once

#include "linalg/linalg.h"
#include <cmath>  // For sqrt
using namespace linalg::aliases;

#define PI 3.14159265358979323846
#define EPSILON 1e-6

// Valeur aléatoire entre [0,1)
static double rand_double() {
	return double(rand()) / double((RAND_MAX));
}

// Valeur aléatoire entre [0,1] pour un vecteur
static double2 rand_double2() {
	return double2{rand_double(),rand_double()};
}

// Valeur aléatoire à l'intérieur d'un disque.
static double2 random_in_unit_disk() {
    while (true) {
        auto p = (2.0 * rand_double2() - 1.0);
        if (length2(p) >= 1) continue;
        return p;
    }
}

// Convertir radian vers degrée
static double rad2deg(double rad) {
	return rad * 360.0 / (2 * PI);
}

// Convertir degrée vers radian
static double deg2rad(double deg) {
	return deg * 2 * PI / 360.0;
}

// Une classe qui représente un rayon
class Ray 
{
public:
    Ray() : origin(0, 0, 0), direction(0, 0, 0) {}
    Ray(double3 origin_, double3 direction_) :
        origin(origin_), direction(normalize(direction_))
    {
    }

    // Fonction pour calculer la réflexion d'un vecteur incident
    // Prend en entrée :
    // - un vecteur incident 'v' : le vecteur qui frappe la surface
    // - un vecteur normal 'n' : le vecteur perpendiculaire à la surface au point d'impact
    // Retourne le vecteur réfléchi basé sur la loi de la réflexion
    Ray reflect(const double3& position, const double3& normal) const {
        double3 reflected_dir = direction - 2 * dot(direction, normal) * normal;
        return Ray(position, normalize(reflected_dir));
    }

    // Fonction pour calculer la réfraction d'un vecteur incident
    // Prend en entrée :
    // - un vecteur incident 'v'
    // - un vecteur normal 'n'
    // - un indice de réfraction 'ni_over_nt' (indice relatif)
    // Retourne le vecteur réfracté si la réfraction est possible
    // Si la réfraction n'est pas possible (cas de réflexion totale interne), retourne un vecteur nul
    Ray refract(const double3& position, const double3& normal, double eta_ratio) const {
        double cos_theta = dot(-direction, normal);
        double3 r_out_perp = eta_ratio * (direction + cos_theta * normal);
        double3 r_out_parallel = -sqrt(fabs(1.0 - length2(r_out_perp))) * normal;
        double3 refracted_dir = r_out_perp + r_out_parallel;
        return Ray(position, normalize(refracted_dir));
    }

    double3 origin;    // Origine du rayon
    double3 direction; // Direction du rayon
};
