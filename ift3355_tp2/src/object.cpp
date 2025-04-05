#include "object.h"

// Fonction retournant soit la valeur v0 ou v1 selon le signe.
int rsign(double value, double v0, double v1) {
    return (int(std::signbit(value)) * (v1-v0)) + v0;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection d'une sphère.
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Sphere::local_intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    // La sphère est centrée à l'origine avec un rayon `radius`
    double3 oc = ray.origin;  // vecteur du centre de la sphère (origine) vers l'origine du rayon

    // Coefficients pour l'équation quadratique at^2 + bt + c = 0
    double a = dot(ray.direction, ray.direction);
    double b = 2.0 * dot(oc, ray.direction);
    double c = dot(oc, oc) - radius * radius;

    // Calculer le discriminant
    double discriminant = b * b - 4 * a * c;

    // Si le discriminant est inférieur à zéro, il n'y a pas d'intersection
    if (discriminant < 0) {
        return false;
    }

    // Calculer les deux solutions possibles pour t
    double sqrt_discriminant = sqrt(discriminant);
    double t1 = (-b - sqrt_discriminant) / (2.0 * a);
    double t2 = (-b + sqrt_discriminant) / (2.0 * a);

    // Trouver le t valide le plus proche dans l'intervalle [t_min, t_max]
    double t = t1;
    if (t < t_min || t > t_max) {
        t = t2;
        if (t < t_min || t > t_max) {
            return false;  // Pas d'intersection valide dans l'intervalle
        }
    }

    // Calculer le point d'intersection et la normale
    double3 intersection_point = ray.origin + t * ray.direction;
    double3 normal = normalize(intersection_point);  // Normale de la sphère au point d'intersection

    // Mettre à jour l'information de l'intersection
    hit->depth = t;
    hit->position = intersection_point;
    hit->normal = dot(ray.direction, normal) > 0 ? -normal : normal;  // Assurer que la normale pointe dans la direction opposée au rayon
    hit->uv = {
        0.5 + atan2(normal.z, normal.x) / (2 * PI),
        0.5 - asin(normal.y) / PI
    };

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour la sphère.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Sphere::compute_aabb() {
    // La sphère est centrée à l'origine avec un rayon `radius`.
    double3 min_point = {-radius, -radius, -radius};
    double3 max_point = { radius,  radius,  radius};
    
    // Les points min et max peuvent être transformés en espace global si nécessaire.
    return AABB{min_point, max_point};
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un quad (rectangle).
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Quad::local_intersect(Ray const ray, double t_min, double t_max, Intersection* hit) {
    // Le quad est situé sur le plan Z = 0 avec une normale dirigée vers +Z
    double3 normal = {0, 0, 1};  // Normale du plan du quad
    double denom = dot(ray.direction, normal);

    // Vérifier si le rayon est parallèle au plan du quad
    if (fabs(denom) < EPSILON) {
        return false;  // Pas d'intersection, car le rayon est parallèle
    }

    // Calculer la distance d'intersection t
    double t = -dot(ray.origin, normal) / denom;
    
    // Vérifier si l'intersection est dans l'intervalle valide
    if (t < t_min || t > t_max) {
        return false;
    }

    // Calculer le point d'intersection
    double3 intersection_point = ray.origin + t * ray.direction;

    // Vérifier si le point d'intersection est dans les limites du quad
    if (fabs(intersection_point.x) > half_size || fabs(intersection_point.y) > half_size) {
        return false;  // En dehors des limites du quad
    }

    // Mettre à jour les informations d'intersection
    hit->depth = t;
    hit->position = intersection_point;
    hit->normal = (denom > 0) ? -normal : normal;  // Assurer que la normale est opposée au rayon
    hit->uv = {(intersection_point.x / half_size + 1) * 0.5, (intersection_point.y / half_size + 1) * 0.5};

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le quad (rectangle).
AABB Quad::compute_aabb() {
    // Supposons que le quad soit sur le plan Z=0 et s'étende de half_size dans les directions X et Y.
    double epsilon = 1e-5; // Petite valeur pour assurer un volume
    double3 min_point = {-half_size, -half_size, -epsilon};
    double3 max_point = { half_size,  half_size,  epsilon};

    return AABB{min_point, max_point};
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un cylindre.
// Référez-vous au PDF pour la paramétrisation des coordonnées UV.
bool Cylinder::local_intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    // Calculer les coefficients pour l'équation quadratique dans le plan XZ
    double a = ray.direction.x * ray.direction.x + ray.direction.z * ray.direction.z;
    double b = 2.0 * (ray.origin.x * ray.direction.x + ray.origin.z * ray.direction.z);
    double c = ray.origin.x * ray.origin.x + ray.origin.z * ray.origin.z - radius * radius;

    // Calculer le discriminant
    double discriminant = b * b - 4 * a * c;

    // Si le discriminant est négatif, il n'y a pas d'intersection
    if (discriminant < 0) {
        return false;
    }

    // Calculer les deux valeurs possibles de t
    double sqrt_discriminant = sqrt(discriminant);
    double t1 = (-b - sqrt_discriminant) / (2.0 * a);
    double t2 = (-b + sqrt_discriminant) / (2.0 * a);

    // Vérifier la validité de t1 et t2 dans les limites
    double y1 = ray.origin.y + t1 * ray.direction.y;
    double y2 = ray.origin.y + t2 * ray.direction.y;

    bool valid_t1 = (t1 >= t_min && t1 <= t_max && y1 >= -half_height && y1 <= half_height);
    bool valid_t2 = (t2 >= t_min && t2 <= t_max && y2 >= -half_height && y2 <= half_height);

    // Choisir l'intersection valide la plus proche
    double t;
    if (valid_t1 && valid_t2) {
        t = (t1 < t2) ? t1 : t2;
    } else if (valid_t1) {
        t = t1;
    } else if (valid_t2) {
        t = t2;
    } else {
        return false;  // Pas d'intersection valide dans les limites de hauteur
    }

    // Calculer le point d'intersection et la normale
    double3 intersection_point = ray.origin + t * ray.direction;
    double3 normal = normalize(double3(intersection_point.x, 0, intersection_point.z));  // Ignorer Y pour la normale

    // Mettre à jour les informations d'intersection
    hit->depth = t;
    hit->position = intersection_point;
    hit->normal = (dot(ray.direction, normal) > 0) ? -normal : normal;  // Assurer que la normale est opposée au rayon
    hit->uv = {
        0.5 + atan2(normal.z, normal.x) / (2 * PI),
        (intersection_point.y + half_height) / (2 * half_height)
    };

    return true;
}


// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le cylindre.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire (comme ici).
AABB Cylinder::compute_aabb() {
    // Approximer le AABB du cylindre en espace local
    double3 min_point = {-radius, -radius, -half_height};
    double3 max_point = { radius,  radius, half_height};
    
    // Transformer les points min et max en espace global si nécessaire.
    return AABB{min_point, max_point};
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un mesh.
// Référez-vous au PDF pour la paramétrisation pour les coordonnées UV.
// Pour plus de d'informations sur la géométrie, référez-vous à la classe object.h.
bool Mesh::local_intersect(Ray const ray, double t_min, double t_max, Intersection* hit) {
    bool has_intersection = false;
    double closest_t = t_max;

    for (const auto& triangle : triangles) {
        Intersection temp_hit;
        if (intersect_triangle(ray, t_min, closest_t, triangle, &temp_hit)) {
            has_intersection = true;
            closest_t = temp_hit.depth;  // Mettre à jour la profondeur de l'intersection la plus proche

            // Interpoler la normale en utilisant les coordonnées barycentriques
            const Vertex& v0 = triangle.v[0];
            const Vertex& v1 = triangle.v[1];
            const Vertex& v2 = triangle.v[2];

            double w0 = temp_hit.uv.x;  // Coordonnée barycentrique w0
            double w1 = temp_hit.uv.y;  // Coordonnée barycentrique w1
            double w2 = 1.0 - w0 - w1;  // Coordonnée barycentrique w2

            double3 interpolated_normal = normalize(
                w0 * normals[v0.ni] +
                w1 * normals[v1.ni] +
                w2 * normals[v2.ni]
            );

            // Interpoler les coordonnées de texture (si applicable)
            double2 interpolated_uv = {
                w0 * tex_coords[v0.ti].x + w1 * tex_coords[v1.ti].x + w2 * tex_coords[v2.ti].x,
                w0 * tex_coords[v0.ti].y + w1 * tex_coords[v1.ti].y + w2 * tex_coords[v2.ti].y
            };

            // Mettre à jour les informations d'intersection
            hit->depth = temp_hit.depth;
            hit->position = temp_hit.position;
            hit->normal = interpolated_normal;
            hit->uv = interpolated_uv;
            hit->key_material = key_material;
        }
    }

    return has_intersection;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de trouver l'intersection avec un triangle.
// S'il y a intersection, remplissez hit avec l'information sur la normale et les coordonnées texture.
bool Mesh::intersect_triangle(Ray ray, double t_min, double t_max, Triangle const tri, Intersection *hit) {
    // Extraire chaque position de sommet des données du mesh
    double3 const &p0 = positions[tri[0].pi]; // Sommet A
    double3 const &p1 = positions[tri[1].pi]; // Sommet B
    double3 const &p2 = positions[tri[2].pi]; // Sommet C

    // Vecteurs des arêtes du triangle
    double3 edge1 = p1 - p0;
    double3 edge2 = p2 - p0;

    // Calculer le déterminant - aussi utilisé pour calculer le paramètre u
    double3 h = cross(ray.direction, edge2);
    double a = dot(edge1, h);

    // Si a est proche de zéro, le rayon est parallèle au triangle
    if (fabs(a) < EPSILON) {
        return false;
    }

    // Calculer la distance du sommet p0 à l'origine du rayon
    double f = 1.0 / a;
    double3 s = ray.origin - p0;

    // Calculer le paramètre u et vérifier la limite
    double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    // Calculer le paramètre v et vérifier la limite
    double3 q = cross(s, edge1);
    double v = f * dot(ray.direction, q);
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    // Calculer t, le rayon intersecte le triangle
    double t = f * dot(edge2, q);
    if (t < t_min || t > t_max) {
        return false;
    }

    // Remplir les données de l'intersection
    hit->depth = t;
    hit->position = ray.origin + t * ray.direction;
    hit->normal = normalize(cross(edge1, edge2));  // Normale du triangle
    hit->uv = { u, v };  // Coordonnées barycentriques en tant que coordonnées UV de texture

    return true;
}

// @@@@@@ VOTRE CODE ICI
// Occupez-vous de compléter cette fonction afin de calculer le AABB pour le Mesh.
// Il faut que le AABB englobe minimalement notre objet à moins que l'énoncé prononce le contraire.
AABB Mesh::compute_aabb() {
    // Initialiser les points min et max avec des valeurs extrêmes
    double3 min_point = {DBL_MAX, DBL_MAX, DBL_MAX};
    double3 max_point = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

    // Itérer à travers chaque position de sommet pour mettre à jour les limites min et max
    for (const auto& pos : positions) {
        min_point.x = std::min(min_point.x, pos.x);
        min_point.y = std::min(min_point.y, pos.y);
        min_point.z = std::min(min_point.z, pos.z);

        max_point.x = std::max(max_point.x, pos.x);
        max_point.y = std::max(max_point.y, pos.y);
        max_point.z = std::max(max_point.z, pos.z);
    }

    // Retourner le AABB en utilisant l'initialisation par accolades
    return AABB{min_point, max_point};
}
