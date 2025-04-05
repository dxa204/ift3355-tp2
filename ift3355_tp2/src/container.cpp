#include "container.h"
#include <stack>

// @@@@@@ VOTRE CODE ICI
// - Parcourir l'arbre DEPTH FIRST SEARCH selon les conditions suivantes:
//       - S'il s'agit d'une feuille, faites l'intersection avec la géométrie.
//       - Sinon, il s'agit d'un noeud altérieur.
//           - Faites l'intersection du rayon avec le AABB gauche et droite. 
//           - S'il y a intersection, ajouter le noeud à ceux à visiter. 
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool BVH::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool has_intersection = false;
    double closest_so_far = t_max; // Suivre l'intersection la plus proche dans l'intervalle [t_min, t_max]

    // Recherche en profondeur basée sur une pile
    std::stack<BVHNode*> stack;
    stack.push(root); // Commencer le parcours à partir de la racine du BVH

    while (!stack.empty()) {
        BVHNode* node = stack.top();
        stack.pop();

        // Vérifier l'intersection avec l'AABB du nœud actuel
        if (!node->aabb.intersect(ray, t_min, closest_so_far)) {
            continue; // Passer ce nœud s'il n'y a pas d'intersection avec son AABB
        }

        if (node->is_leaf()) {
            // Nœud feuille : vérifier l'intersection avec la géométrie réelle
            Intersection temp_hit;
            if (node->object->intersect(ray, t_min, closest_so_far, &temp_hit)) {
                has_intersection = true;
                closest_so_far = temp_hit.depth; // Mettre à jour la profondeur de l'intersection la plus proche
                *hit = temp_hit; // Mettre à jour l'enregistrement de l'intersection la plus proche
            }
        } else {
            // Nœud interne : parcourir les nœuds enfants gauche et droit
            if (node->left) stack.push(node->left);
            if (node->right) stack.push(node->right);
        }
    }

    return has_intersection;
}

// @@@@@@ VOTRE CODE ICI
// - Parcourir tous les objets
//       - Détecter l'intersection avec l'AABB
//           - Si intersection, détecter l'intersection avec la géométrie.
//           - Si intersection, mettre à jour les paramètres.
// - Retourner l'intersection avec la profondeur maximale la plus PETITE.
bool Naive::intersect(Ray ray, double t_min, double t_max, Intersection* hit) {
    bool has_intersection = false;
    double closest_so_far = t_max; // Suivre l'intersection la plus proche dans l'intervalle [t_min, t_max]

    for (const auto& object : objects) {
        // Vérifier l'intersection avec l'AABB de l'objet
        if (object->compute_aabb().intersect(ray, t_min, closest_so_far)) {
            // Si le rayon intersecte l'AABB, vérifier l'intersection avec la géométrie réelle
            Intersection temp_hit;
            if (object->intersect(ray, t_min, closest_so_far, &temp_hit)) {
                has_intersection = true;
                closest_so_far = temp_hit.depth; // Mettre à jour la profondeur de l'intersection la plus proche
                *hit = temp_hit; // Mettre à jour l'enregistrement de l'intersection la plus proche
            }
        }
    }

    return has_intersection;
}
