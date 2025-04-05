#include "raytracer.h"
#include <cstdlib>   // Pour rand() et srand()
#include <ctime>     // Pour time()
#include <algorithm> // Pour std::max

void Raytracer::render(const Scene& scene, Frame* output)
{       
    srand(static_cast<unsigned>(time(0))); // Initialiser avec l'heure actuelle

    // Créer le tampon z (z-buffer)
    double *z_buffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        z_buffer[i] = scene.camera.z_far; // Remplacez l'ancienne valeur DBL_MAX par scene.camera.z_far
    }

    //---------------------------------------------------------------------------------------------------------------
    // Nous vous fournissons ci-bas du code pour une caméra orthographique très simple. Cette caméra peut être utilisée pour tester l’intersection avec la sphère.
    // Vous devez utiliser la scène de test portho.ray pour utiliser cette caméra. 
    // Notez que votre code de caméra ne doit pas être basé sur ce code de caméra. Ce code n’est là que pour prendre en compte le développement initial du test d’intersection.
    // Pour utiliser cette caméra, vous devez supprimer les commentaires qui rendent inactive cette partie du code, et mettre en commentaires la boucle d’image originale.

    /*CameraOrthographic camOrth;
    double3 uVec{ 0,1,0 };
    double3 vVec{ 0,0,1 };
    double y_shift = 2.0 / scene.resolution[1];
    double x_shift = 2.0 / scene.resolution[0];

    for (int y = 0; y < scene.resolution[1]; y++) {
        if (y % 40) {
            std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
        }

        for (int x = 0; x < scene.resolution[0]; x++) {
            double3 color{ 0,0,0 };
            Intersection hit;
            double3 rayOrigin = camOrth.minPosition + uVec * x_shift * x + vVec * y_shift * y;
            double3 rayDirection{ 1,0,0 };
            Ray ray = Ray(rayOrigin, rayDirection);
            double itHits = 0;
            double z_depth = scene.camera.z_far;
            if (scene.container->intersect(ray, EPSILON, z_depth, &hit)) {
                Material& material = ResourceManager::Instance()->materials[hit.key_material];
                color = material.color_albedo;
                itHits = 1.0f;
            }
            output->set_color_pixel(x, y, color);
            output->set_depth_pixel(x, y, itHits);
        }
    }*/

    //---------------------------------------------------------------------------------------------------------------

    // ATTENTION
    // La code commente ci-dessous est la code a utiliser si on ne veut pas incorporer la Profondeur de Champ (partie 2.5)
    // Si vous voulez voir comment la code fonctionne sans la Profondeur de Champ, vous pouvez decommenter la code ci-dessous
    // et commenter la reste de la fonction render en commencant de la ligne 117


    // @@@@@@ VOTRE CODE ICI
    // Calculez les paramètres de la caméra pour les rayons.
    /*
    double3 camW = normalize(scene.camera.position - scene.camera.center);  // direction z
    double3 camU = normalize(cross(camW, scene.camera.up));                 // direction x
    double3 camV = cross(camW, camU);                                       // direction y

    double alpha = scene.camera.z_near / camW[2];
    double3 centerPOV = scene.camera.position - alpha * camW;
    double heightPOV = 2 * alpha * tan(0.5 * deg2rad(scene.camera.fovy));
    double widthPOV = heightPOV * scene.camera.aspect;
    double pixelWidth = widthPOV / scene.resolution[0];
    double pixelHeight = heightPOV / scene.resolution[1];
    double3 bottomLeftCornerPOV = centerPOV - 0.5 * heightPOV * camV - 0.5 * widthPOV * camU;

    double jitter_x = (static_cast<double>(rand()) / RAND_MAX - 0.5) * pixelWidth;
    double jitter_y = (static_cast<double>(rand()) / RAND_MAX - 0.5) * pixelHeight;

    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
        if (y % 40){
            std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
        }

        for(int x = 0; x < scene.resolution[0]; x++) {
            int avg_z_depth = 0;
            double3 avg_ray_color{0,0,0};
            
            for(int iray = 0; iray < scene.samples_per_pixel; iray++) {
                // Utiliser rand() pour générer du jitter pour l'anti-aliasing
                double jitter_x = (static_cast<double>(rand()) / RAND_MAX - 0.5) * pixelWidth;
                double jitter_y = (static_cast<double>(rand()) / RAND_MAX - 0.5) * pixelHeight;
                
                double3 pixel = bottomLeftCornerPOV + (x + 0.5 + jitter_x) * pixelWidth * camU +
                                (y + 0.5 + jitter_y) * pixelHeight * camV;

                Ray ray(scene.camera.position, normalize(pixel - scene.camera.position));
                double3 ray_color{0, 0, 0};
                double z_depth = scene.camera.z_far;
                trace(scene, ray, 0, &ray_color, &z_depth);
                avg_ray_color += ray_color;
                avg_z_depth += z_depth;
            }
            
            avg_z_depth = avg_z_depth / scene.samples_per_pixel;
            avg_ray_color = avg_ray_color / scene.samples_per_pixel;

            // Test de profondeur
            if(avg_z_depth >= scene.camera.z_near && avg_z_depth <= scene.camera.z_far && 
                avg_z_depth < z_buffer[x + y*scene.resolution[0]]) {
                z_buffer[x + y*scene.resolution[0]] = avg_z_depth;
                output->set_color_pixel(x, y, avg_ray_color);
                output->set_depth_pixel(x, y, (avg_z_depth - scene.camera.z_near) / 
                                        (scene.camera.z_far - scene.camera.z_near));
            }
        }
    }*/

    // Paramètres de profondeur de champ
    double focal_distance = 10.0;   // Distance du plan focal
    double aperture_size = 0.1;     // Taille de l'ouverture, valeurs plus grandes augmentent le flou pour les zones non focalisées
    int num_samples = 16;           // Nombre d'échantillons par pixel pour l'effet de flou
    
    double3 direction = normalize(scene.camera.center - scene.camera.position);

    double view_height = 2.0 * tan(scene.camera.fovy * 0.5 * PI / 180.0);
    double view_width = view_height * scene.camera.aspect;
    double pixel_size_x = view_width / scene.resolution[0];
    double pixel_size_y = view_height / scene.resolution[1];

    for (int y = 0; y < scene.resolution[1]; y++) {
        if (y % 40 == 0) {
            std::cout << "\rScanlines completed: " << y << "/" << scene.resolution[1] << '\r';
        }

        for (int x = 0; x < scene.resolution[0]; x++) {
            double3 color{0, 0, 0};
            double z_depth = scene.camera.z_far;

            // Accumuler la couleur de plusieurs échantillons pour simuler la profondeur de champ
            for (int s = 0; s < num_samples; s++) {
                double random_angle = ((double)rand() / RAND_MAX) * 2.0 * PI;
                double random_radius = ((double)rand() / RAND_MAX) * aperture_size;
                double aperture_x = random_radius * cos(random_angle);
                double aperture_y = random_radius * sin(random_angle);

                double3 ray_origin = scene.camera.position + double3{aperture_x, aperture_y, 0};

                double3 target = scene.camera.position + direction * focal_distance;
                double3 focal_point = target + double3(
                    (x - scene.resolution[0] / 2.0) * pixel_size_x,
                    (y - scene.resolution[1] / 2.0) * pixel_size_y,
                    0
                );
                double3 ray_direction = normalize(focal_point - ray_origin);

                Ray ray(ray_origin, ray_direction);
                Intersection hit;

                if (scene.container->intersect(ray, EPSILON, z_depth, &hit)) {
                    color += shade(scene, hit);
                }
            }

            color /= num_samples;
            output->set_color_pixel(x, y, color);
            z_buffer[y * scene.resolution[0] + x] = z_depth;
        }
    }
    delete[] z_buffer;
}

// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
//  - Détermine si le rayon intersecte la géométrie.
//  - Calculer la contribution associée à la réflexion.
//  - Calculer la contribution associée à la réfraction.
//  - Mettre à jour la couleur avec le shading + Ajouter réflexion selon material.reflection +
//    Ajouter réfraction selon material.refraction pour la couleur de sortie.
//  - Mettre à jour la nouvelle profondeur.

void Raytracer::trace(const Scene& scene, Ray ray, int depth, linalg::aliases::double3* out_color, double* out_z_depth)
{
    // Vérifier la profondeur de récursion maximale, arrêter si elle est atteinte.
    if (depth <= 0) {
        *out_color = double3(0, 0, 0);  // Aucun effet de réflexion ou de réfraction supplémentaire
        return;
    }

    Intersection hit;
    // Vérifier si le rayon intersecte un objet dans la scène
    if (scene.container->intersect(ray, EPSILON, DBL_MAX, &hit)) {
        // Calcule la couleur locale en fonction du matériau et de l'éclairage
        double3 local_color = shade(scene, hit);

        // Met à jour la profondeur de l'intersection si elle est demandée
        if (out_z_depth) {
            *out_z_depth = hit.depth;
        }

        // Récupère le matériau de l'objet pour appliquer réflexion et réfraction
        Material& material = ResourceManager::Instance()->materials[hit.key_material];

        // Calcule la couleur de réflexion
        double3 reflection_color = double3(0, 0, 0);
        if (material.k_reflection > 0) {
            // Calcul de la direction de réflexion (loi de réflexion)
            double3 reflect_dir = normalize(ray.direction - 2 * dot(ray.direction, hit.normal) * hit.normal);
            // Créer un rayon réfléchi légèrement décalé pour éviter l'auto-intersection
            Ray reflection_ray(hit.position + hit.normal * EPSILON, reflect_dir);
            // Tracer récursivement le rayon réfléchi en diminuant la profondeur
            trace(scene, reflection_ray, depth - 1, &reflection_color, nullptr);
        }

        // Calcule la couleur de réfraction
        double3 refraction_color = double3(0, 0, 0);
        if (material.k_refraction > 0) {
            double eta = 1.0 / material.refractive_index;  // Indice d'air à celui du matériau
            double cos_i = dot(-ray.direction, hit.normal);
            double sin2_t = eta * eta * (1.0 - cos_i * cos_i);

            // Vérifier s'il y a réflexion totale interne
            if (sin2_t <= 1.0) {
                // Calcul de la direction réfractée si pas de réflexion totale
                double cos_t = sqrt(1.0 - sin2_t);
                double3 refract_dir = normalize(eta * ray.direction + (eta * cos_i - cos_t) * hit.normal);
                // Créer un rayon réfracté, légèrement décalé pour éviter l'auto-intersection
                Ray refraction_ray(hit.position - hit.normal * EPSILON, refract_dir);
                // Tracer récursivement le rayon réfracté en diminuant la profondeur
                trace(scene, refraction_ray, depth - 1, &refraction_color, nullptr);
            }
        }

        // Combine la couleur locale, la réflexion et la réfraction
        *out_color = (1 - material.k_reflection - material.k_refraction) * local_color +
                     material.k_reflection * reflection_color +
                     material.k_refraction * refraction_color;
    } else {
        // Si aucun objet n'est intersecté, utiliser la couleur d'éclairage ambiant de la scène
        *out_color = scene.ambient_light;
    }
}


// @@@@@@ VOTRE CODE ICI
// Veuillez remplir les objectifs suivants:
//  * Calculer la contribution des lumières dans la scène.
//      - Itérer sur toutes les lumières.
//      - Inclure la contribution spéculaire selon le modèle de Blinn en incluant la composante métallique.
//      - Inclure la contribution diffuse. (Faites attention au produit scalaire. >= 0)
//      - Inclure la contribution ambiante
//  * Calculer si le point est dans l'ombre
//      - Itérer sur tous les objets et détecter si le rayon entre l'intersection et la lumière est occulté.
//      - Ne pas considérer les points plus loins que la lumière.
//      - Par la suite, intégrer la pénombre dans votre calcul
//  * Déterminer la couleur du point d'intersection.
//      - Si une texture est présente, prendre la couleur aux coordonnées uv
//      - Sinon, utiliser la couleur associée au matériel.

double3 Raytracer::shade(const Scene& scene, Intersection hit)
{
    // Récupère le matériau de l'objet à l'endroit de l'intersection
    Material& material = ResourceManager::Instance()->materials[hit.key_material];
    double3 color;

    // Vérifie si une texture est associée au matériau en utilisant les dimensions de la texture
    if (material.texture_albedo.width() > 0 && material.texture_albedo.height() > 0) {
        // Texture présente; obtenir la couleur en fonction des coordonnées UV
        double u = hit.uv.x;
        double v = hit.uv.y;

        // Calcule la position dans la texture en fonction des coordonnées UV
        int tex_x = static_cast<int>(u * material.texture_albedo.width());
        int tex_y = static_cast<int>(v * material.texture_albedo.height());

        // Récupère la couleur du pixel dans la texture
        rgb_t tex_color;
        material.texture_albedo.get_pixel(tex_x, tex_y, tex_color);

        // Normalise la couleur de [0,255] à [0,1] pour l'albédo
        color = double3(tex_color.red / 255.0, tex_color.green / 255.0, tex_color.blue / 255.0);
    } else {
        // Utiliser l'albédo du matériau si aucune texture n'est présente
        color = material.color_albedo;
    }

    // Contribution ambiante : ajout de la lumière ambiante multipliée par le coefficient ambiant du matériau
    double3 final_color = material.k_ambient * scene.ambient_light * color;

    // Itère sur toutes les lumières de la scène pour calculer la contribution diffuse et spéculaire
    for (const auto& light : scene.lights) {
        // Calcule la direction de la lumière et la distance de l'intersection à la source lumineuse
        double3 light_dir = normalize(light.position - hit.position);
        double light_dist = length(light.position - hit.position);

        // Crée un rayon d'ombre depuis le point d'intersection vers la lumière pour vérifier l'occlusion
        Ray shadow_ray(hit.position + hit.normal * EPSILON, light_dir);
        Intersection shadow_hit;
        bool in_shadow = scene.container->intersect(shadow_ray, EPSILON, light_dist, &shadow_hit);

        // Si le point n'est pas dans l'ombre de cette lumière, on ajoute ses contributions
        if (!in_shadow) {
            // Composante diffuse (réflexion Lambertienne) en fonction de l'angle entre la normale et la lumière
            double diffuse_intensity = std::max(0.0, dot(hit.normal, light_dir));
            final_color += material.k_diffuse * color * diffuse_intensity * light.emission;

            // Composante spéculaire (modèle de Blinn-Phong) pour simuler les reflets brillants
            double3 view_dir = normalize(scene.camera.position - hit.position);
            double3 halfway_dir = normalize(light_dir + view_dir);
            double spec_intensity = pow(std::max(0.0, dot(hit.normal, halfway_dir)), material.shininess);
            final_color += material.k_specular * color * spec_intensity * light.emission * (1 - material.metallic);
        }
    }

    // Retourne la couleur finale résultant de la somme des contributions ambiantes, diffuses et spéculaires
    return final_color;
}

