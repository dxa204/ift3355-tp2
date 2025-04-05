# Ray Tracer – IFT 3355 TP2

## Project Description

The goal is to implement a basic ray tracer capable of rendering simple 3D scenes composed of geometric primitives like spheres and planes. The ray tracer simulates the interaction of light with objects, including effects such as reflections and shadows, to produce realistic images.

## Features

- **Ray-Sphere and Ray-Plane Intersections:** Calculate intersections between rays and scene objects.
- **Phong Lighting Model:** Implement ambient, diffuse, and specular lighting.
- **Shadow Casting:** Determine shadows from point light sources.
- **Recursive Reflections:** Simulate mirror-like reflections using recursion.
- **Scene Parsing:** Read scene descriptions from a text file.
- **Image Output:** Generate the final rendered image in PPM or PNG format.

## Project Structure

. ├── main.py # Main entry point for the ray tracer ├── raytracer.py # Core ray tracing engine
├── vector.py # Vector and geometric utility functions ├── image.py # Image creation and saving
functions ├── scene.txt # Sample scene description file └── README.md # Project documentation


## How to Run

1. **Install Python 3.x** if not already installed.
2. **Install required dependencies** using pip:

   ```bash
   pip install numpy Pillow
