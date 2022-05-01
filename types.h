#ifndef TYPES_H
#define TYPES_H

#include <vector>

struct vec3 {
    float x;
    float y;
    float z;
};

struct mtlcolor{
    //diffuse color (red,green,blue) component
    float Odr;
    float Odg;
    float Odb;
    //specular color (red,green,blue) component
    float Osr;
    float Osg;
    float Osb;
    /* weights for ambient, diffuse, and specular components of the illumination 
       at the interesection point based on the type of material */
    float ka;
    float kd;
    float ks;
    // weight to specify the width of the specular highlight
    float n;
    //material opacity
    float a;
    //index of relection
    float index_of_reflection;
};

struct Light{
    //light position (x,y,z)
    vec3 position;
    // directional light if w is 0; positional light if w is 1;
    float w;
    // light color (red,green,blue)
    float r;
    float g;
    float b;
};


struct Sphere {
    vec3 center;
    float radius;
    mtlcolor material;
};

struct Scene {
    std::vector<Sphere> spheres;
    std::vector<Light> lights;
    vec3 bkgcolor;
    vec3 eye;
    vec3 viewdir;
    vec3 updir;
    float width;
    float height;
    float vfov;
};
#endif