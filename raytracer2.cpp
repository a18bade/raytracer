#include <iostream>
#include <cmath>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include "types.h"
#define PI 3.14159265

Scene parse_file(std::string filename)
{
    Scene scene;
    FILE *fp = fopen(filename.c_str(), "r");

    if (fp == nullptr)
    {
        std::cout << "Cannot open given file";
        exit(-1);
    }
    size_t max_length = 1024;
    char line[max_length];
    mtlcolor material = {
        .Odr = 0.0,
        .Odg = 0.0,
        .Odb = 0.0,
        .Osr = 0.0,
        .Osg = 0.0,
        .Osb = 0.0,
        .ka = 0.0,
        .kd = 0.0,
        .ks = 0.0,
        .n = 0.0,
        .a = 0.0,
        .index_of_reflection = 0.0,
    };

    while (fgets(line, max_length, fp))
    {
        size_t keywordlength = 32;
        char keyword[keywordlength];
        int32_t num_fields_read = sscanf(line, "%s ", keyword);

        if (num_fields_read < 1)
        {
            continue;
        }
        else if (strcmp(keyword, "sphere") == 0)
        {
            vec3 c;
            float r;

            sscanf(line, "sphere %f %f %f %f", &c.x, &c.y, &c.z, &r);
            Sphere sphere = {
                .center = c,
                .radius = r,
                .material = material,
            };
            scene.spheres.push_back(sphere);
            // printf("sphere (position: %f %f %f )\t (radius: %f)\n",
            // sphere.center.x,sphere.center.y,sphere.center.z,
            // sphere.radius);
        }
        else if (strcmp(keyword, "mtlcolor") == 0)
        {
            sscanf(line, "mtlcolor %f %f %f %f %f %f %f %f %f %f %f %f",
                   &material.Odr, &material.Odg, &material.Odb,
                   &material.Osr, &material.Osg, &material.Osb,
                   &material.ka, &material.kd, &material.ks,
                   &material.n, &material.a, &material.index_of_reflection);
            // printf("material - diffuse color :(%f,%f,%f)\t specular color :(%f,%f,%f)\t weights: (%f,%f,%f)\t, n: %f \t a: %f index_ref: %f\n",
            //        material.Odr, material.Odg, material.Odb,
            //        material.Osr, material.Osg, material.Osb,
            //        material.ka, material.kd, material.ks,
            //        material.n, material.a, material.index_of_reflection);
        }
        else if (strcmp(keyword, "light") == 0)
        {
            Light light;
            sscanf(line, "light %f %f %f %f %f %f %f",
                   &light.position.x, &light.position.y, &light.position.z,
                   &light.w, &light.r, &light.g, &light.b);
            scene.lights.push_back(light);

            // //printf("light (position: %f, %f, %f)\t (w: %f)\t (color: %f, %f, %f)\t\n",
            // light.position.x,light.position.y,light.position.z,
            // light.w, light.r, light.g, light.b);
        }
        else if (strcmp(keyword, "eye") == 0)
        {
            vec3 c;
            sscanf(line, "eye %f %f %f", &c.x, &c.y, &c.z);
            scene.eye = c;
            // printf("eye %f %f %f\n",scene.eye.x,scene.eye.y,scene.eye.z);
        }
        else if (strcmp(keyword, "viewdir") == 0)
        {
            vec3 c;
            sscanf(line, "viewdir %f %f %f", &c.x, &c.y, &c.z);
            scene.viewdir = c;
            //printf("viewdir %f %f %f\n",scene.viewdir.x,scene.viewdir.y,scene.viewdir.z);
        }
        else if (strcmp(keyword, "updir") == 0)
        {
            vec3 c;
            sscanf(line, "updir %f %f %f", &c.x, &c.y, &c.z);
            scene.updir = c;
            //printf("updir %f %f %f\n",scene.updir.x,scene.updir.y,scene.updir.z);
        }
        else if (strcmp(keyword, "vfov") == 0)
        {
            float v;
            sscanf(line, "vfov %f", &v);
            scene.vfov = v;
            //printf("vfov %f\n",scene.vfov);
        }
        else if (strcmp(keyword, "imsize") == 0)
        {
            sscanf(line, "imsize %f %f", &scene.width, &scene.height);
            //printf("imsize %f %f\n",scene.width,scene.height);
        }
        else if (strcmp(keyword, "bkgcolor") == 0)
        {
            vec3 c;
            sscanf(line, "bkgcolor %f %f %f", &c.x, &c.y, &c.z);
            scene.bkgcolor = c;
            //printf("bkgcolor %f %f %f\n",scene.bkgcolor.x,scene.bkgcolor.y,scene.bkgcolor.z);
        }
    }
    return scene;
}

//basic vector operations:

vec3 add(vec3 a, vec3 b)
{
    vec3 sum = {
        .x = a.x + b.x,
        .y = a.y + b.y,
        .z = a.z + b.z,
    };
    return sum;
}

vec3 sub(vec3 a, vec3 b)
{
    vec3 difference = {
        .x = a.x - b.x,
        .y = a.y - b.y,
        .z = a.z - b.z,
    };
    return difference;
}

vec3 scale(float scalar, vec3 a)
{
    vec3 newVec = {
        .x = scalar * a.x,
        .y = scalar * a.y,
        .z = scalar * a.z,
    };
    return newVec;
}

float dot(vec3 a, vec3 b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

vec3 cross(vec3 a, vec3 b)
{
    vec3 newVec = {
        .x = (a.y * b.z) - (a.z * b.y),
        .y = (a.z * b.x) - (a.x * b.z),
        .z = (a.x * b.y) - (a.y * b.x),
    };
    return newVec;
}

float length(vec3 a)
{
    return sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2));
}

vec3 normalize(vec3 a)
{
    float len = length(a);
    vec3 newVec = {
        .x = a.x / len,
        .y = a.y / len,
        .z = a.z / len,
    };
    return newVec;
}

float DegreesToRadians(float degrees)
{
    return degrees * PI / 180.0;
}

float RadiansToDegrees(float radians)
{
    return radians * 180 / PI;
}

std::string toString(vec3 a)
{
    return "[ " + std::to_string(a.x) + "," + std::to_string(a.y) + "," + std::to_string(a.z) + " ]";
}

void debug()
{

    printf("Testing\n");

    vec3 vector1 = {
        .x = 1,
        .y = 2,
        .z = 2,
    };

    vec3 vector2 = {
        .x = 3,
        .y = 4,
        .z = 5,
    };

    vec3 vector3 = {
        .x = 0,
        .y = 3,
        .z = 4,
    };
    printf("\n");

    printf("30 degrees to radians is: %f\n", DegreesToRadians(30));
    printf("Correct Answer: 0.523599\n");

    printf("30 radians to degrees is: %f\n", RadiansToDegrees(30));
    printf("Correct Answer: 1718.87\n");

    printf("Length of (1,2,2) is %f\n", length(vector1));
    printf("Correct Answer: 3\n");

    printf("Sum of (1,2,2) and (3,4,5) is %s\n", toString(add(vector1, vector2)).c_str());
    printf("Correct Answer: (4,6,7)\n");

    printf("Diff of (1,2,2) and (3,4,5) is %s\n", toString(sub(vector1, vector2)).c_str());
    printf("Correct Answer: (-2,-2,-3)\n");

    printf("Scale of (1,2,2) by 6 is %s\n", toString(scale(6, vector1)).c_str());
    printf("Correct Answer: (6,12,12)\n");

    printf("Dot of (1,2,2) and (0,3,4) is %f\n", dot(vector1, vector3));
    printf("Correct Answer: 14\n");

    printf("Cross of (1,2,2) and (0,3,4) is %s\n", toString(cross(vector1, vector3)).c_str());
    printf("Correct Answer: (2,-4,3)\n");
}

float Shadow_Trace_Ray(vec3 intersectionPoint, vec3 rayTowardsLight, Scene scene, float distanceTowardsLight, float typeofLight, int sphereindex)
{
    float S = 1; //assume the object isnt in shadow
    float x0 = intersectionPoint.x;
    float y0 = intersectionPoint.y;
    float z0 = intersectionPoint.z;

    float xD = rayTowardsLight.x;
    float yD = rayTowardsLight.y;
    float zD = rayTowardsLight.z;

    float t = -1;
    float tminus;
    float tplus;

    for (int s = 0; s < scene.spheres.size(); s++)
    {
        if (s != sphereindex)
        {

            float xC = scene.spheres[s].center.x;
            float yC = scene.spheres[s].center.y;
            float zC = scene.spheres[s].center.z;
            float r = scene.spheres[s].radius;

            float a = (xD * xD) + (yD * yD) + (zD * zD);
            float b = 2 * (xD * (x0 - xC) + yD * (y0 - yC) + zD * (z0 - zC));
            float c = pow((x0 - xC), 2) + pow((y0 - yC), 2) + pow((z0 - zC), 2) - pow(r, 2);

            float discriminant = pow(b, 2) - (4 * a * c);

            if (discriminant > 0.0)
            {
                tminus = (-b - sqrt(discriminant)) / (2 * a);
                tplus = (-b + sqrt(discriminant)) / (2 * a);
                if (tminus > 0)
                {
                    t = tminus;
                }
                if (tplus > 0)
                {

                    t = tplus;
                }
            }
        }
    }
    if (typeofLight == 0)
    {
        if (t < 0)
        {
            S = 1.0;
        }
        else
        {
            S = 0.0;
        }
    }
    if (typeofLight == 1)
    {
        if (t < distanceTowardsLight && t > 0)
        {
            S = 0.0;
        }
        else
        {
            S = 1.0;
        }
    }
    return S;
}

vec3 Shade_Ray(Sphere sphere, float dist, vec3 ray, Scene scene)
{

    vec3 intersectionPoint = add(scene.eye, scale(dist, ray));

    vec3 N = normalize(sub(intersectionPoint, sphere.center));
    vec3 L;
    vec3 H;
    vec3 V = normalize(sub(scene.eye, intersectionPoint));

    vec3 color = {
        .x = 0.0,
        .y = 0.0,
        .z = 0.0,
    };

    float sumColorX;
    float sumColorY;
    float sumColorZ;
    float sflag;
    color.x += sphere.material.ka * sphere.material.Odr; //correct
    color.y += sphere.material.ka * sphere.material.Odg; //correct
    color.z += sphere.material.ka * sphere.material.Odb; //correct

    for (int i = 0; i < scene.lights.size(); i++)
    {
        if (scene.lights[i].w == 0)
        { //directional light
            L = normalize(scale(-1, scene.lights[i].position));
        }
        else
        { //positional light
            L = normalize(sub(scene.lights[i].position, intersectionPoint));
        }
        H = normalize(add(L, V));
        float distanceTowardsLight = length(sub(scene.lights[i].position, intersectionPoint));
        sflag = Shadow_Trace_Ray(intersectionPoint, L, scene, distanceTowardsLight, scene.lights[i].w, i);

        ///ignore shadow
        sflag = 1;

        //diffuse component
        sumColorX = sflag * (scene.lights[i].r) * (sphere.material.kd * sphere.material.Odr * fmax(0, dot(N, L)));
        sumColorY = sflag * (scene.lights[i].g) * (sphere.material.kd * sphere.material.Odg * fmax(0, dot(N, L)));
        sumColorZ = sflag * (scene.lights[i].b) * (sphere.material.kd * sphere.material.Odb * fmax(0, dot(N, L)));

        //specular component
        sumColorX += sflag * (scene.lights[i].r) * sphere.material.ks * sphere.material.Osr * pow(fmax(0, dot(N, H)), sphere.material.n);
        sumColorY += sflag * (scene.lights[i].g) * sphere.material.ks * sphere.material.Osg * pow(fmax(0, dot(N, H)), sphere.material.n);
        sumColorZ += sflag * (scene.lights[i].b) * sphere.material.ks * sphere.material.Osb * pow(fmax(0, dot(N, H)), sphere.material.n);

        color.x += sumColorX;
        color.y += sumColorY;
        color.z += sumColorZ;

        sumColorX = 0;
        sumColorY = 0;
        sumColorZ = 0;
        // printf("specular x is : %f\n", color.x);
        // printf("specular y is : %f\n", color.y);
        // printf("specular z is : %f\n", color.z);
        // printf("ks %f\n", sphere.material.ks);
        // printf("osr %f\n", sphere.material.Osr);
        // printf("n*h is %f\n",fmax(0,dot(N,H)));
        // printf("(n*h)^n %f\n",pow(fmax(0,dot(N,H)),sphere.material.n));

        // printf("color.r: %f\n",color.x);
        // printf("color.g: %f\n",color.y);
        // printf("color.b: %f\n",color.z);
    }

    if (color.x > 1)
    {
        color.x = 1;
    }
    if (color.y > 1)
    {
        color.y = 1;
    }
    if (color.z > 1)
    {
        color.z = 1;
    }
    return color;

    //print statements for testing
    // printf("dist: %f\n",dist);
    // // printf("ray dir %s\n", toString(ray).c_str());
    // // printf("eye %s\n", toString(scene.eye).c_str());
    // printf("intersection_point %s\n", toString(intersectionPoint).c_str());
    // printf("n: %s\n", toString(N).c_str());
    // printf("l: %s\n", toString(L).c_str());
    // printf("v: %s\n", toString(V).c_str());
    // printf("h: %s\n", toString(H).c_str());
    // printf("n dot l: %f\n", dot(N,L));
    // printf("n dot h: %f\n", dot(N,H));
    // printf("color.r: %f\n",color.x);
    // printf("color.g: %f\n",color.y);
    // printf("color.b: %f\n",color.z);
    //exit(-1);
}

vec3 calculateReflectionRay(vec3 ray_origin, vec3 intersectionPoint, vec3 surface_normal)
{
    surface_normal = normalize(surface_normal);
    vec3 incoming_Ray = normalize(sub(intersectionPoint, ray_origin));
    vec3 I = scale(-1, incoming_Ray);
    float a = dot(surface_normal, I);
    vec3 reflectedRay = normalize(sub(scale(2.0 * a, surface_normal), I));
    return reflectedRay;
}

vec3 calculateTransparentRay(vec3 ray_origin,vec3 intersectionPoint,vec3 surface_normal,float indexR_incoming,float indexR_outgoing){
    vec3 T = {
        .x = 0,
        .y = 0,
        .z = 0,
    };
    surface_normal = normalize(surface_normal);
    vec3 incoming_Ray = normalize(sub(intersectionPoint, ray_origin));
    vec3 I = scale(-1, incoming_Ray);
    float ratio = indexR_incoming/indexR_outgoing;
    float a = dot(surface_normal, I);
    float underSquareRoot = (1-(pow(ratio,2)*(1-pow(a,2))));
    if (underSquareRoot < 0){
        //internal reflection: T doesnt exist
        return T;
    }
    vec3 reverseN = scale(-1,surface_normal);
    vec3 firstTerm = scale(sqrt(underSquareRoot),reverseN);
    vec3 secondTerm = scale(ratio,sub(scale(a,surface_normal),I));
    T = add(firstTerm,secondTerm);

    T = normalize(T);
    
    return T;

    

}

float calculateFresnelCoeffSpecular(Sphere sphere, vec3 normal, vec3 I)
{
    float eta = sphere.material.index_of_reflection;
    float Fo = pow((eta - 1) / (eta + 1), 2);
    float cos_angle = dot(I, normal);
    float Fr = Fo + ((1 - Fo) * pow((1 - cos_angle), 5));
    return Fr;
}

float calculateFresnelCoeffTransparent(float nt, float ni, vec3 I, vec3 normal){
    float Fo = pow(((nt-ni)/(nt+ni)),2); //correct
    float cos_angle = dot(I, normal);
    float Fr = Fo + ((1 - Fo) * pow((1 - cos_angle), 5)); //correct
    return Fr;

}
vec3 Recursive_Trace_Ray(vec3 ray, vec3 origin, Scene scene, int depth)
{
    if (depth == 3)
    {
        vec3 color = {
            .x = 0,
            .y = 0,
            .z = 0,
        };
        return color;
    }

    vec3 color = scene.bkgcolor;
    bool hasCollided = false;
    Sphere closestSphere;

    float x0 = origin.x;
    float y0 = origin.y;
    float z0 = origin.z;

    float xD = ray.x;
    float yD = ray.y;
    float zD = ray.z;

    float minimumt = 100000.0;
    float tminus;
    float tplus;

    for (int s = 0; s < scene.spheres.size(); s++)
    {

        float xC = scene.spheres[s].center.x;
        float yC = scene.spheres[s].center.y;
        float zC = scene.spheres[s].center.z;
        float r = scene.spheres[s].radius;

        float a = (xD * xD) + (yD * yD) + (zD * zD);
        float b = 2 * (xD * (x0 - xC) + yD * (y0 - yC) + zD * (z0 - zC));
        float c = pow((x0 - xC), 2) + pow((y0 - yC), 2) + pow((z0 - zC), 2) - pow(r, 2);

        float discriminant = pow(b, 2) - (4 * a * c);

        if (discriminant > 0.0)
        {
            tminus = (-b - sqrt(discriminant)) / (2 * a);
            tplus = (-b + sqrt(discriminant)) / (2 * a);

            if (tminus<minimumt & tminus> 0.001f)
            {
                hasCollided = true;
                minimumt = tminus;
                closestSphere = scene.spheres[s];
            }
            if (tplus<minimumt & tplus> 0.001f)
            {
                hasCollided = true;
                minimumt = tplus;
                closestSphere = scene.spheres[s];
            }
        }
    }
    if (hasCollided)
    {
        color = Shade_Ray(closestSphere, minimumt, ray, scene);
        vec3 intersectionPoint = add(ray, scale(minimumt, ray));
        vec3 N = normalize(sub(intersectionPoint, closestSphere.center));
        vec3 reflectedRay = calculateReflectionRay(ray, intersectionPoint, N);
        //get I
        vec3 incoming_Ray = normalize(sub(intersectionPoint, ray));
        vec3 I = scale(-1, incoming_Ray);
        //calculate fresnel affect
        float fresnelAffect = calculateFresnelCoeffSpecular(closestSphere, N, I);
        vec3 intermediateColor = Recursive_Trace_Ray(reflectedRay, intersectionPoint, scene, depth + 1);
        intermediateColor = scale(fresnelAffect, intermediateColor);

        if (color.x < 1)
        {
            if (color.y < 1)
            {
                if (color.z < 1)
                {
                    color = add(color, intermediateColor);
                }
            }
        }

        //transparent portion
        //assumption no overlapping spheres and the rays only go through spheres and air ,mentioned in class
        float fresnelAffect2 = calculateFresnelCoeffTransparent(1.0,closestSphere.material.index_of_reflection,I,N);
        vec3 T = calculateTransparentRay(ray,intersectionPoint,N,1.0,closestSphere.material.index_of_reflection);
        vec3 origin = intersectionPoint;
        vec3 intersectionPoint2 = add(intersectionPoint,T);
        vec3 surface_normal = scale(-1,normalize(sub(intersectionPoint2, closestSphere.center)));
        T = calculateTransparentRay(origin,intersectionPoint2,surface_normal,closestSphere.material.index_of_reflection,1.0);
        vec3 intermediateColor2 = Recursive_Trace_Ray(T, intersectionPoint2, scene, depth + 1);
        intermediateColor2 = scale(fresnelAffect2, intermediateColor2);
        if (color.x < 1)
        {
            if (color.y < 1)
            {
                if (color.z < 1)
                {
                    color = add(color, intermediateColor2);
                }
            }
        }
    }
    return color;
}

vec3 Trace_Ray(vec3 ray, Scene scene)
{

    vec3 color = scene.bkgcolor;
    bool hasCollided = false;
    Sphere closestSphere;

    float x0 = scene.eye.x;
    float y0 = scene.eye.y;
    float z0 = scene.eye.z;

    float xD = ray.x;
    float yD = ray.y;
    float zD = ray.z;

    float minimumt = 100000.0;
    float tminus;
    float tplus;

    for (int s = 0; s < scene.spheres.size(); s++)
    {

        float xC = scene.spheres[s].center.x;
        float yC = scene.spheres[s].center.y;
        float zC = scene.spheres[s].center.z;
        float r = scene.spheres[s].radius;

        float a = (xD * xD) + (yD * yD) + (zD * zD);
        float b = 2 * (xD * (x0 - xC) + yD * (y0 - yC) + zD * (z0 - zC));
        float c = pow((x0 - xC), 2) + pow((y0 - yC), 2) + pow((z0 - zC), 2) - pow(r, 2);

        float discriminant = pow(b, 2) - (4 * a * c);

        if (discriminant > 0.0)
        {
            tminus = (-b - sqrt(discriminant)) / (2 * a);
            tplus = (-b + sqrt(discriminant)) / (2 * a);

            if (tminus<minimumt & tminus> 0.001f)
            {
                hasCollided = true;
                minimumt = tminus;
                closestSphere = scene.spheres[s];
            }
            if (tplus<minimumt & tplus> 0.001f)
            {
                hasCollided = true;
                minimumt = tplus;
                closestSphere = scene.spheres[s];
            }
        }
    }
    if (hasCollided)
    {
        color = Shade_Ray(closestSphere, minimumt, ray, scene);
        vec3 intersectionPoint = add(scene.eye, scale(minimumt, ray));

        //reflective portion
        vec3 N = normalize(sub(intersectionPoint, closestSphere.center));
        //get Reflection Ray
        vec3 reflectedRay = calculateReflectionRay(scene.eye, intersectionPoint, N);
        //get I
        vec3 incoming_Ray = normalize(sub(intersectionPoint, scene.eye));
        vec3 I = scale(-1, incoming_Ray);
        //calculate fresnel affect
        float fresnelAffect = calculateFresnelCoeffSpecular(closestSphere, N, I);
        vec3 intermediateColor = Recursive_Trace_Ray(reflectedRay, intersectionPoint, scene, 0);
        intermediateColor = scale(fresnelAffect, intermediateColor);

        if (color.x < 1)
        {
            if (color.y < 1)
            {
                if (color.z < 1)
                {
                    color = add(color, intermediateColor);
                }
            }
        }

        //transparent portion
        //assumption no overlapping spheres and the rays only go through spheres and air ,mentioned in class
        //calculate fresnel affect
        float fresnelAffect2 = calculateFresnelCoeffTransparent(1.0,closestSphere.material.index_of_reflection,I,N);
        vec3 T = calculateTransparentRay(scene.eye,intersectionPoint,N,1.0,closestSphere.material.index_of_reflection);
        vec3 origin = intersectionPoint;
        vec3 intersectionPoint2 = add(intersectionPoint,T);
        vec3 surface_normal = scale(-1,normalize(sub(intersectionPoint2, closestSphere.center)));
        T = calculateTransparentRay(origin,intersectionPoint2,surface_normal,closestSphere.material.index_of_reflection,1.0);
        vec3 intermediateColor2 = Recursive_Trace_Ray(T, intersectionPoint2, scene, 0);
        intermediateColor2 = scale(fresnelAffect2, intermediateColor2);

        if (color.x < 1)
        {
            if (color.y < 1)
            {
                if (color.z < 1)
                {
                    color = add(color, intermediateColor2);
                }
            }
        }
    }
    return color;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "No file provided!";
        exit(-1);
    }

    std::string filename = argv[1];
    Scene scene = parse_file(filename);
    //also changed this too probs doesnt have anything to do with it
    vec3 imagePPM[(int)scene.height][(int)scene.width];

    //debug();

    vec3 u = normalize(cross(scene.viewdir, scene.updir));
    vec3 v = normalize(cross(u, scene.viewdir));
    vec3 n = scene.viewdir;
    float aspect_ratio = scene.width / scene.height;
    float d = 5.0;
    float he = 2.0 * d * tan(DegreesToRadians(scene.vfov / 2.0));
    float wi = he * aspect_ratio;

    //4 corners of the viewing window

    vec3 ul = add(sub(add(scene.eye, scale(d, normalize(n))), scale((wi / 2.0), u)), scale(he / 2.0, v));
    vec3 ur = add(add(add(scene.eye, scale(d, normalize(n))), scale((wi / 2.0), u)), scale(he / 2.0, v));
    vec3 ll = sub(sub(add(scene.eye, scale(d, normalize(n))), scale((wi / 2.0), u)), scale(he / 2.0, v));
    vec3 lr = sub(add(add(scene.eye, scale(d, normalize(n))), scale((wi / 2.0), u)), scale(he / 2.0, v));

    vec3 deltaH = scale(1 / (2 * scene.width), sub(ur, ul));
    vec3 deltaV = scale(1 / (2 * scene.height), sub(ll, ul));

    for (int i = 0; i < scene.height; i++)
    {
        for (int j = 0; j < scene.width; j++)
        {
            vec3 dividedURUL = scale(1 / scene.width, sub(ur, ul));
            vec3 dividedLLUL = scale(1 / scene.height, sub(ll, ul));
            vec3 index = add(deltaH, deltaV);
            vec3 pointOnViewWindow = add(add(add(ul, scale(i, dividedURUL)), scale(j, dividedLLUL)), index);
            vec3 raydir = normalize(sub(pointOnViewWindow, scene.eye));
            vec3 color = Trace_Ray(raydir, scene);
            //weird fix on this line?? to stop swap of images switched i and j.
            imagePPM[j][i] = {
                .x = color.x,
                .y = color.y,
                .z = color.z,
            };
        }
    }
    std::string outputfile = "default.ppm";
    std::ofstream output_stream(outputfile, std::ios::out | std::ios::binary);
    output_stream << "P3\n"
                  << scene.width << "\n"
                  << scene.height << "\n"
                  << 255 << "\n";
    for (int i = 0; i < scene.width; i++)
    {
        for (int j = 0; j < scene.height; j++)
        {

            output_stream << (imagePPM[i][j].x * 255) << " ";
            output_stream << (imagePPM[i][j].y * 255) << " ";
            output_stream << (imagePPM[i][j].z * 255) << " "
                          << "\n";
        }
        output_stream << "\n";
    }
    output_stream.close();
}
