/*
    A ray-casting renderer based on code found in chapter 15
    of "Computer Graphics: Principles and Practice (3rd Edition)"
    by Hughes, van Dam, McGuire, Sklar, Foley, Feiner, Akeley.
    
    Chapter 15 is available at http://cgpp.net/file/ppcg3e_ch15.pdf
    
    This code is written by Alessandro Gentilini, September 2013.
*/

#include <limits>
#include <cmath>

#undef INFINITY
#define INFINITY (std::numeric_limits<float>::infinity())

#define PI 3.1415926f

template <class T>
T square( const T &a)
{
    return a * a;
}

class Vector2
{
public: float x, y;
};

class Vector3
{
public:
    float x, y, z;
    Vector3(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz)
    {}
    Vector3()
        : x(0), y(0), z(0)
    {}

    float length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    Vector3 direction() const
    {
        const float inverse_magnitude = 1 / length();
        return Vector3(x * inverse_magnitude, y * inverse_magnitude, z * inverse_magnitude);
    }

    Vector3 cross(const Vector3 &v) const
    {
        // Formula 2 in
        // Weisstein, Eric W. "Cross Product." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/CrossProduct.html
        Vector3 u_cross_v;
        const Vector3 &u = *this;
        u_cross_v.x = u.y * v.z - u.z * v.y;
        u_cross_v.y = u.z * v.x - u.x * v.z;
        u_cross_v.z = u.x * v.y - u.y * v.x;
        return u_cross_v;
    }

    float dot(const Vector3 &A) const
    {
        // Formula 7 in
        // Weisstein, Eric W. "Dot Product." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/DotProduct.html
        const Vector3 &B = *this;
        return A.x * B.x + A.y * B.y + A.z * B.z;
    }
};

Vector3 operator-(const Vector3 &a, const Vector3 &b)
{
    Vector3 r(a);
    r.x -= b.x;
    r.y -= b.y;
    r.z -= b.z;
    return r;
}

Vector3 operator+(const Vector3 &a, const Vector3 &b)
{
    Vector3 r(a);
    r.x += b.x;
    r.y += b.y;
    r.z += b.z;
    return r;
}

Vector3 operator*(const Vector3 &v, float a)
{
    return Vector3(v.x * a, v.y * a, v.z * a);
}

Vector3 operator/(const Vector3 &v, float a)
{
    return Vector3(v.x / a, v.y / a, v.z / a);
}

Vector3 operator-(const Vector3 &v)
{
    return Vector3(-v.x, -v.y, -v.z);
}

float dot(const Vector3 &a, const Vector3 &b)
{
    return a.dot(b);
}

typedef Vector2 Point2;

typedef Vector3 Point3;

class Color3
{
public:
    float r, g, b;
    Color3()
        : r(0), g(0), b(0)
    {}
    Color3(float rr, float gg, float bb)
        : r(rr), g(gg), b(bb)
    {}
};

Color3 operator/(const Color3 &c, float a)
{
    return Color3(c.r / a, c.g / a, c.b / a);
}

Color3 operator*(const Color3 &c, float a)
{
    return Color3(c.r * a, c.g * a, c.b * a);
}

Color3 operator+(const Color3 &c1, const Color3 &c2)
{
    return Color3(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b);
}



typedef Color3 Radiance3;
typedef Color3 Power3;

// ???
Color3 operator*(const Radiance3 &r, const Color3 &c)
{
    return Color3(r.r * c.r, r.g * c.g, r.b * c.b);
}

class Ray
{
private:
    Point3 m_origin;
    Vector3 m_direction;
public:
    Ray(const Point3 &org, const Vector3 &dir) :
        m_origin(org), m_direction(dir) {}
    const Point3 &origin() const
    {
        return m_origin;
    }
    const Vector3 &direction() const
    {
        return m_direction;
    }
};

#include <vector>
#include <string>
class Image
{
private:
    int    m_width;
    int    m_height;
    std::vector<Radiance3> m_data;
    int PPMGammaEncode(float radiance, float displayConstant) const;
public:
    Image(int width, int height) :
        m_width(width), m_height(height), m_data(width *height) {}
    int width() const
    {
        return m_width;
    }
    int height() const
    {
        return m_height;
    }
    void set(int x, int y, const Radiance3 &value)
    {
        m_data[x + y * m_width] = value;
    }
    const Radiance3 &get(int x, int y) const
    {
        return m_data[x + y * m_width];
    }
    void save(const std::string &filename, float displayConstant = 15.0f) const;
};

//#include <cmath>
int Image::PPMGammaEncode(float radiance, float d) const
{
    return int(pow(std::min(1.0f, std::max(0.0f, radiance * d)),
                   1.0f / 2.2f) * 255.0f);
}

void Image::save(const std::string &filename, float d) const
{
    FILE *file = fopen(filename.c_str(), "wt");

    fprintf(file, "P3 %d %d 255\n", m_width, m_height);

    for (int y = 0; y < m_height; ++y)
    {
        fprintf(file, "\n# y = %d\n", y);
        for (int x = 0; x < m_width; ++x)
        {
            const Radiance3 &c(get(x, y));
            fprintf(file, "%d %d %d\n",
                    PPMGammaEncode(c.r, d),
                    PPMGammaEncode(c.g, d),
                    PPMGammaEncode(c.b, d));
        }
    }
    fclose(file);
}

class BSDF
{
public:
    Color3 k_L;
    Color3 k_G;
    float s;

    BSDF()
        : k_L(Color3(0.0f, 0.0f, 0.8f)), k_G(Color3(0.2f, 0.2f, 0.2f)), s(100.0f)
    {}


    // changed
    //Vector3 n;

    // changed
    Color3 evaluateFiniteScatteringDensity( const Vector3 &n,
                                            const Vector3 &w_i,
                                            const Vector3 &w_o) const
    {

        const Vector3 &w_h = (w_i + w_o).direction();

        return

            (k_L + k_G * ((s + 8.0f) *

                          powf(std::max(0.0f, w_h.dot(n)), s) / 8.0f)) /

            PI;

    }
};



class Triangle
{
private:
    Point3    m_vertex[3];
    Vector3    m_normal[3];
    BSDF    m_bsdf;
public:

    Triangle(const Point3 &v1, const Point3 &v2, const Point3 &v3,
             const Vector3 &n1, const Vector3 &n2, const Vector3 &n3
            )
    {
        m_vertex[0] = v1;
        m_vertex[1] = v2;
        m_vertex[2] = v3;
        m_normal[0] = n1.direction();
        m_normal[1] = n2.direction();
        m_normal[2] = n3.direction();
    }

    Triangle(const Point3 &v1, const Point3 &v2, const Point3 &v3        )
    {
        m_vertex[0] = v1;
        m_vertex[1] = v2;
        m_vertex[2] = v3;

        m_normal[0] = ((v2 - v1).cross(v3 - v1)).direction();
        m_normal[1] = ((v3 - v2).cross(v1 - v2)).direction();
        m_normal[2] = ((v1 - v3).cross(v2 - v3)).direction();
    }

    const Point3 &vertex(int i) const
    {
        return m_vertex[i];
    }

    const Vector3 &normal(int i) const
    {
        return m_normal[i];
    }

    const BSDF &bsdf() const
    {
        return m_bsdf;
    }
};

class Light
{
public:
    Point3 position;
    /** Over the entire sphere. */
    Power3 power;
};

class Scene
{
public:
    std::vector<Triangle> triangleArray;
    std::vector<Light>    lightArray;
};

class Camera
{
public:
    float zNear;
    float zFar;
    float fieldOfViewX;
    Camera() : zNear(-0.1f), zFar(-100.0f), fieldOfViewX(PI / 2.0f) {}
};

float intersect(const Ray &R, const Triangle &T, float weight[3])
{
    const Vector3 &e1 = T.vertex(1) - T.vertex(0);
    const Vector3 &e2 = T.vertex(2) - T.vertex(0);
    const Vector3 &q = R.direction().cross(e2);
    const float a = e1.dot(q);
    const Vector3 &s = R.origin() - T.vertex(0);
    const Vector3 &r = s.cross(e1);
    // Barycentric vertex weights
    weight[1] = s.dot(q) / a;
    weight[2] = R.direction().dot(r) / a;
    weight[0] = 1.0f - (weight[1] + weight[2]);

    const float dist = e2.dot(r) / a;
    static const float epsilon = 1e-7f;
    static const float epsilon2 = 1e-10;
    if ((a <= epsilon) || (weight[0] < -epsilon2) ||
            (weight[1] < -epsilon2) || (weight[2] < -epsilon2) ||
            (dist <= 0.0f))
    {
        // The ray is nearly parallel to the triangle, or the
        // intersection lies outside the triangle or behind
        // the ray origin: "infinite" distance until intersection.
        return INFINITY;
    }
    else
    {
        return dist;
    }
}

bool visible(const Vector3 &P, const Vector3 &direction, float
             distance, const Scene &scene)
{
    static const float rayBumpEpsilon = 1e-4;

    const Ray shadowRay(P + direction * rayBumpEpsilon, direction);

    distance -= rayBumpEpsilon;

    // Test each potential shadow caster to see if it lies between P and the light
    float ignore[3];
    for (unsigned int s = 0; s < scene.triangleArray.size(); ++s)
    {
        if (intersect(shadowRay, scene.triangleArray[s], ignore) < distance)
        {

            // This triangle is closer than the light

            return false;

        }
    }

    return true;
}


void shade(const Scene &scene, const Triangle &T, const Point3 &P,
           const Vector3 &n, const Vector3 &w_o, Radiance3 &L_o)
{
    L_o = Color3(0.0f, 0.0f, 0.0f);
    // For each direction (to a light source)
    for (unsigned int i = 0; i < scene.lightArray.size(); ++i)
    {
        const Light &light = scene.lightArray[i];
        const Vector3 &offset = light.position - P;
        const float distanceToLight = offset.length();

        const Vector3 &w_i = offset / distanceToLight;
        if (visible(P, w_i, distanceToLight, scene))
        {
            const Radiance3 &L_i = light.power / (4 * PI * square(distanceToLight));
            // Scatter the light
            L_o = L_o + //changed
                  L_i *
                  // changed
                  T.bsdf().evaluateFiniteScatteringDensity(n, w_i, w_o) *
                  std::max(0.0f, dot(w_i, n));
        }
    }
}

bool sampleRayTriangle(const Scene &scene, int x, int y,
                       const Ray &R, const Triangle &T,
                       Radiance3 &radiance, float &distance);

bool sampleRayTriangle(const Scene &scene, int x, int y, const Ray &R,
                       const Triangle &T, Radiance3 &radiance, float &distance)
{
    float weight[3];
    const float d = intersect(R, T, weight);
    if (d >= distance)
    {
        return false;
    }
    // This intersection is closer than the previous one
    distance = d;
    // Intersection point
    const Point3 &P = R.origin() + R.direction() * d;
    // Find the interpolated vertex normal at the intersection
    const Vector3 &n = (T.normal(0) * weight[0] +
                        T.normal(1) * weight[1] +
                        T.normal(2) * weight[2]).direction();
    const Vector3 &w_o = -R.direction();
    shade(scene, T, P, n, w_o, radiance);
    // Debugging intersect: set to white on any intersection
    //radiance = Radiance3(1, 1, 1);
    // Debugging barycentric
    //radiance = Radiance3(weight[0], weight[1], weight[2]) / 15;
    return true;
}

Ray computeEyeRay(float x, float y, int width, int height, const Camera &camera)
{
    const float aspect = float(height) / width;
    // Compute the side of a square at z = -1 based on our
    // horizontal left-edge-to-right-edge field of view
    const float s = -2.0f * tan(camera.fieldOfViewX * 0.5f);
    const Vector3 &start =
        Vector3( (x / width - 0.5f) * s, -(y / height - 0.5f) * s * aspect, 1.0f) * camera.zNear;

    return Ray(start, start.direction());
}

#include <iostream>

/** Trace eye rays with origins in the box from [x0, y0] to (x1, y1).*/
void rayTrace(Image &image, const Scene &scene,
              const Camera &camera, int x0, int x1, int y0, int y1)
{
    // For each pixel
    for (int y = y0; y < y1; ++y)
    {
        for (int x = y0; x < x1; ++x)
        {
            // Ray through the pixel
            const Ray &R = computeEyeRay(x + 0.5f, y + 0.5f, image.width(),
                                         image.height(), camera);
            // Distance to closest known intersection
            float distance = INFINITY;
            Radiance3 L_o;
            // For each triangle
            for (unsigned int t = 0; t < scene.triangleArray.size(); ++t)
            {
                const Triangle &T = scene.triangleArray[t];
                if (sampleRayTriangle(scene, x, y, R, T, L_o, distance))
                {
                    image.set(x, y, L_o);
                }
            }
            if ( y % 10 == 0 && x % 10 == 0)
            {
                std::cout << y << "," << x << "\n";
            }
        }
    }
}

#include <fstream>


int main(int, char **)
{
    {
        Image image(800, 500);
        Scene scene;
        Triangle t(Point3(0, 1, -2), Point3(-1.9, -1, -2), Point3(1.6, -0.5, -2),
                   Vector3( 0.0f, 0.6f, 1.0f), Vector3(-0.4f, -0.4f, 1.0f), Vector3( 0.4f, -0.4f, 1.0f));
        scene.triangleArray.push_back(t);
        Light light;
        light.position = Point3(1.0f, 3.0f, 1.0f);
        light.power = Power3(10, 10, 10);
        scene.lightArray.push_back(light);
        Camera camera;  
        rayTrace(image, scene, camera, 0, image.width(), 0, image.height());
        image.save("first_light.ppm", 1);
    }

    {
        Image image(800, 500);
        Scene scene;

        Light light;
        light.position = Point3(3.0f, .0f, .0f);
        light.power = Power3(1000, 1000, 1000);
        scene.lightArray.push_back(light);
        Camera camera;

        std::ifstream v("./vertices.txt");
        std::vector< Point3 > vertices;
        // dummy:
        vertices.push_back(Point3());
        float x, y, z;
        Point3 centroid;
        float max_x = -INFINITY;
        float max_y = -INFINITY;
        float max_z = -INFINITY;
        float min_x = INFINITY;
        float min_y = INFINITY;
        float min_z = INFINITY;        
        while (v >> x >> y >> z)
        {
            vertices.push_back(Point3(x, y-50, z-100));
            centroid = centroid + Point3(x, y, z);
            max_x = std::max(max_x,x);
            max_y = std::max(max_y,y);
            max_z = std::max(max_z,z);
            min_x = std::min(min_x,x);
            min_y = std::min(min_y,y);
            min_z = std::min(min_z,z);
        }
        centroid = centroid/vertices.size();
        std::cerr << "centroid=" << centroid.x << " " << centroid.y << " " << centroid.z << "\n";
        std::cerr << "x=" << min_x << ".." << max_x << "\n";
        std::cerr << "y=" << min_y << ".." << max_y << "\n";
        std::cerr << "z=" << min_z << ".." << max_z << "\n";

        std::ifstream t("./triangles.txt");
        int i1, i2, i3;
        const size_t sz = vertices.size();
        size_t cnt = 0;
        while ( t >> i1 >> i2 >> i3)
        {
            if (i1 >= sz || i2 >= sz || i3 >= sz)
            {
                std::cerr << "out of range indexes ("<< sz <<"): " << i1 << " " << i2 << " " << i3 << "\n";
            }
            else
            {
                if ( cnt++ % 1 == 0 )
                {
                    scene.triangleArray.push_back(Triangle(vertices[i1], vertices[i2], vertices[i3]));
                }
            }
        }

        rayTrace(image, scene, camera, 0, image.width(), 0, image.height());
        image.save("teapot.ppm", 1);
    }

    return 0;
}
