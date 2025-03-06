#include<bits/stdc++.h>
using namespace std;



class Point
{

public:
    double x, y, z;

    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }

    Point operator+(Point b)
    {
        return Point(x + b.x, y + b.y, z + b.z);
    }
    Point operator-(Point b)
    {
        return Point(x - b.x, y - b.y, z - b.z);
    }
    Point operator-(double b)
    {
        return Point(x - b, y - b, z - b);
    }
    Point operator*(double b)
    {
        return Point(x * b, y * b, z * b);
    }
    Point operator/(double b)
    {
        return Point(x / b, y / b, z / b);
    }

};




Point normalized_point(Point p)
{
    double length = p.x * p.x + p.y * p.y + p.z * p.z;

    length = sqrt(length);

    Point result;

    result.x = p.x / length;
    result.y = p.y / length;
    result.z = p.z / length;

    return result;

}
Point cross_product(Point p1, Point p2)
{
    return Point(p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x);

}
double dot_product(Point p1, Point p2)
{
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
double getLength(Point p)
{
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

class Ray{
    public:
    Point start;
    Point dir;
};


void drawSphere(double radius,int slices,int stacks)
{
    struct Point points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0; i<=stacks; i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        for(j=0; j<=slices; j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0; i<stacks; i++)
    {
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);

        for(j=0; j<slices; j++)
        {

            glBegin(GL_QUADS);
            {

                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            }
            glEnd();
        }
    }
}

class PointLight{
    public:
        Point light_pos;
        double color[3];

    void draw()
    {
        glPushMatrix();
        glPointSize(5);
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_POINTS);
        glVertex3f(light_pos.x,light_pos.y,light_pos.z);
        glEnd();
        glPopMatrix();

    }

};
class SpotLight
{
public:
    PointLight point_light;
    Point light_direction;
    double cutoff_angle;

    void draw()
    {
        glPushMatrix();
        glPointSize(10);
        glColor3f(point_light.color[0],point_light.color[1],point_light.color[2]);
        glBegin(GL_POINTS);
        glVertex3f(point_light.light_pos.x,point_light.light_pos.y,point_light.light_pos.z);
        glEnd();
        glPopMatrix();

    }

};

class Object{
public:
    Point reference_point; // should have x, y, z
    double height, width, length;
    double color[3];
    double coEfficients[4]; // ambient, diffuse, specular, reflection coeff
    int shine; // exponent term of specular component
    Object() {}
    virtual void draw() {

    }
    virtual double intersect(Ray ray,double* color,int level) {
        return -1.0;
    }
    void setColor() {}
    void setShine() {}
    void setCoEfficients() {}
};

extern vector<Object *> objects;
extern vector<PointLight *> pointLights;
extern vector<SpotLight *> spotLights;
extern int recursion_level;
double ephsilon = 0.0001;

class Sphere : public Object
{
public:
    Sphere(Point center,double radius)
    {
        reference_point = center;
        length = radius;
    }
    void draw()
    {
        // write codes for drawing sphere
        glColor3f(color[0],color[1],color[2]);
        glPushMatrix();
        glTranslatef(reference_point.x,reference_point.y,reference_point.z);
        drawSphere(length,30,20);
        glPopMatrix();
    }
    void getColorAt(Point p,double* color)
    {
        color[0] = this->color[0];
        color[1] = this->color[1];
        color[2] = this->color[2];
    }
    double intersect(Ray ray,double* color,int level)
    {
        //code for finding tmin_start
        ray.start = ray.start - reference_point; // if sphere were at origin then ray was at,

        double a = 1;
        double b = 2 * dot_product(ray.dir,ray.start);
        double c = dot_product(ray.start,ray.start) - length*length;

        double d = b*b - 4 * a * c;


        if(d < 0) return -1;

        d = sqrt(d);

        double t1 = (-b + d) /(2 * a);
        double t2 = (-b - d) /(2 * a);


        double t = -1;

        if(t1 > 0 && t2 > 0)
        {
            t = min(t1,t2);
        }
        else if(t1 > 0)
        {
            t = t1;
        }
        else if(t2 > 0)
        {
            t = t2;
        }
        else
        {
            t = -1;
        }
        //code for finding tmin_end
        if(t>0)
        {
            if(level == 0) return t;


            ray.start = ray.start + reference_point;

            Point intersection_point = ray.start + ray.dir * t;

            double intersection_color[3]={-1,-1,-1};

            getColorAt(intersection_point, intersection_color);

            color[0] = intersection_color[0] * coEfficients[0];
            color[1] = intersection_color[1] * coEfficients[0];
            color[2] = intersection_color[2] * coEfficients[0];


            ///calc normal at intersection point
            Point normal_direction = intersection_point - reference_point;

            normal_direction = normalized_point(normal_direction);

            for(int i=0;i<pointLights.size();i++)
            {

                double dist_t = getLength(pointLights[i]->light_pos - intersection_point);

                if(dist_t < ephsilon)
                {
                    continue;
                }
                //construct ray of light
                Ray light_ray;
                light_ray.start = pointLights[i]->light_pos;
                light_ray.dir = normalized_point(intersection_point - pointLights[i]->light_pos);

                bool shadow = false;
                for(int i=0;i<objects.size();i++)
                {
                    double color[3];
                    double temp_t = objects[i]->intersect(light_ray,color,0);

                    if(temp_t != -1 && temp_t + ephsilon < dist_t)
                    {
                        shadow = true;
                        break;
                    }
                }

                if(shadow == false)
                {
                    Point s = light_ray.dir;
                    Point v = ray.dir;
                    Point m = normal_direction;  ///normal

                    s = normalized_point(s);

                    m = normalized_point(m);

                    double lambart_val = max(0.0,dot_product(s*-1,m));

                    //printf("%lf\n",lambart_val);

                    double temp = 2 * dot_product(s,m);
                    Point reflection_dir = s - normal_direction * temp;
                    reflection_dir = normalized_point(reflection_dir);

                    double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));

                  //  printf("before: %lf %lf %lf\n",dist_t,light_ray.dir.y,light_ray.dir.z);

                    color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                    color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];
                   // printf("after: %lf %lf %lf\n",color[0],color[1],color[2]);
                }
            }
            for(int i=0;i<spotLights.size();i++)
            {
                Ray spot_point_light_ray;
                spot_point_light_ray.start = spotLights[i]->point_light.light_pos;
                spot_point_light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);

                Point spot_light_dir = spotLights[i] -> light_direction;
                spot_light_dir = normalized_point(spot_light_dir);

                double cos_theta = dot_product(spot_light_dir,spot_point_light_ray.dir);
                double theta = acos(cos_theta) * 180.0 / acos(-1);

                if(theta < spotLights[i]->cutoff_angle)
                {

                    double dist_t = getLength(spotLights[i]->point_light.light_pos - intersection_point);

                    if(dist_t < ephsilon)
                    {
                        continue;
                    }
                    Ray light_ray;
                    light_ray.start = spotLights[i]->point_light.light_pos;
                    light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);


                    ///calc normal
                    Point normal_direction = intersection_point - reference_point;
                    normal_direction = normalized_point(normal_direction);

                    bool shadow = false;
                    for(int i=0;i<objects.size();i++)
                    {
                        double color[3];
                        double temp_t = objects[i]->intersect(light_ray,color,0);

                        if(temp_t != -1 && temp_t + ephsilon < dist_t)
                        {
                            shadow = true;
                            break;
                        }
                    }


                    if(shadow == false)
                    {
                        Point s = light_ray.dir;
                        Point v = ray.dir;
                        Point m = normal_direction;

                        s = normalized_point(s);
                        m = normalized_point(m);

                        double lambart_val = max(0.0,dot_product(s*-1,m));

                        double temp = 2 * dot_product(s,m);
                        Point reflection_dir = s - normal_direction * temp;
                        reflection_dir = normalized_point(reflection_dir);

                        double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));



                        color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                        color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];


                    }
                }
            }



            if(level >= recursion_level)
                return t;

            //construct reflected ray
            Point s = ray.dir;

            Point m = normal_direction;

            s = normalized_point(s);
            m = normalized_point(m);

            double temp = 2 * dot_product(s,m);
            Point reflection_dir = s - normal_direction * temp;
            reflection_dir = normalized_point(reflection_dir);

            // actually slightly forward from the point (by moving the
            // start a little bit towards the reflection direction)
            // to avoid self intersection

            Ray desired_reflected_ray;
            desired_reflected_ray.start = intersection_point + reflection_dir * ephsilon;
            desired_reflected_ray.dir = reflection_dir;

            //printf("%lf %lf %lf\n",desired_reflected_ray.start.x,desired_reflected_ray.start.y,desired_reflected_ray.start.z);
            double tmpcolor[3];
			double mini_t = INT_MAX;
			int index = -1;
			for(int k=0;k<objects.size();k++)
			{
				double t = objects[k]->intersect(desired_reflected_ray,tmpcolor, 0);

                if(t > 0 && t < mini_t)
                {
                    mini_t = t;
                    index = k;
                }

			}
			// if nearest object is found, then shade the pixel
			if(index != -1)
			{

				double getcolor[3];

				double t = objects[index]->intersect(desired_reflected_ray,getcolor,level+1);

				color[0] += getcolor[0] * coEfficients[3];
                color[1] += getcolor[1] * coEfficients[3];
                color[2] += getcolor[2] * coEfficients[3];

               // printf("%lf %lf %lf\n",color[0],color[1],color[2]);
			}

            return t;

        }
        else
        {
            return -1.0;
        }

    }
};

double getDet(double matrix[3][3])
{
    double sum = 0;

    sum = sum + matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
    sum = sum - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]);
    sum = sum + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    return sum;
}

class Triangle : public Object
{

public:
    Point p1;
    Point p2;
    Point p3;
    Triangle(Point p1,Point p2, Point p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }
    void draw()
    {

        glColor3f(color[0],color[1],color[2]);
        glPushMatrix();

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(p1.x,p1.y,p1.z);
            glVertex3f(p2.x,p2.y,p2.z);
            glVertex3f(p3.x,p3.y,p3.z);
        }
        glEnd();
        glPopMatrix();
    }
    void getColorAt(Point p,double* color)
    {
        color[0] = this->color[0];
        color[1] = this->color[1];
        color[2] = this->color[2];
    }
    double intersect(Ray ray,double* color,int level)
    {
        Point aa = this->p3 - this->p1;
        Point bb = this->p2 - this->p1;

        Point normal_direction = cross_product(bb,aa);
        normal_direction = normalized_point(normal_direction);

        if(dot_product(normal_direction,ray.dir) < 0)
        {
            normal_direction = normal_direction * -1;
        }

        Point a = this->p1;
        Point b = this->p2;
        Point c = this->p3;

        double mat_A[3][3]={
            {
                a.x-b.x,
                a.x-c.x,
                ray.dir.x
            },
            {
                a.y-b.y,
                a.y-c.y,
                ray.dir.y
            },
            {
                a.z-b.z,
                a.z-c.z,
                ray.dir.z
            }

        };

        double mat_beta_calc[3][3]={
            {
                a.x-ray.start.x,
                a.x-c.x,
                ray.dir.x
            },
            {
                a.y-ray.start.y,
                a.y-c.y,
                ray.dir.y
            },
            {
                a.z-ray.start.z,
                a.z-c.z,
                ray.dir.z
            }

        };
        double mat_gamma_calc[3][3]={
            {
                a.x-b.x,
                a.x-ray.start.x,
                ray.dir.x
            },
            {
                a.y-b.y,
                a.y-ray.start.y,
                ray.dir.y
            },
            {
                a.z-b.z,
                a.z-ray.start.z,
                ray.dir.z
            }

        };

        double mat_t_calc[3][3] = {

            {
                a.x-b.x,
                a.x-c.x,
                a.x - ray.start.x
            },
            {
                a.y-b.y,
                a.y-c.y,
                a.y - ray.start.y
            },
            {
                a.z-b.z,
                a.z-c.z,
                a.z - ray.start.z
            }

        };

        double det_a = getDet(mat_A);
        double det_beta = getDet(mat_beta_calc);
        double det_gamma = getDet(mat_gamma_calc);
        double det_t = getDet(mat_t_calc);


        double beta = det_beta / det_a;
        double gamma = det_gamma / det_a;
        double t = det_t / det_a;


        if(beta + gamma < 1 && beta > 0 && gamma > 0 && t>0)
        {

            if(level == 0) return t;


            Point intersection_point = ray.start + ray.dir * t;

            double intersection_color[3]={-1,-1,-1};

            getColorAt(intersection_point, intersection_color);

            color[0] = intersection_color[0] * coEfficients[0];
            color[1] = intersection_color[1] * coEfficients[0];
            color[2] = intersection_color[2] * coEfficients[0];

            ///calc normal at intersection point

            normal_direction = normalized_point(normal_direction);

            for(int i=0;i<pointLights.size();i++)
            {

                double dist_t = getLength(pointLights[i]->light_pos - intersection_point);

                if(dist_t < ephsilon)
                {
                    continue;
                }
                //construct ray of light
                Ray light_ray;
                light_ray.start = pointLights[i]->light_pos;
                light_ray.dir = normalized_point(intersection_point - pointLights[i]->light_pos);

                bool shadow = false;
                for(int i=0;i<objects.size();i++)
                {
                    double color[3];
                    double temp_t = objects[i]->intersect(light_ray,color,0);
                    //printf("%lf\n",temp_t);
                    if(temp_t != -1 && temp_t + ephsilon < dist_t)
                    {
                        shadow = true;
                        break;
                    }
                }

                if(shadow == false)
                {
                    Point s = light_ray.dir;
                    Point v = ray.dir;
                    Point m = normal_direction;  ///normal

                    s = normalized_point(s);

                    m = normalized_point(m);

                    double lambart_val = max(0.0,dot_product(s*-1,m));



                    double temp = 2 * dot_product(s,m);
                    Point reflection_dir = s - normal_direction * temp;
                    reflection_dir = normalized_point(reflection_dir);

                    double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));


                    color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                    color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];

                }
            }
            for(int i=0;i<spotLights.size();i++)
            {
                Ray spot_point_light_ray;
                spot_point_light_ray.start = spotLights[i]->point_light.light_pos;
                spot_point_light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);

                Point spot_light_dir = spotLights[i] -> light_direction;
                spot_light_dir = normalized_point(spot_light_dir);

                double cos_theta = dot_product(spot_light_dir,spot_point_light_ray.dir);
                double theta = acos(cos_theta) * 180.0 / acos(-1);

                if(theta < spotLights[i]->cutoff_angle)
                {

                    double dist_t = getLength(spotLights[i]->point_light.light_pos - intersection_point);

                    if(dist_t < ephsilon)
                    {
                        continue;
                    }

                    Ray light_ray;
                    light_ray.start = spotLights[i]->point_light.light_pos;
                    light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);


                    ///calc normal
                    //Point normal_direction = intersection_point - reference_point;
                    normal_direction = normalized_point(normal_direction);

                    bool shadow = false;
                    for(int i=0;i<objects.size();i++)
                    {
                        double color[3];
                        double temp_t = objects[i]->intersect(light_ray,color,0);

                        if(temp_t != -1 && temp_t + ephsilon < dist_t)
                        {
                            shadow = true;
                            break;
                        }
                    }
                    if(shadow == false)
                    {
                        Point s = light_ray.dir;
                        Point v = ray.dir;
                        Point m = normal_direction;

                        s = normalized_point(s);
                        m = normalized_point(m);

                        double lambart_val = max(0.0,dot_product(s*-1,m));


                        double temp = 2 * dot_product(s,m);
                        Point reflection_dir = s - normal_direction * temp;
                        reflection_dir = normalized_point(reflection_dir);

                        double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));



                        color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                        color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];

                        //printf("after: %lf %lf %lf\n",intersection_color[0],intersection_color[1],intersection_color[2]);

                    }
                }
            }
            if(level >= recursion_level)
                return t;

            //construct reflected ray
            Point s = ray.dir;

            Point m = normal_direction;

            s = normalized_point(s);
            m = normalized_point(m);

            double temp = 2 * dot_product(s,m);
            Point reflection_dir = s - normal_direction * temp;
            reflection_dir = normalized_point(reflection_dir);

            // actually slightly forward from the point (by moving the
            // start a little bit towards the reflection direction)
            // to avoid self intersection

            Ray desired_reflected_ray;
            desired_reflected_ray.start = intersection_point + reflection_dir * ephsilon;
            desired_reflected_ray.dir = reflection_dir;

            //printf("%lf %lf %lf\n",desired_reflected_ray.start.x,desired_reflected_ray.start.y,desired_reflected_ray.start.z);
            double tmpcolor[3];
			double mini_t = INT_MAX;
			int index = -1;
			for(int k=0;k<objects.size();k++)
			{
				double t = objects[k]->intersect(desired_reflected_ray,tmpcolor, 0);

                if(t > 0 && t < mini_t)
                {
                    mini_t = t;
                    index = k;
                }

			}
			// if nearest object is found, then shade the pixel
			if(index != -1)
			{

				double getcolor[3];

				double t = objects[index]->intersect(desired_reflected_ray,getcolor,level+1);

				color[0] += getcolor[0] * coEfficients[3];
                color[1] += getcolor[1] * coEfficients[3];
                color[2] += getcolor[2] * coEfficients[3];

               // printf("%lf %lf %lf\n",color[0],color[1],color[2]);
			}

            return t;

        }
        else
        {
            return -1.0;
        }

    }
};

class Floor : public Object
{

public:

    Floor(double floorWidth,double tileWidth)
    {
        reference_point = Point(-floorWidth/2.0,-floorWidth/2.0,0);
        length = tileWidth;

        this->width = floorWidth;
        this->length = tileWidth;
    }
    void draw()
    {
        int ci = 0;
        int cj = 0;
        for(double i=0;i<width;i+=length,ci++)
        {

            for(double j=0;j<width;j+=length,cj++)
            {

                if((ci + cj)%2 == 0)
                {
                    glColor3f(1,1,1);
                }
                else
                {
                    glColor3f(0,0,0);
                }
                glPushMatrix();
                glBegin(GL_QUADS);
                {
                    glVertex3f(reference_point.x+j,reference_point.y+i,0);
                    glVertex3f(reference_point.x+j,reference_point.y+i+length,0);
                    glVertex3f(reference_point.x+j+length,reference_point.y+i+length,0);
                    glVertex3f(reference_point.x+j+length,reference_point.y+i,0);

                }
                glEnd();
                glPopMatrix();
            }
        }
    }
    void getColorAt(Point p,double* color)
    {
        int ci = (p.x - reference_point.x) / length;
        int cj = (p.y - reference_point.y) / length;


        if((ci+cj)%2 == 0)
        {
            color[0] = 1;
            color[1] = 1;
            color[2] = 1;
            return;

        }
        else
        {
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
            return;
        }
    }

    double intersect(Ray ray,double* color,int level)
    {
        Point normal_direction;

        if(ray.dir.z > 0)
            normal_direction = Point(0, 0, 1);
        else
            normal_direction = Point(0, 0, -1);

        double dot_pro = dot_product(normal_direction , ray.dir);

        double t = -1;

        if(dot_pro != 0)
        {
            t = -dot_product(normal_direction , ray.start) / dot_pro;
        }

        if(t < 0) return -1;

        Point p = ray.start + ray.dir * t;

        if(p.x < reference_point.x || p.x > reference_point.x + width || p.y < reference_point.y || p.y > reference_point.y + width){
            return -1;
        }

        if(t>0)
        {

            if(level == 0) return t;


            Point intersection_point = ray.start + ray.dir * t;

            double intersection_color[3]={-1,-1,-1};

            getColorAt(intersection_point, intersection_color);

            color[0] = intersection_color[0] * coEfficients[0];
            color[1] = intersection_color[1] * coEfficients[0];
            color[2] = intersection_color[2] * coEfficients[0];

            ///calc normal at intersection point

            normal_direction = normalized_point(normal_direction);

            for(int i=0;i<pointLights.size();i++)
            {

                double dist_t = getLength(pointLights[i]->light_pos - intersection_point);

                if(dist_t < ephsilon)
                {
                    continue;
                }
                //construct ray of light
                Ray light_ray;
                light_ray.start = pointLights[i]->light_pos;
                light_ray.dir = normalized_point(intersection_point - pointLights[i]->light_pos);

                bool shadow = false;
                for(int i=0;i<objects.size();i++)
                {
                    double color[3];
                    double temp_t = objects[i]->intersect(light_ray,color,0);
                    //printf("%lf\n",temp_t);
                    if(temp_t != -1 && temp_t + ephsilon < dist_t)
                    {
                        shadow = true;
                        break;
                    }
                }

                if(shadow == false)
                {
                    Point s = light_ray.dir;
                    Point v = ray.dir;
                    Point m = normal_direction;  ///normal

                    s = normalized_point(s);

                    m = normalized_point(m);

                    double lambart_val = max(0.0,dot_product(s*-1,m));



                    double temp = 2 * dot_product(s,m);
                    Point reflection_dir = s - normal_direction * temp;
                    reflection_dir = normalized_point(reflection_dir);

                    double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));


                    color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                    color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];

                }
            }
            for(int i=0;i<spotLights.size();i++)
            {
                Ray spot_point_light_ray;
                spot_point_light_ray.start = spotLights[i]->point_light.light_pos;
                spot_point_light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);

                Point spot_light_dir = spotLights[i] -> light_direction;
                spot_light_dir = normalized_point(spot_light_dir);

                double cos_theta = dot_product(spot_light_dir,spot_point_light_ray.dir);
                double theta = acos(cos_theta) * 180.0 / acos(-1);

                if(theta < spotLights[i]->cutoff_angle)
                {

                    double dist_t = getLength(spotLights[i]->point_light.light_pos - intersection_point);

                    if(dist_t < ephsilon)
                    {
                        continue;
                    }

                    Ray light_ray;
                    light_ray.start = spotLights[i]->point_light.light_pos;
                    light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);


                    ///calc normal
                    //Point normal_direction = intersection_point - reference_point;
                    normal_direction = normalized_point(normal_direction);

                    bool shadow = false;
                    for(int i=0;i<objects.size();i++)
                    {
                        double color[3];
                        double temp_t = objects[i]->intersect(light_ray,color,0);

                        if(temp_t != -1 && temp_t + ephsilon < dist_t)
                        {
                            shadow = true;
                            break;
                        }
                    }
                    if(shadow == false)
                    {
                        Point s = light_ray.dir;
                        Point v = ray.dir;
                        Point m = normal_direction;

                        s = normalized_point(s);
                        m = normalized_point(m);

                        double lambart_val = max(0.0,dot_product(s*-1,m));


                        double temp = 2 * dot_product(s,m);
                        Point reflection_dir = s - normal_direction * temp;
                        reflection_dir = normalized_point(reflection_dir);

                        double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));



                        color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                        color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];

                        //printf("after: %lf %lf %lf\n",intersection_color[0],intersection_color[1],intersection_color[2]);

                    }
                }
            }
            if(level >= recursion_level)
                return t;

            //construct reflected ray
            Point s = ray.dir;

            Point m = normal_direction;

            s = normalized_point(s);
            m = normalized_point(m);

            double temp = 2 * dot_product(s,m);
            Point reflection_dir = s - normal_direction * temp;
            reflection_dir = normalized_point(reflection_dir);

            // actually slightly forward from the point (by moving the
            // start a little bit towards the reflection direction)
            // to avoid self intersection

            Ray desired_reflected_ray;
            desired_reflected_ray.start = intersection_point + reflection_dir * ephsilon;
            desired_reflected_ray.dir = reflection_dir;

            //printf("%lf %lf %lf\n",desired_reflected_ray.start.x,desired_reflected_ray.start.y,desired_reflected_ray.start.z);
            double tmpcolor[3];
			double mini_t = INT_MAX;
			int index = -1;
			for(int k=0;k<objects.size();k++)
			{
				double t = objects[k]->intersect(desired_reflected_ray,tmpcolor, 0);

                if(t > 0 && t < mini_t)
                {
                    mini_t = t;
                    index = k;
                }

			}
			// if nearest object is found, then shade the pixel
			if(index != -1)
			{

				double getcolor[3];

				double t = objects[index]->intersect(desired_reflected_ray,getcolor,level+1);

				color[0] += getcolor[0] * coEfficients[3];
                color[1] += getcolor[1] * coEfficients[3];
                color[2] += getcolor[2] * coEfficients[3];

               // printf("%lf %lf %lf\n",color[0],color[1],color[2]);
			}

            return t;

        }
        else
        {
            return -1.0;
        }

    }
};

class General : public Object
{
    double a,b,c,d,e,f,g,h,i,j;

    public:
    General(double a,double b,double c,
    double d,double e,double f,double g,double h,double i,double j)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        this->e = e;
        this->f = f;
        this->g = g;
        this->h = h;
        this->i = i;
        this->j = j;
    }
    void getColorAt(Point p,double* color)
    {
        color[0] = this->color[0];
        color[1] = this->color[1];
        color[2] = this->color[2];
    }
    void draw()
    {

    }
    bool check_clip(double t,Ray ray)
    {
        Point p = ray.start + ray.dir * t;

        if(abs(length) > ephsilon)
        {
            if(p.x >= reference_point.x && p.x <= reference_point.x + length)
            {
                return true;
            }
            else
                return false;
        }
        if(abs(width) > ephsilon)
        {
            if(p.y >=reference_point.y && p.y <= reference_point.y + width)
            {
                return true;
            }
            else
                return false;
        }
        if(abs(height) > ephsilon)
        {
            if(p.z >=reference_point.z && p.z <= reference_point.z + height)
            {
                return true;
            }
            else
                return false;
        }

        return true;
    }
    double intersect(Ray ray,double* color,int level)
    {
        double R_ox = ray.start.x;
        double R_oy = ray.start.y;
        double R_oz = ray.start.z;

        double R_dx = ray.dir.x;
        double R_dy = ray.dir.y;
        double R_dz = ray.dir.z;

        double A = a*(R_dx * R_dx) + b*(R_dy * R_dy) + c*(R_dz * R_dz) + d * (R_dx * R_dy) + e * (R_dx * R_dz) + f * (R_dy * R_dz);
        double B = 2 * a * R_ox * R_dx + 2 * b * R_oy * R_dy + 2 * c * R_oz * R_dz + d * (R_ox * R_dy + R_oy * R_dx) + e * (R_ox * R_dz + R_oz * R_dx) +
                    f * (R_oy * R_dz + R_oz * R_dy) + g * R_dx + h * R_dy + i * R_dz;

        double C = a * R_ox * R_ox + b * R_oy * R_oy + c * R_oz * R_oz + d * R_ox * R_oy + e * R_ox * R_oz + f * R_oy * R_oz + g * R_ox + h * R_oy + i * R_oz + j;

        double D = B*B - 4 * A * C;

        double t = -1;

        if(D < 0) return -1;
        else if(A < ephsilon)
        {
            t = -C/B;

            if(check_clip(t,ray) == false) return -1;
        }
        else
        {
            D = sqrt(D);

            double t1 = (-B + D)/(2 * A);
            double t2 = (-B - D)/(2 * A);

            if(check_clip(t1,ray) == false)
            {
                t1 = -1;
            }
            if(check_clip(t2,ray) == false)
            {
                t2 = -1;
            }


            if(t1 > 0 && t2 > 0)
            {
                t = min(t1,t2);
            }
            else if(t1 > 0)
            {
                t = t1;
            }
            else if(t2 > 0)
            {
                t = t2;
            }


        }
        if(t>0)
        {


            if(level == 0) return t;


            Point intersection_point = ray.start + ray.dir * t;

            double intersection_color[3]={-1,-1,-1};

            getColorAt(intersection_point, intersection_color);

            color[0] = intersection_color[0] * coEfficients[0];
            color[1] = intersection_color[1] * coEfficients[0];
            color[2] = intersection_color[2] * coEfficients[0];

            ///calc normal at intersection point

            Point normal_direction;
            normal_direction.x = 2 * a * intersection_point.x + d * intersection_point.y + e * intersection_point.z + g;
            normal_direction.y = 2 * b * intersection_point.y + d * intersection_point.x + f * intersection_point.z + h;
            normal_direction.z = 2 * c * intersection_point.z + e * intersection_point.x + f * intersection_point.y + i;

            normal_direction = normalized_point(normal_direction);

            for(int i=0;i<pointLights.size();i++)
            {

                double dist_t = getLength(pointLights[i]->light_pos - intersection_point);

                if(dist_t < ephsilon)
                {
                    continue;
                }
                //construct ray of light
                Ray light_ray;
                light_ray.start = pointLights[i]->light_pos;
                light_ray.dir = normalized_point(intersection_point - pointLights[i]->light_pos);

                bool shadow = false;
                for(int i=0;i<objects.size();i++)
                {
                    double color[3];
                    double temp_t = objects[i]->intersect(light_ray,color,0);
                    //printf("%lf\n",temp_t);
                    if(temp_t != -1 && temp_t + ephsilon < dist_t)
                    {
                        shadow = true;
                        break;
                    }
                }

                if(shadow == false)
                {
                    Point s = light_ray.dir;
                    Point v = ray.dir;
                    Point m = normal_direction;  ///normal

                    s = normalized_point(s);

                    m = normalized_point(m);

                    double lambart_val = max(0.0,dot_product(s*-1,m));



                    double temp = 2 * dot_product(s,m);
                    Point reflection_dir = s - normal_direction * temp;
                    reflection_dir = normalized_point(reflection_dir);

                    double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));


                    color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                    color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                    color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                    color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];

                }
            }
            for(int i=0;i<spotLights.size();i++)
            {
                Ray spot_point_light_ray;
                spot_point_light_ray.start = spotLights[i]->point_light.light_pos;
                spot_point_light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);

                Point spot_light_dir = spotLights[i] -> light_direction;
                spot_light_dir = normalized_point(spot_light_dir);

                double cos_theta = dot_product(spot_light_dir,spot_point_light_ray.dir);
                double theta = acos(cos_theta) * 180.0 / acos(-1);

                if(theta < spotLights[i]->cutoff_angle)
                {

                    double dist_t = getLength(spotLights[i]->point_light.light_pos - intersection_point);

                    if(dist_t < ephsilon)
                    {
                        continue;
                    }

                    Ray light_ray;
                    light_ray.start = spotLights[i]->point_light.light_pos;
                    light_ray.dir = normalized_point(intersection_point - spotLights[i]->point_light.light_pos);


                    ///calc normal
                    //Point normal_direction = intersection_point - reference_point;
                    normal_direction = normalized_point(normal_direction);

                    bool shadow = false;
                    for(int i=0;i<objects.size();i++)
                    {
                        double color[3];
                        double temp_t = objects[i]->intersect(light_ray,color,0);

                        if(temp_t != -1 && temp_t + ephsilon < dist_t)
                        {
                            shadow = true;
                            break;
                        }
                    }
                    if(shadow == false)
                    {
                        Point s = light_ray.dir;
                        Point v = ray.dir;
                        Point m = normal_direction;

                        s = normalized_point(s);
                        m = normalized_point(m);

                        double lambart_val = max(0.0,dot_product(s*-1,m));


                        double temp = 2 * dot_product(s,m);
                        Point reflection_dir = s - normal_direction * temp;
                        reflection_dir = normalized_point(reflection_dir);

                        double phongValue = max(0.0, dot_product(reflection_dir,ray.dir*-1));



                        color[0] += pointLights[i]->color[0] * coEfficients[1] * lambart_val * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[1] * lambart_val * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[1] * lambart_val * intersection_color[2];

                        color[0] += pointLights[i]->color[0] * coEfficients[2] * pow(phongValue,shine) * intersection_color[0];
                        color[1] += pointLights[i]->color[1] * coEfficients[2] * pow(phongValue,shine) * intersection_color[1];
                        color[2] += pointLights[i]->color[2] * coEfficients[2] * pow(phongValue,shine) * intersection_color[2];

                        //printf("after: %lf %lf %lf\n",intersection_color[0],intersection_color[1],intersection_color[2]);

                    }
                }
            }
            if(level >= recursion_level)
                return t;

            //construct reflected ray
            Point s = ray.dir;

            Point m = normal_direction;

            s = normalized_point(s);
            m = normalized_point(m);

            double temp = 2 * dot_product(s,m);
            Point reflection_dir = s - normal_direction * temp;
            reflection_dir = normalized_point(reflection_dir);

            // actually slightly forward from the point (by moving the
            // start a little bit towards the reflection direction)
            // to avoid self intersection

            Ray desired_reflected_ray;
            desired_reflected_ray.start = intersection_point + reflection_dir * ephsilon;
            desired_reflected_ray.dir = reflection_dir;

            double tmpcolor[3];
			double mini_t = INT_MAX;
			int index = -1;
			for(int k=0;k<objects.size();k++)
			{
				double t = objects[k]->intersect(desired_reflected_ray,tmpcolor, 0);

                if(t > 0 && t < mini_t)
                {
                    mini_t = t;
                    index = k;
                }

			}
			// if nearest object is found, then shade the pixel
			if(index != -1)
			{

				double getcolor[3];

				double t = objects[index]->intersect(desired_reflected_ray,getcolor,level+1);

				color[0] += getcolor[0] * coEfficients[3];
                color[1] += getcolor[1] * coEfficients[3];
                color[2] += getcolor[2] * coEfficients[3];

               // printf("%lf %lf %lf\n",color[0],color[1],color[2]);
			}

            return t;

        }
        else
        {
            return -1.0;
        }

    }

};
