#ifndef MESH_H
#define MESH_H

#include <QGLWidget>
#include <iterator>

// TO MODIFY
// Advice => You should create a new file for each module but not necessarily for each class


class vertex_it;
class face_around_it;

class Point
{
    double _x;
    double _y;
    double _z;
    int face_index;
    double courbure;
    double Laplacien[3]={0.0,0.0,0.0};

public:
    Point():_x(),_y(),_z() {}
    Point(float x_, float y_, float z_):_x(x_),_y(y_),_z(z_),face_index(0) {}
    // get
    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }
    int get_face_index() const { return face_index; }
    void set_face_index(int index) { face_index=index; }
    void set_courbure(double value_courbure) { courbure=value_courbure; }
    double get_courbure() { return courbure; }
    double* get_laplacien() { return Laplacien; }
    void set_laplacien(double x, double y, double z) { Laplacien[0]=x;Laplacien[1]=y;Laplacien[2]=z; }

};

//** TO MODIFY
class Mesh
{
public:
    QVector<Point> vertexTab;
    QVector<QVector<int>> facesTab;
    QVector<QVector<int>> adjFacesTab;

public:
    Mesh();
    
    void drawMesh();
    void drawMeshWireFrame();
    void drawMeshCurvature();

    vertex_it begin();
    vertex_it end();

    face_around_it begin(int vertex_index);
    face_around_it end(int vertex_index);

    void compute_normal_and_mean();
    QVector<int> split_triangle(const int face_index,const int vertex_index);
    void flip_edge(int face_1_index,int face_2_index);

    int nb_faces;
};



class vertex_it
{
public:
    int m_index=0;
    Mesh* m_mesh_ptr=nullptr;

public:
    vertex_it(Mesh* mesh_ptr,int index);
    Point& operator*() const {return (m_mesh_ptr->vertexTab[m_index]);}
    void operator++() {m_index+=1;}
    bool operator==(const vertex_it& other) const { return m_index==other.m_index;}
    bool operator!=(const vertex_it& other) const { return !(*this == other);}
};



class face_around_it
{
    int m_vertex_index=0;
    int m_face_index=0;
    Mesh* m_mesh_ptr=nullptr;

public:
    face_around_it(Mesh* mesh_ptr,int vertex_index,int face_index);
    QVector<int>& operator*() const {return (m_mesh_ptr->facesTab[m_face_index]);}
    void operator++();
    bool operator==(const face_around_it& other) const { return (m_vertex_index==other.m_vertex_index)&&(m_face_index==other.m_face_index);}
    bool operator!=(const face_around_it& other) const { return !(*this == other);}
};


class Triangulation :public Mesh
{
public:
    int P_inf_index;
    double x_moy;
    double y_moy;
    double z_moy;

    //Constructeur
    Triangulation();
    QVector<int> add_point_in_triangulation(int pt_index);
    void triangulate();
    //Réécriture du drawmesh pour ne pas prendre en compte les triangles avec pt infini
    void drawMeshWireFrame();

};

// Fonctions annexes non methodes de classe
double segment_intersection_test(Point& A, Point& P, Point* C, Point* D);
double orientation_test(Point& pt1,Point& pt2,Point& pt3);
double in_triangle_test(QVector<Point*> triangle, Point& pt);
double in_circumscribed_circle_test(QVector<Point*> triangle, Point& pt);



#endif // MESH_H
