#include "mesh.h"
#include <QTextStream>
#include <QFile>
#include <iostream>
#include <map>
#include <cmath>
#include <queue>

Mesh::Mesh()
{
    //vertexTab.push_back(Point(-0.5,-0.5,-0.5)); //0
    //vertexTab.push_back(Point(0.5,-0.5,-0.5)); // 1
    //vertexTab.push_back(Point(0,0.5,-0.5)); // 2
    //vertexTab.push_back(Point(0,-0.5,0.5)); // 3

    //Lecture du fichier OFF
    QString fileName="C:/empty.off";
   // QString fileName="C:/queen.off";
    QFile fichier(fileName);
    fichier.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream flux(&fichier);

    QString ligne;
    //Premiere ligne
    ligne=flux.readLine();
    QStringList mots = ligne.split(" ");

    int nb_vertices=mots.at(0).toInt();
    nb_faces=mots.at(1).toInt();

    int indice_face =0;
    //Création de la map pour l'adjacence
    std::map<std::pair<int,int>,std::pair<int,int>> myMap;

    while(! flux.atEnd())
    {
        ligne=flux.readLine();
        QStringList mots = ligne.split(" ");
        //Les sommets
        if(mots.length() == 3)
        {
            //Lecture du sommet
            double x = mots.at(0).toDouble();
            double y = mots.at(1).toDouble();
            double z = mots.at(2).toDouble();
            //Ajout du sommet à la liste
            vertexTab.push_back(Point(x,y,z));
        }
        //Les faces
        if(mots.length() == 4)
        {
            //Lecture de la face
            int k_simplex = mots.at(0).toInt();
            int indice_0 = mots.at(1).toInt();
            int indice_1 = mots.at(2).toInt();
            int indice_2 = mots.at(3).toInt();


            adjFacesTab.push_back({0,0,0});
            //Premiere paire : 1-2
            if(myMap.find(std::make_pair(indice_2,indice_1)) != myMap.end())
            {
                //On a trouvé la face adjacente
                auto value=myMap.find(std::make_pair(indice_2,indice_1));
                int Face_adjacente=value->second.first;
                int place_dans_face_adjacente =value->second.second;
                //On modifie la face actuelle et adjacente
                adjFacesTab[indice_face][0]=Face_adjacente;
                adjFacesTab[Face_adjacente][place_dans_face_adjacente]=indice_face;
            }
            else
            {
                //On n'a pas trouvé la face adjacente
                //On ajoute notre face à la map
                myMap.insert(std::make_pair(std::make_pair(indice_1,indice_2),std::make_pair(indice_face,0)));
            }

            //Deuxieme paire : 2-0
            if(myMap.find(std::make_pair(indice_0,indice_2)) != myMap.end())
            {
                //On a trouvé la face adjacente
                auto value=myMap.find(std::make_pair(indice_0,indice_2));
                int Face_adjacente=value->second.first;
                int place_dans_face_adjacente =value->second.second;
                //On modifie la face actuelle et adjacente
                adjFacesTab[indice_face][1]=Face_adjacente;
                adjFacesTab[Face_adjacente][place_dans_face_adjacente]=indice_face;
            }
            else
            {
                //On n'a pas trouvé la face adjacente
                //On ajoute notre face à la map
                myMap.insert(std::make_pair(std::make_pair(indice_2,indice_0),std::make_pair(indice_face,1)));
            }


            //Troisieme paire : 0-1
            if(myMap.find(std::make_pair(indice_1,indice_0)) != myMap.end())
            {
                //On a trouvé la face adjacente
                auto value=myMap.find(std::make_pair(indice_1,indice_0));
                int Face_adjacente=value->second.first;
                int place_dans_face_adjacente =value->second.second;
                //On modifie la face actuelle et adjacente
                adjFacesTab[indice_face][2]=Face_adjacente;
                adjFacesTab[Face_adjacente][place_dans_face_adjacente]=indice_face;
            }
            else
            {
                //On n'a pas trouvé la face adjacente
                //On ajoute notre face à la map
                myMap.insert(std::make_pair(std::make_pair(indice_0,indice_1),std::make_pair(indice_face,2)));
            }

            //Ajout de la face à la liste
            QVector<int> indices_sommets {indice_0,indice_1,indice_2};
            facesTab.push_back(indices_sommets);
            vertexTab[indice_0].set_face_index(indice_face);
            vertexTab[indice_1].set_face_index(indice_face);
            vertexTab[indice_2].set_face_index(indice_face);

            indice_face+=1;
        }
    }

    compute_normal_and_mean();

}

// The following functions could be displaced into a module OpenGLDisplayMesh that would include Mesh
// Draw a vertex
void glVertexDraw(const Point & p) {
    glVertex3f(p.x(), p.y(), p.z());
}

void Mesh::drawMesh() {
    for(int i = 0; i < nb_faces; i+=1) {

        glColor3d(1,1,0);

        //glBegin(GL_TRIANGLES);
        //glVertexDraw(vertexTab[faces[i]]);
        //glVertexDraw(vertexTab[faces[i+1]]);
        //glVertexDraw(vertexTab[faces[i+2]]);
        //glEnd();

        /* Visualisation faces adjacentes
        int j=148;
        if(i==j || i==adjFacesTab[j][0] || i==adjFacesTab[j][2] || i==adjFacesTab[j][2])
        {
            glColor3d(1,0,0);
        }
        else
        {
            glColor3d(0,1,0);
        }
        */

        glBegin(GL_TRIANGLES);
        glVertexDraw(vertexTab[facesTab.at(i).at(0)]);
        glVertexDraw(vertexTab[facesTab.at(i).at(1)]);
        glVertexDraw(vertexTab[facesTab.at(i).at(2)]);
        glEnd();
    }
}

//Example with a wireframe tedraedra
void Mesh::drawMeshWireFrame() {
    for(int i = 0; i < nb_faces; i+=1) {
        glBegin(GL_LINE_STRIP);
            glVertexDraw(vertexTab[facesTab.at(i).at(0)]);
            glVertexDraw(vertexTab[facesTab.at(i).at(1)]);
        glEnd();
        glBegin(GL_LINE_STRIP);
            glVertexDraw(vertexTab[facesTab.at(i).at(1)]);
            glVertexDraw(vertexTab[facesTab.at(i).at(2)]);
        glEnd();
        glBegin(GL_LINE_STRIP);
            glVertexDraw(vertexTab[facesTab.at(i).at(2)]);
            glVertexDraw(vertexTab[facesTab.at(i).at(0)]);
        glEnd();
    }

}



//Example with a wireframe tedraedra
void Mesh::drawMeshCurvature() {

    //Calcule courbure max:
    double max_curv=-10000.0;
    double min_curv=+10000.0;
    for(int i=0; i<vertexTab.size();i++)
    {
        if(vertexTab[i].get_courbure()>max_curv)
        {
            max_curv=vertexTab[i].get_courbure();
        }
        if(vertexTab[i].get_courbure()<max_curv)
        {
            min_curv=vertexTab[i].get_courbure();
        }
    }


    for(int i = 0; i < nb_faces; i+=1) {

        double mean_curv_face=0.0;
        mean_curv_face+=vertexTab[facesTab[i][0]].get_courbure();
        mean_curv_face+=vertexTab[facesTab[i][1]].get_courbure();
        mean_curv_face+=vertexTab[facesTab[i][2]].get_courbure();
        mean_curv_face=mean_curv_face/3.0;

        double normalized_curv=(log(mean_curv_face)-log(min_curv))/(log(max_curv)-log(min_curv));
        if(normalized_curv*10>0.5) {glColor3d(normalized_curv*10,0.0,0.0);}
        else if(normalized_curv*100>0) {glColor3d(0.0,normalized_curv*100,0.0);}
        else {glColor3d(0.0,0.0,normalized_curv*1000000);}
        /*
        if(normalized_curv<0.33){glColor3d(1,0,0);}
        else if(normalized_curv<0.33){glColor3d(0,1,0);}
        else {glColor3d(0,0,1);}
        */

        glBegin(GL_TRIANGLES);
        glVertexDraw(vertexTab[facesTab.at(i).at(0)]);
        glVertexDraw(vertexTab[facesTab.at(i).at(1)]);
        glVertexDraw(vertexTab[facesTab.at(i).at(2)]);
        glEnd();
    }

}



void Mesh::compute_normal_and_mean()
{
    for(vertex_it vertex_iterator=this->begin();vertex_iterator!=this->end();++vertex_iterator)
    {
        float A=0.0f;
        float laplacien[3] = {0.0f,0.0f,0.0f};
        int vertex_index= vertex_iterator.m_index;

        for(face_around_it face_around_iterator=this->begin(vertex_iterator.m_index);face_around_iterator!=this->end(vertex_iterator.m_index);++face_around_iterator)
        {
             QVector<int>& this_triangle = (*face_around_iterator);
             face_around_it iterator_copy(face_around_iterator);
             ++iterator_copy;
             QVector<int>& next_triangle = (*iterator_copy);






             Point i(0.0f,0.0f,0.0f);
             Point j(0.0f,0.0f,0.0f);
             Point k(0.0f,0.0f,0.0f);

             Point Pt_i(0.0f,0.0f,0.0f);
             Point Pt_j(0.0f,0.0f,0.0f);

             if(this_triangle[0]==vertex_index)
             {
                 i=vertexTab[this_triangle[0]];
                 j=vertexTab[this_triangle[2]];
                 k=vertexTab[this_triangle[1]];

                 Pt_i=vertexTab[this_triangle[0]];
                 Pt_j=vertexTab[this_triangle[2]];
             }
             else if(this_triangle[1]==vertex_index)
             {
                 i=vertexTab[this_triangle[1]];
                 j=vertexTab[this_triangle[0]];
                 k=vertexTab[this_triangle[2]];

                 Pt_i=vertexTab[this_triangle[1]];
                 Pt_j=vertexTab[this_triangle[0]];
             }
             else
             {
                 i=vertexTab[this_triangle[2]];
                 j=vertexTab[this_triangle[1]];
                 k=vertexTab[this_triangle[0]];

                 Pt_i=vertexTab[this_triangle[2]];
                 Pt_j=vertexTab[this_triangle[1]];
             }

             Point vec_ki(i.x()-k.x(),i.y()-k.y(),i.z()-k.z());
             Point vec_kj(j.x()-k.x(),j.y()-k.y(),j.z()-k.z());

             float norme_ki=sqrt(vec_ki.x()*vec_ki.x()+vec_ki.y()*vec_ki.y()+vec_ki.z()*vec_ki.z());
             float norme_kj=sqrt(vec_kj.x()*vec_kj.x()+vec_kj.y()*vec_kj.y()+vec_kj.z()*vec_kj.z());

             Point vec_ki_normalise(vec_ki.x()/norme_ki,vec_ki.y()/norme_ki,vec_ki.z()/norme_ki);
             Point vec_kj_normalise(vec_kj.x()/norme_kj,vec_kj.y()/norme_kj,vec_kj.z()/norme_kj);

             float cos_alpha_ij=0.0f;
             cos_alpha_ij+=vec_ki_normalise.x()*vec_kj_normalise.x();
             cos_alpha_ij+=vec_ki_normalise.y()*vec_kj_normalise.y();
             cos_alpha_ij+=vec_ki_normalise.z()*vec_kj_normalise.z();

             float sin_alpha_ij=sqrt(1-cos_alpha_ij*cos_alpha_ij);
             float cotan_alpha_ij=cos_alpha_ij/sin_alpha_ij;







             if(next_triangle[0]==vertex_index)
             {
                 i=vertexTab[next_triangle[0]];
                 j=vertexTab[next_triangle[1]];
                 k=vertexTab[next_triangle[2]];
             }
             else if(next_triangle[1]==vertex_index)
             {
                 i=vertexTab[next_triangle[1]];
                 j=vertexTab[next_triangle[2]];
                 k=vertexTab[next_triangle[0]];
             }
             else
             {
                 i=vertexTab[next_triangle[2]];
                 j=vertexTab[next_triangle[0]];
                 k=vertexTab[next_triangle[1]];
             }

             Point new_vec_ki(i.x()-k.x(),i.y()-k.y(),i.z()-k.z());
             Point new_vec_kj(j.x()-k.x(),j.y()-k.y(),j.z()-k.z());

             float new_norme_ki=sqrt(new_vec_ki.x()*new_vec_ki.x()+new_vec_ki.y()*new_vec_ki.y()+new_vec_ki.z()*new_vec_ki.z());
             float new_norme_kj=sqrt(new_vec_kj.x()*new_vec_kj.x()+new_vec_kj.y()*new_vec_kj.y()+new_vec_kj.z()*new_vec_kj.z());

             Point new_vec_ki_normalise(new_vec_ki.x()/new_norme_ki,new_vec_ki.y()/new_norme_ki,new_vec_ki.z()/new_norme_ki);
             Point new_vec_kj_normalise(new_vec_kj.x()/new_norme_kj,new_vec_kj.y()/new_norme_kj,new_vec_kj.z()/new_norme_kj);

             float cos_beta_ij=0.0f;
             cos_beta_ij+=new_vec_ki_normalise.x()*new_vec_kj_normalise.x();
             cos_beta_ij+=new_vec_ki_normalise.y()*new_vec_kj_normalise.y();
             cos_beta_ij+=new_vec_ki_normalise.z()*new_vec_kj_normalise.z();

             float sin_beta_ij=sqrt(1-cos_beta_ij*cos_beta_ij);
             float cotan_beta_ij=cos_beta_ij/sin_beta_ij;


             laplacien[0]+=(cotan_alpha_ij+cotan_beta_ij)*(Pt_j.x()-Pt_i.x());
             laplacien[1]+=(cotan_alpha_ij+cotan_beta_ij)*(Pt_j.y()-Pt_i.y());
             laplacien[2]+=(cotan_alpha_ij+cotan_beta_ij)*(Pt_j.z()-Pt_i.z());







             float prod_vect_x=vec_ki_normalise.y()*vec_kj_normalise.z()-vec_ki_normalise.z()*vec_kj_normalise.y();
             float prod_vect_y=vec_ki_normalise.z()*vec_kj_normalise.x()-vec_ki_normalise.x()*vec_kj_normalise.z();
             float prod_vect_z=vec_ki_normalise.x()*vec_kj_normalise.y()-vec_ki_normalise.y()*vec_kj_normalise.x();

             A+=0.5f*(prod_vect_x*prod_vect_x+prod_vect_y*prod_vect_y+prod_vect_z*prod_vect_z);
        }
        A=A/6;
        Point Laplacien(laplacien[0]/A,laplacien[1]/A,laplacien[2]/A);

        face_around_it face_around_iterator=this->begin(vertex_iterator.m_index);
        QVector<int> one_triangle=(*face_around_iterator);
        Point i=vertexTab[one_triangle[0]];
        Point j=vertexTab[one_triangle[2]];
        Point k=vertexTab[one_triangle[1]];
        Point vec_ki(i.x()-k.x(),i.y()-k.y(),i.z()-k.z());
        Point vec_kj(j.x()-k.x(),j.y()-k.y(),j.z()-k.z());
        double norme_ki=sqrt(vec_ki.x()*vec_ki.x()+vec_ki.y()*vec_ki.y()+vec_ki.z()*vec_ki.z());
        double norme_kj=sqrt(vec_kj.x()*vec_kj.x()+vec_kj.y()*vec_kj.y()+vec_kj.z()*vec_kj.z());
        Point vec_ki_normalise(vec_ki.x()/norme_ki,vec_ki.y()/norme_ki,vec_ki.z()/norme_ki);
        Point vec_kj_normalise(vec_kj.x()/norme_kj,vec_kj.y()/norme_kj,vec_kj.z()/norme_kj);
        double normale_x=vec_ki_normalise.y()*vec_kj_normalise.z()-vec_ki_normalise.z()*vec_kj_normalise.y();
        double normale_y=vec_ki_normalise.z()*vec_kj_normalise.x()-vec_ki_normalise.x()*vec_kj_normalise.z();
        double normale_z=vec_ki_normalise.x()*vec_kj_normalise.y()-vec_ki_normalise.y()*vec_kj_normalise.x();


        double courbure= Laplacien.x()*Laplacien.x()+Laplacien.y()*Laplacien.y()+Laplacien.z()*Laplacien.z();
        if(normale_x*Laplacien.x()+normale_y*Laplacien.y()+normale_z*Laplacien.z()>0)
        {
            (*vertex_iterator).set_courbure(courbure);
        }
        else
        {
            (*vertex_iterator).set_courbure(-1.0*courbure);
        }
        (*vertex_iterator).set_laplacien(Laplacien.x(),Laplacien.y(),Laplacien.z());
    }
}






void Mesh::split_triangle(const int face_index, const int vertex_index)
{
    // LES FACES
    int A_index(facesTab[face_index][0]);
    int B_index(facesTab[face_index][1]);
    int C_index(facesTab[face_index][2]);
    //Modification de la face principale:
    facesTab[face_index][2]=vertex_index;
    //Creation face f1:
    QVector<int> f1 = {B_index,C_index,vertex_index};
    int f1_index=facesTab.size();
    facesTab.push_back(f1);
    //Creation face f2:
    QVector<int> f2 = {C_index,A_index,vertex_index};
    int f2_index=facesTab.size();
    facesTab.push_back(f2);

    // LES FACES ADJACENTES
    //Pour face principale:
    int fA_index(adjFacesTab[face_index][0]);
    int fB_index(adjFacesTab[face_index][1]);
    int fC_index(adjFacesTab[face_index][2]);

    adjFacesTab[face_index][0]=f1_index;
    adjFacesTab[face_index][1]=f2_index;
    //Pour f1:
    QVector<int> adj_f1 = {f2_index,face_index,fA_index};
    adjFacesTab.push_back(adj_f1);
    //Pour f2:
    QVector<int> adj_f2 = {face_index,f1_index,fB_index};
    adjFacesTab.push_back(adj_f2);
    //Pour fA:
    for (int i=0;i<3;i++)
    {
        if(adjFacesTab[fA_index][i]==face_index)
        {
            adjFacesTab[fA_index][i]=f1_index;
            break;
        }
    }
    //Pour fB:
    for (int i=0;i<3;i++)
    {
        if(adjFacesTab[fB_index][i]==face_index)
        {
            adjFacesTab[fB_index][i]=f2_index;
            break;
        }
    }
    // GESTION DES SOMMETS
    vertexTab[A_index].set_face_index(face_index);
    vertexTab[B_index].set_face_index(face_index);
    vertexTab[C_index].set_face_index(fA_index);
    vertexTab[vertex_index].set_face_index(face_index);

    //std::cout << "f:"<< face_index << " f1:" << f1_index << " f2:" << f2_index << std::endl;
    //std::cout << "voisins de f --- " << adjFacesTab[face_index][0] << " " <<  adjFacesTab[face_index][1] << " " <<  adjFacesTab[face_index][2] << std::endl;
    //std::cout << "voisins de f1 --- " << adjFacesTab[f1_index][0] << " " <<  adjFacesTab[f1_index][1] << " " <<  adjFacesTab[f1_index][2] << std::endl;
    //std::cout << "voisins de f2 --- " << adjFacesTab[f2_index][0] << " " <<  adjFacesTab[f2_index][1] << " " <<  adjFacesTab[f2_index][2] << std::endl;

}





void Mesh::flip_edge(int face_1_index,int face_2_index)
{
    QVector<int> triangle_f = facesTab[face_1_index];
    QVector<int> triangle_g = facesTab[face_2_index];
    //
    //Rehcherche de l'arete
    int arete_1;
    int arete_2;

    if(triangle_g[0]==triangle_f[0] || triangle_g[1]==triangle_f[0] || triangle_g[2]==triangle_f[0])
    {
        if(triangle_g[0]==triangle_f[1] || triangle_g[1]==triangle_f[1] || triangle_g[2]==triangle_f[1])
        {
            arete_1=triangle_f[0];
            arete_2=triangle_f[1];
        }
        else
        {
            arete_1=triangle_f[2];
            arete_2=triangle_f[0];
        }
    }
    else
    {
        arete_1=triangle_f[1];
        arete_2=triangle_f[2];
    }

    //
    /*
    //Recherche des sommets pour f
    if((triangle_f[0]==arete_1 && triangle_f[1]==arete_2)||(triangle_f[1]==arete_1 && triangle_f[2]==arete_2)||(triangle_f[2]==arete_1 && triangle_f[0]==arete_2))
    {
        //Alors bien orienté
    }
    else
    {
        int tampon = arete_1;
        arete_1=arete_2;
        arete_2=tampon;
    }
    */
    int Af= arete_1;
    int Bf= arete_2;
    int Cf;
    if(triangle_f[0]!=arete_1 && triangle_f[0]!=arete_2) {Cf=triangle_f[0];}
    else if(triangle_f[1]!=arete_1 && triangle_f[1]!=arete_2) {Cf=triangle_f[1];}
    else {Cf=triangle_f[2];}

    //Recherche des sommets pour g
    /*
    if((triangle_g[0]==arete_2 && triangle_g[1]==arete_1)||(triangle_g[1]==arete_2 && triangle_g[2]==arete_1)||(triangle_g[2]==arete_2 && triangle_g[0]==arete_1))
    {
        //Alors bien orienté
    }
    else
    {
        int tampon = arete_1;
        arete_1=arete_2;
        arete_2=tampon;
    }
    */
    int Ag= arete_2;
    int Bg= arete_1;
    int Cg;
    if(triangle_g[0]!=arete_1 && triangle_g[0]!=arete_2) {Cg=triangle_g[0];}
    else if(triangle_g[1]!=arete_1 && triangle_g[1]!=arete_2) {Cg=triangle_g[1];}
    else {Cg=triangle_g[2];}

    //
    // Recherche des indices des faces adjacentes pour f
    int fA;
    int fB;
    int fC=face_2_index;
    /*
    std::cout << " ---- " << std::endl;
    std::cout << face_1_index << " -- ";
    std::cout <<facesTab[face_1_index][0] << " ";
    std::cout <<facesTab[face_1_index][1] << " ";
    std::cout <<facesTab[face_1_index][2] << " " << std::endl;
    */
    for(int i=0;i<3;i++)
    {
        /*
        std::cout << adjFacesTab[face_1_index][i] << " - ";
        std::cout <<facesTab[adjFacesTab[face_1_index][i]][0] << " ";
        std::cout <<facesTab[adjFacesTab[face_1_index][i]][1] << " ";
        std::cout <<facesTab[adjFacesTab[face_1_index][i]][2] << " " << std::endl;
        */

        if(adjFacesTab[face_1_index][i]!=face_2_index)
        {
            QVector<int> face = facesTab[adjFacesTab[face_1_index][i]];
            if(face[0]!=Af && face[1]!=Af && face[2]!=Af)
            {
                fA=adjFacesTab[face_1_index][i];
            }
            else
            {
                fB=adjFacesTab[face_1_index][i];
            }
        }
    }

    // Recherche des indices des faces adjacentes pour g
    int gA;
    int gB;
    int gC=face_1_index;
    for(int i=0;i<3;i++)
    {
        if(adjFacesTab[face_2_index][i]!=face_1_index)
        {
            QVector<int> face = facesTab[adjFacesTab[face_2_index][i]];
            if(face[0]!=Ag && face[1]!=Ag && face[2]!=Ag)
            {
                gA=adjFacesTab[face_2_index][i];
            }
            else
            {
                gB=adjFacesTab[face_2_index][i];
            }
        }
    }

    //
    // Modification des faces
    QVector<int> new_f={Cg,Bf,Cf};
    QVector<int> new_g={Cf,Bg,Cg};
    facesTab[face_1_index]=new_f;
    facesTab[face_2_index]=new_g;
    //Modification des adjacences
    QVector<int> new_f_adj={fA,face_2_index,gB};
    QVector<int> new_g_adj={gA,face_1_index,fB};
    adjFacesTab[face_1_index]=new_f_adj;
    adjFacesTab[face_2_index]=new_g_adj;

    for(int i=0;i<3;i++)
    {
        if(adjFacesTab[fB][i]==face_1_index)
        {
            adjFacesTab[fB][i]=face_2_index;
        }
    }
    for(int i=0;i<3;i++)
    {
        if(adjFacesTab[gB][i]==face_2_index)
        {
            adjFacesTab[gB][i]=face_1_index;
        }
    }

    //
    //Modification des sommets
    vertexTab[Cg].set_face_index(face_1_index);
    vertexTab[Bf].set_face_index(face_1_index);
    vertexTab[Cf].set_face_index(face_1_index);
    vertexTab[Cf].set_face_index(face_2_index);
    vertexTab[Bg].set_face_index(face_2_index);
    vertexTab[Cg].set_face_index(face_2_index);
}




vertex_it Mesh::begin()
{
    return vertex_it(this,0);
}
vertex_it Mesh::end()
{
    return vertex_it(this,vertexTab.size());
}

face_around_it Mesh::begin(int vertex_index)
{
    face_around_it last(this,vertex_index,vertexTab[vertex_index].get_face_index());
    ++last;
    return last;
}

face_around_it Mesh::end(int vertex_index)
{
    return face_around_it(this,vertex_index,vertexTab[vertex_index].get_face_index());
}




/////////////////////////// ITERATEURS /////////////////////

vertex_it::vertex_it(Mesh* mesh_ptr,int index)
{
    m_index=index;
    m_mesh_ptr=mesh_ptr;
}

face_around_it::face_around_it(Mesh* mesh_ptr,int vertex_index,int face_index)
{
    m_vertex_index=vertex_index;
    m_face_index=face_index;
    m_mesh_ptr=mesh_ptr;
}

void face_around_it::operator++()
{
    if(m_mesh_ptr->facesTab[m_face_index][0]==m_vertex_index)
    {
        m_face_index=m_mesh_ptr->adjFacesTab[m_face_index][1];
    }
    else if(m_mesh_ptr->facesTab[m_face_index][1]==m_vertex_index)
    {
        m_face_index=m_mesh_ptr->adjFacesTab[m_face_index][2];
    }
    else
    {
        m_face_index=m_mesh_ptr->adjFacesTab[m_face_index][0];
    }
}





///////////////////////////////////////////////////////////////
/////////////////////// TRIANGULATION /////////////////////////
/// ///////////////////////////////////////////////////////////

// CLASSE TRAINGULATION

//Constructeur
Triangulation::Triangulation()
{
    //reinitialisation après constructeur de mesh
    vertexTab={};
    facesTab={};
    adjFacesTab={};
    nb_faces=0;

    //Lecture du fichier OFF
    //QString fileName="C:/test.off";
    QString fileName="C:/franke5.off";
    QFile fichier(fileName);
    fichier.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream flux(&fichier);

    QString ligne;
    //Premiere ligne
    ligne=flux.readLine();
    //Deuxieme ligne
    ligne=flux.readLine();
    QStringList mots = ligne.split(" ");

    int nb_vertices=mots.at(0).toInt();
    nb_faces=mots.at(1).toInt();

    int indice_face =0;

    while(! flux.atEnd())
    {
        ligne=flux.readLine();
        QStringList mots = ligne.split(" ");
        //Les sommets
        if(mots.length() == 3)
        {
            //Lecture du sommet
            double x = mots.at(0).toDouble();
            double y = mots.at(1).toDouble();
            double z = mots.at(2).toDouble();
            //Ajout du sommet à la liste
            vertexTab.push_back(Point(x,y,z));
        }
    }

    Point A(0,0,0);
    Point B(0.5,1,0);
    Point C(1,0,0);
    Point D(0.5,0.5,0);

    triangulate();
}


// pt est  un point 2D (x,y,z=0)
void Triangulation::add_point_in_triangulation(int pt_index)
{
    //Chercher les triangles qui contiennent le pt (mais pas parmi les 3 triangles infinis)
    int nb_faces=facesTab.size();
    for(int i=3;i<nb_faces;i++)
    {
        QVector<Point*> triangle;
        triangle.push_back(&vertexTab[facesTab[i][0]]);
        triangle.push_back(&vertexTab[facesTab[i][1]]);
        triangle.push_back(&vertexTab[facesTab[i][2]]);
        double test=in_triangle_test(triangle,vertexTab[pt_index]);
        if(test>0.0)
        {
            //Alors split le triangle
            split_triangle(i,pt_index);
            break;
        }
        else if(test==0.0)
        {
            std::cout << "Vertex " << pt_index << "sur bord d'un triangle"<< std::endl;
        }
    }
}


//Reecriture de la fonction drawmesh specialement pour la triangulation
//Les triangles avec sommet infini n'apparaissent pas
void Triangulation::drawMeshWireFrame()
{
    nb_faces=facesTab.size();
    for(int i = 0; i < nb_faces; i+=1) {
        //Cas d'un triangle avec pt infini
        if(facesTab.at(i).at(0)==P_inf_index ||facesTab.at(i).at(1)==P_inf_index ||facesTab.at(i).at(2)==P_inf_index)
        {
            continue;
        }
        //Affichage classique
        glBegin(GL_LINE_STRIP);
            glVertexDraw(vertexTab[facesTab.at(i).at(0)]);
            glVertexDraw(vertexTab[facesTab.at(i).at(1)]);
        glEnd();
        glBegin(GL_LINE_STRIP);
            glVertexDraw(vertexTab[facesTab.at(i).at(1)]);
            glVertexDraw(vertexTab[facesTab.at(i).at(2)]);
        glEnd();
        glBegin(GL_LINE_STRIP);
            glVertexDraw(vertexTab[facesTab.at(i).at(2)]);
            glVertexDraw(vertexTab[facesTab.at(i).at(0)]);
        glEnd();
    }

}



void Triangulation::triangulate()
{
    //Calcul des coords min et max en x et en y
    double x_min=vertexTab[0].x();
    double x_max=vertexTab[0].x();
    double y_min=vertexTab[0].y();
    double y_max=vertexTab[0].y();
    double z_min=vertexTab[0].z();
    double z_max=vertexTab[0].z();

    for(int i=0;i<vertexTab.size();i++)
    {
        double x=vertexTab[i].x();
        double y=vertexTab[i].y();
        double z=vertexTab[i].z();
        if(x<x_min){x_min=x;}
        if(x>x_max){x_max=x;}
        if(y<y_min){y_min=y;}
        if(y>y_max){y_max=y;}
        if(z<z_min){z_min=z;}
        if(z>z_max){z_max=z;}
    }


    //Creation du triangle fictif et du point infini
    int P1_index=vertexTab.size();
    int P2_index=P1_index+1;
    int P3_index=P2_index+1;
    int Pinfini_index=P3_index+1;

    P_inf_index=Pinfini_index;

    double epsilon=0.1;
    vertexTab.push_back(Point(0.5*(x_min+x_max),y_max+(y_max-y_min),0.0)); //P1
    vertexTab.push_back(Point(x_min-(x_max-x_min)-epsilon,y_min-epsilon,0.0)); //P2
    vertexTab.push_back(Point(x_max+(x_max-x_min)+epsilon,y_min-epsilon,0.0)); //P3
    vertexTab.push_back(Point(0.5*(x_min+x_max),0.5*(y_min+y_max),z_min-1.0)); //Pinfini


    int inf_face_1=facesTab.size();
    int inf_face_2=inf_face_1+1;
    int inf_face_3=inf_face_2+1;
    int fake_face_index=inf_face_3+1;
    facesTab.push_back({P1_index,Pinfini_index,P2_index});//F_infini_1
    facesTab.push_back({P2_index,Pinfini_index,P3_index});//F_infini_2
    facesTab.push_back({P3_index,Pinfini_index,P1_index});//F_infini_3
    facesTab.push_back({P1_index,P2_index,P3_index});//Fake face
    adjFacesTab.push_back({inf_face_2,fake_face_index,inf_face_3});//F_infini_1
    adjFacesTab.push_back({inf_face_3,fake_face_index,inf_face_1});//F_infini_2
    adjFacesTab.push_back({inf_face_1,fake_face_index,inf_face_2});//F_infini_3
    adjFacesTab.push_back({inf_face_2,inf_face_3,inf_face_1});//Fake face

    //Parcourt les sommets et ajout dans la triangulation
    int nb_true_vertices=vertexTab.size()-4;
    for(int i=0;i<nb_true_vertices;i++)
    {
        add_point_in_triangulation(i);
        QVector<QVector<int>> pile;
        bool again=true;
        while(again==true)
        {
            //Parcourt des faces (sauf les faces infinies)
            int nombre_faces=facesTab.size();
            for(int f=3;f<nombre_faces;f++)
            {
                //La face en question
                QVector<int> face = facesTab[f];
                QVector<QVector<int>> sommets_adjacents={};

                //Parcourt des faces adjacentes pour recuperer les 3 sommets adjacents (ceux pas sur le triangle)
                for(int adj=0;adj<3;adj++)
                {
                    //Sommet 0
                    int sommet_0=facesTab[adjFacesTab[f][adj]][0];
                    if(sommet_0!=face[0] && sommet_0!=face[1] && sommet_0!=face[2]) { sommets_adjacents.push_back({sommet_0,adjFacesTab[f][adj]});}
                    //Sommet 1
                    int sommet_1=facesTab[adjFacesTab[f][adj]][1];
                    if(sommet_1!=face[0] && sommet_1!=face[1] && sommet_1!=face[2]) { sommets_adjacents.push_back({sommet_1,adjFacesTab[f][adj]});}
                    //Sommet 2
                    int sommet_2=facesTab[adjFacesTab[f][adj]][2];
                    if(sommet_2!=face[0] && sommet_2!=face[1] && sommet_2!=face[2]) { sommets_adjacents.push_back({sommet_2,adjFacesTab[f][adj]});}
                }

                //Verifie si un des sommets ajacents est dans le cercle circonscrit
                int nb_sommets_adjacents=sommets_adjacents.size(); //Normalement tjs 3
                for(int k=0;k<nb_sommets_adjacents;k++)
                {
                    int indice_sommet=sommets_adjacents[k][0];
                    int indice_face_adj=sommets_adjacents[k][1];

                    //Si le sommet est le point infini, ne pas le prendre en compte
                    if(indice_sommet!=Pinfini_index)
                    {
                        //Si le sommet est strictement dans le cercle circonscrit
                        if(in_circumscribed_circle_test({&vertexTab[facesTab[f][0]],&vertexTab[facesTab[f][1]],&vertexTab[facesTab[f][2]]},vertexTab[indice_sommet])>0.0)
                        {
                            bool already_in_pile=false;
                            //Parcourt de la pile pour voir si la paire y est déjà
                            int pile_size=pile.size();
                            for(int p=0;p<pile_size;p++)
                            {
                                if(pile[p][0]==indice_face_adj && pile[p][1]==f)
                                {
                                    already_in_pile=true;
                                }
                                if(pile[p][0]==f && pile[p][1]==indice_face_adj)
                                {

                                }
                            }
                            //Si elle n'y est pas, on l'ajoute
                            if(already_in_pile==false){pile.push_back({f,indice_face_adj});}
                        }
                    }
                }
            }

            int taille_pile = pile.size();
            //Cas ou on est de Delaunay
            if(taille_pile==0){again=false;}
            //Sinon on vide la pile
            else
            {
                for(int elt=0; elt<taille_pile; elt++)
                {
                    //std::cout << "TRIANGLE 289 " << facesTab[289][0] << " "<< facesTab[289][1] << " "<< facesTab[289][2] << std::endl;
                    //On flipe l'arete commune
                    QVector<int> paire_faces=pile[elt];
                    flip_edge(paire_faces[0],paire_faces[1]);
                }
                pile={};
            }
            //again=false;
            // Affichage de l'avancement
            std::cout << "Nombre de sommets ajoutes : " << i << "/" << nb_true_vertices << std::endl;
        }
    }

    //Puis on redirige les 3 sommets du triangle fictif vers le sommet infini
    int nombre_faces = facesTab.size();
    for(int f=0;f<nombre_faces;f++)
    {
        for(int sommet=0;sommet<3;sommet ++)
        {
            if(facesTab[f][sommet]==P1_index ||facesTab[f][sommet]==P2_index ||facesTab[f][sommet]==P3_index)
            {
                facesTab[f][sommet]=Pinfini_index;
            }
        }
    }

    nb_faces=facesTab.size();
    std::cout << "Nombre sommets triangulation :" << vertexTab.size() << std::endl;
    std::cout << "Nombre triangles triangulation :" <<  facesTab.size() << std::endl;
}










// FONCTIONS ANNEXES

// Pt1, Pt2 et Pt3 sont des points 2D (x,y,z=0)
// Renvoie +1 si pt1,pt2 et pt3 dans le sens trigo
// Renvoie -1 si pt1,pt2 et pt3 dans le sens horaire
// Renvoie 0 si pt1,pt2 et pt3 alignés
double orientation_test(Point& pt1,Point& pt2,Point& pt3)
{
    Point vec_21(pt2.x()-pt1.x(),pt2.y()-pt1.y(),pt2.z()-pt1.z());
    Point vec_23(pt3.x()-pt1.x(),pt3.y()-pt1.y(),pt3.z()-pt1.z());
    Point prod_vect(0.0,0.0,vec_23.x()*vec_21.y()-vec_23.y()*vec_21.x());

    if(prod_vect.z()>0.0){return 1.0;}
    else if (prod_vect.z()<0.0){return -1.0;}
    else {return 0.0;}
}


// Triangle est un triangle 2D (x,y,z=0) trié bien orienté
// Renvoie +1 si pt dans triangle
// Renvoie -1 si pt hors du triangle
// Renvoie 0 si pt au bord du triangle
double in_triangle_test(QVector<Point*> triangle, Point& pt)
{
    Point* A = triangle[0];
    Point* B = triangle[1];
    Point* C = triangle[2];
    Point* P = &pt;

    //(AB vect AP).z
    double AB_AP= (B->x()-A->x())*(P->y()-A->y()) - (B->y()-A->y())*(P->x()-A->x());
    //(BC vect BP).z
    double BC_BP= (C->x()-B->x())*(P->y()-B->y()) - (C->y()-B->y())*(P->x()-B->x());
    //(CA vect CP).z
    double CA_CP= (A->x()-C->x())*(P->y()-C->y()) - (A->y()-C->y())*(P->x()-C->x());

    //Resultat test
    if(AB_AP<0.0 || BC_BP<0.0 || CA_CP<0.0) {return -1.0;}
    else if(AB_AP==0.0 || BC_BP==0.0 || CA_CP==0.0) {return 0.0;}
    else {return 1.0;}
}


// Renvoie +1 si pt dans cercle circonscrit
// Renvoie -1 si pt hors du cercle circonscrit ou au bord
double in_circumscribed_circle_test(QVector<Point*> triangle, Point& pt)
{
    //if(pt.z()<-0.5){return -1.0;}
    Point* A = triangle[0];
    Point* B = triangle[1];
    Point* C = triangle[2];
    Point* P = &pt;

    Point Ap(A->x(),A->y(), A->x()*A->x()+A->y()*A->y());
    Point Bp(B->x(),B->y(), B->x()*B->x()+B->y()*B->y());
    Point Cp(C->x(),C->y(), C->x()*C->x()+C->y()*C->y());
    Point Pp(P->x(),P->y(), P->x()*P->x()+P->y()*P->y());

    Point ApBp(Bp.x()-Ap.x(),Bp.y()-Ap.y(),Bp.z()-Ap.z());
    Point ApCp(Cp.x()-Ap.x(),Cp.y()-Ap.y(),Cp.z()-Ap.z());
    Point ApPp(Pp.x()-Ap.x(),Pp.y()-Ap.y(),Pp.z()-Ap.z());

    //ApBp vect ApCp
    double nx=ApBp.y()*ApCp.z()-ApBp.z()*ApCp.y();
    double ny=ApBp.z()*ApCp.x()-ApBp.x()*ApCp.z();
    double nz=ApBp.x()*ApCp.y()-ApBp.y()*ApCp.x();
    Point n(nx,ny,nz);

    //Produit scalaire avec ApPp
    double scal = n.x()*ApPp.x()+n.y()*ApPp.y()+n.z()*ApPp.z();
    if(scal<0.0){return 1.0;}
    else{return -1.0;}
}
