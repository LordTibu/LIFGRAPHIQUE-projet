
#ifndef VIEWER_ETUDIANT_H
#define VIEWER_ETUDIANT_H

#include "Viewer.h"



class ViewerEtudiant : public Viewer
{
public:
    ViewerEtudiant();

    int init();
    int render();
    int update( const float time, const float delta );

protected:

    /// Declaration des Mesh
    //Mesh des formes de bases
    Mesh m_cubo_q;
    Mesh m_cubo_s;
    Mesh m_cylinder;
    Mesh m_cone;
    Mesh m_sphere;
    Mesh m_disque;
    //Mesh du terrain
    Mesh m_terrain;
    //Mesh de l'extruction
    Mesh m_extru;
    //Mesh pour les textures animee
    Mesh m_iman;
    Mesh m_iman2;
    //Mesh pour le cube-map
    Mesh m_cubo_map;
    //Mesh pour les objets importer
    Mesh obj1;
    Mesh obj2;
    Mesh obj3;

    /// Declaration des Textures
    //Texture des formes de base
    GLuint m_cubo_texture;
    GLuint m_globo_texture;
    //Textures pour le terrain
    Image m_terrainAlti;
    GLuint m_terrain_texture;
    //Textures pour le billboard
    GLuint m_arbre_texture;
    GLuint m_palmera_texture;
    //Textures pour le cube-map
    GLuint m_map_texture;
    //Textures pour l'extruction
    GLuint m_tent_texture;
    //Texture pour animation (requin)
    GLuint m_iman_texture1;
     GLuint m_iman_texture2;
      GLuint m_iman_texture3;
       GLuint m_iman_texture4;
        GLuint m_iman_texture5;
         GLuint m_iman_texture6;
          GLuint m_iman_texture7;
           GLuint m_iman_texture8;
            GLuint m_iman_texture9;
             GLuint m_iman_texture10;
    //Textures pour animation (cascade)
    GLuint m_iman2_texture1;
     GLuint m_iman2_texture2;
      GLuint m_iman2_texture3;
       GLuint m_iman2_texture4;
        GLuint m_iman2_texture5;
         GLuint m_iman2_texture6;
          GLuint m_iman2_texture7;
           GLuint m_iman2_texture8;
            GLuint m_iman2_texture9;
    //Textures pour les objets
    GLuint m_ele_texture;
    GLuint m_turtle_texture;
    GLuint m_ufo_texture;

    /// Declaration des fonction de creation de Mesh du type init_votreObjet()
    //init formes de bases
    void init_cubo_q();
    void init_cubo_s();
    void init_cylinder();
    void init_cone();
    void init_sphere();
    void init_disque();
    //init terrain
    void init_terrain(Mesh& m_terrain, const Image& im);
    //init billboard
    void init_billboard();
    //init cube-map
    void init_map();
    //init objet par extruction
    void init_obj_extru();
    void creation_forme();//creation de la silhouette 2D
    //init des textures animees, requin et cascade
    void init_iman();
    void init_iman2();

    ///Declaration des variables et tableaux
    int compteur=0;//Compteur pour les textures animees
    vec2 tab[200];//Tableau pour les coordones des arbres
    vec2 tab2[200];//Tableau pour les coordones des palmiers
    int  NBPT;//Declaration des nombres de points pour la sihoulette 2D
    Point objt_v[100][100];//Tableau pour les vertex de l'objet cree par extruction
    Vector objt_vn[100][100];//Tableau pour les normales de l'objet cree par extruction
    int Narbre=200;//Nombre d'arbes a positionner
    int Npalmera=100;//Nombre de palmiers a positionne

    /// Declaration des fonctions draw_votreObjet(const Transform& T)
    //Draw formes de base
    void draw_cubo_q(const Transform& T);
    void draw_cubo_s(const Transform& T);
    void draw_cylinder(const Transform& T);
    void draw_cone(const Transform& T);
    void draw_sphere(const Transform& T);
    //Draw bombe
    void dessineObjet(const Transform& T);
    //Draw terrain
    void draw_terrain(const Transform& T);
    //Draw billboards, arbre et palmier
    void draw_billboard(const Transform& T);
    void draw_billboard2(const Transform& T);
    //Draw cube-map
    void draw_map(const Transform& T);
    //Draw objet par extruction
    void draw_obj_extru(const Transform& T);
    //Draw texture animee, requin et cascade
    void draw_iman(const Transform& T);
    void draw_iman2(const Transform& T);
    //Draw  objet importer, UFO et elephant sur tortue
    void draw_obj1(const Transform& T);
    void draw_obj2(const Transform& T);
};



#endif
