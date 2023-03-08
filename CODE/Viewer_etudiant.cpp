#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "draw.h"
#include "Viewer_etudiant.h"
#include "wavefront.h"


using namespace std;


///Procedure pour initialiser un cube avec 6 quadrilateres
void ViewerEtudiant::init_cubo_q()
{
    //On choisit les primitives
    m_cubo_q=Mesh(GL_TRIANGLE_STRIP);
    ///On cree chaqu'une des differentes faces
    //1 face
    //normale face
    m_cubo_q.normal(0,-1,0);
    //sommets de la face
    m_cubo_q.vertex(-1,-1,-1);
    m_cubo_q.vertex(1,-1,-1);
    m_cubo_q.vertex(-1,-1,1);
    m_cubo_q.vertex(1,-1,1);
    //On demande un nouveau strip
    m_cubo_q.restart_strip();
    //2 face
    m_cubo_q.normal(0,1,0);

    m_cubo_q.vertex(1,1,-1);
    m_cubo_q.vertex(-1,1,-1);
    m_cubo_q.vertex(1,1,1);
    m_cubo_q.vertex(-1,1,1);

    m_cubo_q.restart_strip();

    //3 face
    m_cubo_q.normal(1,0,0);

    m_cubo_q.vertex(1,-1,1);
    m_cubo_q.vertex(1,-1,-1);
    m_cubo_q.vertex(1,1,1);
    m_cubo_q.vertex(1,1,-1);

    m_cubo_q.restart_strip();

    //4 face
    m_cubo_q.normal(-1,0,0);

    m_cubo_q.vertex(-1,-1,-1);
    m_cubo_q.vertex(-1,-1,1);
    m_cubo_q.vertex(-1,1,-1);
    m_cubo_q.vertex(-1,1,1);

    m_cubo_q.restart_strip();

    //5 face
    m_cubo_q.normal(0,0,1);

    m_cubo_q.vertex(-1,-1,1);
    m_cubo_q.vertex(1,-1,1);
    m_cubo_q.vertex(-1,1,1);
    m_cubo_q.vertex(1,1,1);

    m_cubo_q.restart_strip();

    //6 face
    m_cubo_q.normal(0,0,-1);

    m_cubo_q.vertex(1,-1,-1);
    m_cubo_q.vertex(-1,-1,-1);
    m_cubo_q.vertex(1,1,-1);
    m_cubo_q.vertex(-1,1,-1);

    m_cubo_q.restart_strip();
	}

///Procedure pour initialiser un cube avec une structure indexée
void ViewerEtudiant::init_cubo_s()
{
    //On choisit les primitives
    m_cubo_s=Mesh(GL_TRIANGLE_STRIP);
    //Sommets du cube               1           2          3         4         5         6          7        8
    static float point[8][3] = {{-1,-1,-1}, {1,-1,-1}, {1,-1,1}, {-1,-1,1}, {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1,1,1}};
    //Faces du cube
    static int face[6][4] = {{0,1,2,3}, {5,4,7,6}, {2,1,5,6}, {0,3,7,4}, {3,2,6,7}, {1,0,4,5}};
    //Normales des faces
    static float normal[6][3] = {{0,-1,0}, {0,1,0}, {1,0,0}, {-1,0,0}, {0,0,1}, {0,0,-1}};

    ///Maintenant on va placé la normale associé a chaqu'une des faces et les faces
    for (int i=0; i<6; i++)
        {
            //Normal de la face i
            m_cubo_s.normal(normal[i][0], normal[i][1], normal[i][2]);
            //Sommets de la face i
            m_cubo_s.texcoord(0,0);
            m_cubo_s.vertex(point[face[i][0]][0], point[face[i][0]][1], point[face[i][0]][2]);
            m_cubo_s.texcoord(1,0);
            m_cubo_s.vertex(point[face[i][1]][0], point[face[i][1]][1], point[face[i][1]][2]);
            m_cubo_s.texcoord(0,1);
            m_cubo_s.vertex(point[face[i][3]][0], point[face[i][3]][1], point[face[i][3]][2]);
            m_cubo_s.texcoord(1,1);
            m_cubo_s.vertex(point[face[i][2]][0], point[face[i][2]][1], point[face[i][2]][2]);
            m_cubo_s.restart_strip();
        }
}

///Procedure pour initialiser un cylindre
void ViewerEtudiant::init_cylinder()
{
    //Choix des primitives
    m_cylinder=Mesh(GL_TRIANGLE_STRIP);

    const int div = 25.0;   //Nombre de fois que le cylindre va etre coupe
    float alpha;            //Angle qu'on va faire varie
    float step = 2.0*M_PI/div;//Le pas qu'on va effectuer pour la creation des differents points

    ///Maintenant on va placé les normales et les points
    for (int i=0;i<=div;i++)
    {
        //Variation de l'angle alpha de 0 a 2PI
        alpha = i * step;
        //Disque du bas
        m_cylinder.normal(Vector(cos(alpha), 0, sin(alpha)));
        m_cylinder.texcoord(float(i)/div,0);
        m_cylinder.vertex(Point(cos(alpha), -1,sin(alpha)));
        //Disque du haut
        m_cylinder.normal(Vector(cos(alpha), 0, sin(alpha)));
        m_cylinder.texcoord(float(i)/div,1);
        m_cylinder.vertex(Point(cos(alpha),1,sin(alpha)));
    }
}

///Procedure pour initialiser un cone
void ViewerEtudiant::init_cone()
{
    //Choix des primitives
    m_cone=Mesh(GL_TRIANGLE_STRIP);

    const int div = 25;     //Nombre de fois que le cone va etre coupe
    float alpha;            //Angle qu'on va faire varie
    float step = 2.0 * M_PI / (div);//Le pas qu'on va faire entre chaqu'un des points

    ///Maintenant on va placé les normales et les points
    for (int i=0;i<=div;i++)
    {
        //Variation de l'angle alpha de 0 a 2PI
        alpha = i * step;
        //Cercle du bas
        m_cone.normal(Vector( cos(alpha)/sqrtf(2.f),1.f/sqrtf(2.f),sin(alpha)/sqrtf(2.f) ));
        m_cone.texcoord(float(i)/div,0);
        m_cone.vertex( Point( cos(alpha), 0 ,sin(alpha)));
        //Pointe du cône
        m_cone.normal(Vector( cos(alpha)/sqrtf(2.f),1.f/sqrtf(2.f),sin(alpha)/sqrtf(2.f) ));
        m_cone.texcoord(float(i)/div,1);
        m_cone.vertex( Point(0, 1, 0) );
    }
}

///Procedure pour initialiser un disque
void ViewerEtudiant::init_disque()
{
    //Choix des primitives
    m_disque=Mesh(GL_TRIANGLE_STRIP);

    const int div = 40;     //Nombre de fois que le disque va etre coupe
    float alpha;            //Angle qu'on va faire varie
    float step = 2.0 * M_PI / (div);//Le pas qu'on va faire entre chaque point

    ///Maintenant on va placé les normales et les points
    for (int i=0;i<=div;i++)
    {
        //Variation de l'angle de 0 a 2PI
        alpha = float(i) * step;
        //Contour du disque
        m_disque.normal(Vector( cos(alpha)/sqrtf(2.f),1.f/sqrtf(2.f),sin(alpha)/sqrtf(2.f) ));
        m_disque.texcoord(float(i)/div,0);
        m_disque.vertex( Point( sin(alpha), 0,cos(alpha) ));
        //Centre du disque
        m_disque.normal(Vector( cos(alpha)/sqrtf(2.f),1.f/sqrtf(2.f),sin(alpha)/sqrtf(2.f) ));
        m_disque.texcoord(float(i)/div,1);
        m_disque.vertex( Point(0, 0, 0) );

    }
}

///Procedure pour initialiser une sphere
void ViewerEtudiant::init_sphere()
{
    //Choix des primitives
    m_sphere = Mesh(GL_TRIANGLE_STRIP);

    float beta, alpha1, alpha2;     //Angles qu'on va faire varie
    const int divBeta = 16;         //Nombre de fois qu'on va coupe Beta
    const int divAlpha = divBeta/2;//Nombre de fois qu'on va coupe Alpha

    ///On fait varie alpha1 et alpha2, on fait la superposition des cercles
    for(int i=0; i<divAlpha; i++)
    {
        alpha1 = -0.5f * M_PI + float(i) * M_PI / divAlpha;
        alpha2 = -0.5f * M_PI + float(i+1) * M_PI / divAlpha;
        ///Variation de beta, on fait le dessin des sommets d’un cercle
        for(int j=0; j<=divBeta; j++)
        {
            beta = float(j) * 2.f * M_PI / (divBeta);
            m_sphere.normal(Vector(cos(alpha1)*cos(beta), sin(alpha1), cos(alpha1)*sin(beta)));
            m_sphere.texcoord(beta/(2.0*M_PI),0.5 + alpha1/M_PI);
            m_sphere.vertex(Point(cos(alpha1)*cos(beta), sin(alpha1), cos(alpha1)*sin(beta)));

            m_sphere.normal(Vector(cos(alpha2)*cos(beta), sin(alpha2), cos(alpha2)*sin(beta)));
            m_sphere.texcoord(beta/(2.0*M_PI),0.5 + alpha2/M_PI);
            m_sphere.vertex(Point(cos(alpha2)*cos(beta), sin(alpha2), cos(alpha2)*sin(beta)));
        }
        //A chaque fois qu'on finit un cercle on demande de nouveau un strip
        m_sphere.restart_strip();
    }
}

///Fonction qui retourne la normale au point (i,j) d'une image
Vector terNormal(const Image& im, const int i, const int j)
{
    int ip = i-1;
    int in = i+1;
    int jp = j-1;
    int jn = j+1;
    //On calcule les points autour de (i,j)
    Vector a( ip, im(ip, j).r, j );
    Vector b( in, im(in, j).r, j );
    Vector c( i, im(i, jp).r, jp );
    Vector d( i, im(i, jn).r, jn );
    //On calcule les vecteurs unitaires
    Vector ab = normalize(b - a);
    Vector cd = normalize(d - c);
    //Puis on renvoi le produit de ces deux vecteurs qui va etre la normal au point (i,j)
    Vector n = cross(ab,cd);
    return n;
}

///Procedure qui initialise un terrain a partir d'un height map en niveau de gris
void ViewerEtudiant::init_terrain(Mesh& m_terrain, const Image& im)
{
    //Choix des primitives
    m_terrain = Mesh(GL_TRIANGLE_STRIP);

    ///On parcour la largeur et longeur de l'image pixel par pixel pour place les normales et les points
    for(int i=1;i<im.width()-2;i++)
    {
        for(int j=1;j<im.height()-1;j++)
        {
            //Calcule de la normale a (i+1,j)
            m_terrain.normal(terNormal(im,i+1,j));
            m_terrain.texcoord(float(i+1)/im.width(),float(j)/im.height());
            //On definie y grace a la couleur de l'image
            m_terrain.vertex(Point(i+1, 25.f*im(i+1,j).r,j));
            //Calcule de la normale a (i,j)
            m_terrain.normal(terNormal(im,i,j));
            m_terrain.texcoord(float(i)/im.width(),float(j)/im.height());
            //Meme chose ici
            m_terrain.vertex(Point(i, 25.f*im(i,j).r,j));
        }
        //On demande un nouveau strip au final de chaque iteration
        m_terrain.restart_strip();
    }

}

///Procedure qui initialise un billboard
void ViewerEtudiant::init_billboard()
{
    //Choix des primitives
    m_quad = Mesh(GL_TRIANGLE_STRIP);
    ///Construction des quad de tal facon que l'image soit toujours face a la camera
    //1 quad
    m_quad.texcoord(0,0);
    m_quad.vertex(-1,-1,0);
    m_quad.texcoord(0,1);
    m_quad.vertex(-1,1,0);
    m_quad.texcoord(1,0);
    m_quad.vertex(1,-1,0);
    m_quad.texcoord(1,1);
    m_quad.vertex(1,1,0);
    //On demande un nouveau strip apres la creation du premier quad
    m_quad.restart_strip();
    //2 quad
    m_quad.texcoord(0,0);
    m_quad.vertex(0,-1,-1);
    m_quad.texcoord(0,1);
    m_quad.vertex(0,1,-1);
    m_quad.texcoord(1,0);
    m_quad.vertex(0,-1,1);
    m_quad.texcoord(1,1);
    m_quad.vertex(0,1,1);
}

///Procedure qui initialise un cubemap
void ViewerEtudiant::init_map()
{
    //Choix des primitives
    m_cubo_map=Mesh(GL_TRIANGLE_STRIP);
    //Sommets du cube               1           2           3       4           5          6        7         8
    static float point[8][3] = {{-1,-1,-1}, {1,-1,-1}, {1,-1,1}, {-1,-1,1}, {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1,1,1}};
    //Faces du cube
    static int face[6][4] = {{0,1,2,3}, {5,4,7,6}, {2,1,5,6}, {0,3,7,4}, {3,2,6,7}, {1,0,4,5}};
    //Normales des faces
    static float normal[6][3] = {{0,-1,0}, {0,1,0}, {1,0,0}, {-1,0,0}, {0,0,1}, {0,0,-1}};
    //Coordonees des textcoord pour chaqu'une des faces
    static vec2 txc[6][4]={
                            //Bas
                            {{1/4.f,0.f},{2/4.f,0.f},{2/4.f,1/3.f},{1/4.f,1/3.f}},

                            //Haut
                            {{2/4.f,1.f},{1/4.f,1.f},{1/4.f,2/3.f},{2/4.f,2/3.f}},

                            //Right
                            {{1/4.f,1/3.f},{3/4.f,1/3.f},{3/4.f,2/3.f},{2/4.f,2/3.f}},

                            //Left
                            {{0.f,1/3.f},{1/4.f,1/3.f},{1/4.f,2/3.f},{0.f,2/3.f}},


                            //Front
                            {{1/4.f,1/3.f},{2/4.f,1/3.f},{2/4.f,2/3.f},{1/4.f,2/3.f}},


                            //Back
                            {{1/4.f,1/3.f},{1.f,1/3.f},{1.f,2/3.f},{3/4.f,2/3.f}}

                            };
    ///Maintenant on va placé la normale associé a chaqu'une des faces et les faces
    for (int i=0; i<6; i++)
    {
        //Normal de la face i
        m_cubo_map.normal(normal[i][0],-1*normal[i][1], normal[i][2]);
        //Sommets de la face i
        m_cubo_map.texcoord(txc[i][0].x,txc[i][0].y);
        //Ici on inverse la coordonee z pour qu'on voit l'interieur du cube
        m_cubo_map.vertex(point[face[i][0]][0], point[face[i][0]][1],-point[face[i][0]][2]);
        m_cubo_map.texcoord(txc[i][1].x,txc[i][0].y);
        m_cubo_map.vertex(point[face[i][1]][0], point[face[i][1]][1], -point[face[i][1]][2]);
        m_cubo_map.texcoord(txc[i][3].x,txc[i][3].y);
        m_cubo_map.vertex(point[face[i][3]][0], point[face[i][3]][1],-point[face[i][3]][2]);
        m_cubo_map.texcoord(txc[i][2].x,txc[i][2].y);
        m_cubo_map.vertex(point[face[i][2]][0], point[face[i][2]][1], -point[face[i][2]][2]);
        //On demande un nouveau strip
        m_cubo_map.restart_strip();
    }
}

///Procedure qui cree la forme de la figure qui va subir l'extruction
void ViewerEtudiant::creation_forme()
{
    // Nombre de points de la silhouette 2D
    NBPT = 6;
    /// Points de la silhouette 2D
    Point objt_p[NBPT];
    objt_p[0] =Point(0.0,0.0,0.0);
    objt_p[1] =Point(0.0,1.0,0.0);
    objt_p[2] = Point(1.0,2.0,0.0);
    objt_p[3] = Point(2.0, 1.0, 0.0);
    objt_p[4] = Point(2.0, 0.0, 0.0);
    objt_p[5] = Point(0.0, 0.0, 0.0);
    //Taille de l'extruction
    float prof = 5.0;

    //On applique la matrice de translation pour la coordonne z
    for(int i=0; i < 2; i++){
        for(int j=0; j < NBPT; j++){
            objt_v[i][j].x =objt_p[j].x;
            objt_v[i][j].y =objt_p[j].y;
            objt_v[i][j].z =objt_p[j].z+i*prof;
        }
    }

    //On initialise tout les normales nulles
    for(int i=0; i<2; i++)
        for(int j=0; j<NBPT; j++)
            objt_vn[i][j] = Vector(0, 0, 0);

    //On initialise les normales
    for(int i=0; i<NBPT-1; i++){
        Vector a, b, vntmp;
        a = normalize(objt_v[0][i] - objt_v[1][i+1]);
        b = normalize(objt_v[0][i] - objt_v[1][i]);
        vntmp = cross(a, b);
        //On transmet cette normal au quatres sommets de la face
        objt_vn[0][i] = vntmp + objt_vn[0][i];
        objt_vn[1][i] = vntmp + objt_vn[1][i];
        objt_vn[1][i+1] = vntmp + objt_vn[1][i+1];
        objt_vn[0][i+1] = vntmp + objt_vn[0][i+1];
    }
    //Normale a un sommet est la moyenne des 4 sommets de la face
    for(int i=0; i<2; i++){
        for(int j=0; j<NBPT; j++){
            float q = 4.0f;
            if (j == NBPT-1)
            q = 2.0f;
            objt_vn[i][j] = objt_vn[i][j] / q;
            }
        }
}

///Procedure qui initialise l'objet
void ViewerEtudiant::init_obj_extru()
{
    //Choix des primitives
    m_extru = Mesh(GL_TRIANGLES);

    ///On parcours chaque point pour place la normale et les coordonnes des points
    for(int j=0; j<NBPT-1; j++)
    {

        //1er triangle
        m_extru.texcoord(0,0);
        m_extru.normal(objt_vn[1][j]);
        m_extru.vertex(objt_v[1][j].x, objt_v[1][j].y,objt_v[1][j].z);

        m_extru.texcoord(1,1);
        m_extru.normal(objt_vn[0][j+1]);
        m_extru.vertex(objt_v[0][j+1].x, objt_v[0][j+1].y,objt_v[0][j+1].z);

        m_extru.texcoord(1,0);
        m_extru.normal(objt_vn[0][j]);
        m_extru.vertex(objt_v[0][j].x,objt_v[0][j].y, objt_v[0][j].z);

        //2eme triangle

        m_extru.texcoord(0,0);
        m_extru.normal(objt_vn[1][j]);
        m_extru.vertex(objt_v[1][j].x, objt_v[1][j].y, objt_v[1][j].z);

        m_extru.texcoord(0,1);
        m_extru.normal(objt_vn[1][j+1]);
        m_extru.vertex(objt_v[1][j+1].x, objt_v[1][j+1].y, objt_v[1][j+1].z);

        m_extru.texcoord(1,1);
        m_extru.normal(objt_vn[0][j+1]);
        m_extru.vertex(objt_v[0][j+1].x, objt_v[0][j+1].y, objt_v[0][j+1].z);
    }
}

///Procedure qui initialise un quad ou va etre present une texture animee
void ViewerEtudiant::init_iman()
{

    m_iman = Mesh(GL_TRIANGLE_STRIP);

    m_iman.normal(  0, 0, 1 );

    m_iman.texcoord(0,0 );
    m_iman.vertex(-1, -1, 0 );

    m_iman.texcoord(1,0);
    m_iman.vertex(  1, -1, 0 );

    m_iman.texcoord(0,1);
    m_iman.vertex( -1, 1, 0 );

    m_iman.texcoord( 1,1);
    m_iman.vertex(  1,  1, 0 );
}

///Procedure qui initialise un quad ou va etre present une texture animee
void ViewerEtudiant::init_iman2()
{

   m_iman2 = Mesh(GL_TRIANGLE_STRIP);

    m_iman2.normal(  0, 0, 1 );

    m_iman2.texcoord(0,0 );
    m_iman2.vertex(-1, -20, 0 );

    m_iman2.texcoord(1,0);
    m_iman2.vertex(  m_terrainAlti.width()/2, -20, 0 );

    m_iman2.texcoord(0,1);
    m_iman2.vertex( -1, 3, 0 );

    m_iman2.texcoord( 1,1);
    m_iman2.vertex(  m_terrainAlti.width()/2,  3, 0 );;
}

/*
 * Fonction dans laquelle les initialisations sont faites.
 */
int ViewerEtudiant::init()
{
    Viewer::init();

    m_camera.lookat( Point(0,0,0), 150 );


    /// Chargement des textures
    //Textures pour les formes de base
    m_cubo_texture=read_texture(0,"data/nuclear.jpg") ;
    m_globo_texture=read_texture(0,"data/metal.jpg") ;
    //Textures pour le terrain
    m_terrainAlti = read_image("data/terrain/terreno.jpg");
    m_terrain_texture= read_texture(0,"data/terrain/terrenotex.png");
    //Texture pour les billboards et le cube-map
    m_arbre_texture= read_texture(0,"data/billboard/arbol.png");
    m_palmera_texture= read_texture(0,"data/billboard/palmera.png");
    m_map_texture=read_texture(0,"data/cubemap/map.jpg");
    //Textures pour la premiere textures animee
    m_iman_texture1=read_texture(0,"data/shark/shark1.png");
    m_iman_texture2=read_texture(0,"data/shark/shark2.png");
    m_iman_texture3=read_texture(0,"data/shark/shark3.png");
    m_iman_texture4=read_texture(0,"data/shark/shark4.png");
    m_iman_texture5=read_texture(0,"data/shark/shark5.png");
    m_iman_texture6=read_texture(0,"data/shark/shark6.png");
    m_iman_texture7=read_texture(0,"data/shark/shark7.png");
    m_iman_texture8=read_texture(0,"data/shark/shark8.png");
    m_iman_texture9=read_texture(0,"data/shark/shark9.png");
    m_iman_texture10=read_texture(0,"data/shark/shark10.png");
    //Textures pour la seconde texture animee
    m_iman2_texture1=read_texture(0,"data/cascada/cascada1.png");
    m_iman2_texture2=read_texture(0,"data/cascada/cascada2.png");
    m_iman2_texture3=read_texture(0,"data/cascada/cascada3.png");
    m_iman2_texture4=read_texture(0,"data/cascada/cascada4.png");
    m_iman2_texture5=read_texture(0,"data/cascada/cascada5.png");
    m_iman2_texture6=read_texture(0,"data/cascada/cascada6.png");
    m_iman2_texture7=read_texture(0,"data/cascada/cascada7.png");
    m_iman2_texture8=read_texture(0,"data/cascada/cascada8.png");
    m_iman2_texture9=read_texture(0,"data/cascada/cascada9.png");
    //Texture pour l'extruction
    m_tent_texture=read_texture(0,"data/tent.png");
    //Texture pour les objets importer
    m_ele_texture=read_texture(0,"data/elefantefull.png");
    m_turtle_texture=read_texture(0,"data/turtle.jpg");
    m_ufo_texture=read_texture(0,"data/ufo.jpg");

    //// Appel des fonctions init_votreObjet pour creer les Mesh
    ///Initialisation des formes de base
    init_sphere();
    init_cone();
    init_cube();
    init_disque();
    init_cylinder();
    ///Initialisation du terrain
    init_terrain(m_terrain, m_terrainAlti);
    ///Initialisation du cube-map et des billboards
    init_map();
    init_billboard();
    ///Initialisation pour l'extruction
    creation_forme();
    init_obj_extru();
    ///Initialisation pour les textures animees
    init_iman();
    init_iman2();
    ///Initialisation pour les objets importer
    obj1= read_mesh("data/UFO.obj");
    obj2= read_mesh("data/elefante.obj");
    obj3= read_mesh("data/turtle.obj");

    ///Initialisation des coordonnes des billboards
    float i,j,x,y;
    //Arbre
    for (int k=0; k<Narbre; k++)
    {
        do{
        //On prend des nombres au hasard entre 0 et la longeur et larguer de l'image
        i =rand() % m_terrainAlti.width();

        j = rand() % m_terrainAlti.height();
        //Pour positioner les arbres seulement dans une zone, ou 25.f*m_terrainAlti(i, j).r+1 est la position Y des bilboards
        }while(i>192/4 || j>(-i+48) || 25.f*m_terrainAlti(i, j).r+1<=5);
        //On garde les coordonees dans un tableau
        tab[k]=vec2(i,j);

    }
    //Palmiers
    for (int w=0; w<Npalmera; w++)
    {
        do{
        //Meme chose que les arbres
        x =rand() % m_terrainAlti.width();

        y = rand() % m_terrainAlti.height();
        }while(25.f*m_terrainAlti(x, y).r+1<=5 || x>180 || y>95);

        tab2[w]=vec2(x,y);

    }

    return 0;
}

/*
 * Exemple de definition de fonction permettant l affichage
 * de 'votreObjet' subissant la Transform T
 */
 ///Procedure pour dessiner un cube avec 6 quadrilateres
 void ViewerEtudiant::draw_cubo_q(const Transform& T)
{
    gl.model( T );
    gl.draw( m_cubo_q);
}

///Procedure pour dessiner un cube avec structure indexée
void ViewerEtudiant::draw_cubo_s(const Transform& T)
{
    gl.alpha(0.5f);
    //Appel au texture
    gl.texture(m_cubo_texture);
    //Appel a la transformation
    gl.model( T );
    //Appel au mesh
    gl.draw( m_cubo_s);
}

///Procedure pour dessiner un cylindre
void ViewerEtudiant::draw_cylinder(const Transform& T)
{
    gl.texture(m_cubo_texture);
    gl.model( T );
    gl.draw( m_cylinder);
}

///Procedure pour dessiner un cone
void ViewerEtudiant::draw_cone(const Transform& T)
{
    gl.texture(m_cubo_texture);
    gl.model( T );
    gl.draw( m_cone);
}

///Procedure pour dessiner une sphere
void ViewerEtudiant::draw_sphere(const Transform& T)
{
    gl.texture(m_globo_texture);
    gl.model( T );
    gl.draw( m_sphere);
}

///Procedure pour dessiner une Bombe
void ViewerEtudiant::dessineObjet(const Transform& T)
{
    gl.alpha(0.5f);

    gl.texture(m_globo_texture);
    gl.model(T*Scale(3,3,6));
    gl.draw( m_sphere);

    gl.texture(m_cubo_texture);
    gl.model(T*Scale(1,2.3,1)*Translation(0,0,4));
    gl.draw(m_cube);

    gl.model(T*Scale(2.3,1,1)*Translation(0,0,4));
    gl.draw(m_cube);

    gl.texture(m_globo_texture);
    gl.model(T*Scale(2,2,4)*Translation(0,0,2)*RotationX(270));
    gl.draw(m_cone);

    gl.model(T*Scale(4,4,3)*Translation(0,0,3)*RotationX(270));
    gl.draw(m_cone);

    gl.model(T*Scale(4,4,1.8)*Translation(0,0,6)*RotationX(270));
    gl.draw(m_cylinder);

    gl.model(T*Scale(4,4,4)*Translation(0,0,3.1)*RotationX(270));
    gl.draw(m_disque);

}

///Procedure pour dessiner un terrain
void ViewerEtudiant::draw_terrain(const Transform &T)
{
    gl.alpha(0.5f);
    gl.model( T );
    gl.texture(m_terrain_texture);
    gl.draw( m_terrain);
}

///Procedure pour dessiner un arbre
void ViewerEtudiant::draw_billboard(const Transform &T)
{
    gl.alpha(0.5f);
    Transform tt=T;
    gl.model(tt );
    gl.texture(m_arbre_texture);
    gl.draw( m_quad);

    gl.model(tt *RotationY(180));
    gl.draw( m_quad);
}

///Procedure pour dessiner un palmier
void ViewerEtudiant::draw_billboard2(const Transform &T)
{
    gl.alpha(0.5f);
    Transform tt=T*Scale(3,3,3);
    gl.model(tt );
    gl.texture(m_palmera_texture);
    gl.draw( m_quad);

    gl.model(tt *RotationY(180));
    gl.draw( m_quad);
}

///Procedure pour dessiner le cube-map
void ViewerEtudiant::draw_map(const Transform &T)
{
    gl.alpha(0.5f);
    gl.model(T);
    gl.texture(m_map_texture);
    gl.draw(m_cubo_map);

}

///Procedure pour dessiner une tente
void ViewerEtudiant::draw_obj_extru(const Transform& T)
{
    gl.alpha(0.5f);
    gl.texture(m_tent_texture);
    gl.model(T*Scale(0.8,0.8,0.8));
    gl.draw(m_extru);
}

///Procedure pour dessiner une texture animee (requin)
void ViewerEtudiant::draw_iman(const Transform& T)
{
    //Nombre d'images
    int x=10;
    //Le compteur est une valeur entier qui va augmenter donc selon le resultat
    //du compteur modulo x, les images vont varier et vont changer au cours du temps
    if(compteur%x==0)gl.texture(m_iman_texture1);

    if(compteur%x==1)gl.texture(m_iman_texture2);

    if(compteur%x==2)gl.texture(m_iman_texture3);

    if(compteur%x==3)gl.texture(m_iman_texture4);

    if(compteur%x==4)gl.texture(m_iman_texture5);

    if(compteur%x==5)gl.texture(m_iman_texture6);

    if(compteur%x==6)gl.texture(m_iman_texture7);

    if(compteur%x==7)gl.texture(m_iman_texture8);

    if(compteur%x==8)gl.texture(m_iman_texture9);

    if(compteur%x==9)gl.texture(m_iman_texture10);


    gl.alpha(0.5f);

    gl.model(T);
    gl.draw(m_iman);

    gl.model(T*RotationY(180));
    gl.draw(m_iman);
}

///Procedure pour dessiner un UFO
void  ViewerEtudiant::draw_obj1(const Transform& T)
{
     gl.alpha(0.5f);
    gl.model(T);
    gl.texture(m_ufo_texture);
    gl.draw(obj1);
}

///Procedure pour dessiner des elephants au sommet d'une tortue
void  ViewerEtudiant::draw_obj2(const Transform& T)
{
    gl.alpha(0.5f);
    gl.model(T*Scale(5,5,5)*Translation(5,-5,0));
    gl.texture(m_ele_texture);
    gl.draw(obj2);

    gl.model(T*RotationY(90)*Scale(5,5,5)*Translation(5,-5,0));
    gl.draw(obj2);

    gl.model(T*RotationY(180)*Scale(5,5,5)*Translation(5,-5,0));
    gl.draw(obj2);

    gl.model(T*RotationY(270)*Scale(5,5,5)*Translation(5,-5,0));
    gl.draw(obj2);

    gl.model(T*RotationX(270)*Scale(1,1,1)*Translation(0,0,-50));
    gl.texture(m_turtle_texture);
    gl.draw(obj3);

}

///Procedure qui dessine une texture animee (cascades)
void ViewerEtudiant::draw_iman2(const Transform& T)
{
    int x=9;
    if(compteur%x==0)gl.texture(m_iman2_texture1);

    if(compteur%x==1)gl.texture(m_iman2_texture2);

    if(compteur%x==2)gl.texture(m_iman2_texture3);

    if(compteur%x==3)gl.texture(m_iman2_texture4);

    if(compteur%x==4)gl.texture(m_iman2_texture5);

    if(compteur%x==5)gl.texture(m_iman2_texture6);

    if(compteur%x==6)gl.texture(m_iman2_texture7);

    if(compteur%x==7)gl.texture(m_iman2_texture8);

    if(compteur%x==8)gl.texture(m_iman2_texture9);


    Transform tt=T;

    gl.alpha(0.5f);

    gl.model(tt*Translation(-47,-11,49));
    gl.draw(m_iman2);

    gl.model(tt*Translation(49,-11,49)*RotationY(90));
    gl.draw(m_iman2);

    gl.model(tt*Translation(-45,-9,-47)*RotationY(270));
    gl.draw(m_iman2);

    gl.model(tt*Translation(49,-9,-45)*RotationY(180));
    gl.draw(m_iman2);
}

/*
 * Fonction dans laquelle les appels pour les affichages sont effectues.
 */
int ViewerEtudiant::render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    manageCameraLight();

    gl.camera(m_camera);

    /// Appel des fonctions du type 'draw_votreObjet'
    Transform TT=Scale(50, 50, 50);
    draw_map(TT);
    Transform T=Scale(0.5,0.5,0.5)*Translation(-90,-20,-90);


    draw_terrain(T);
    //On parcours le tableau avec les coordonees des arbres
    for (int k=0; k<Narbre; k++){

        //On applique une translation pour pouvoir place l'arbre a la bonne place par rapport au terrain
        Transform Ttt =T*Translation(tab[k].x, 25.f*m_terrainAlti(tab[k].x, tab[k].y).r+1, tab[k].y);
        draw_billboard(Ttt);
    }
    //De meme
    for (int w=0; w<Npalmera; w++){

         Transform Tttt =T*Translation(tab2[w].x, 25.f*m_terrainAlti(tab2[w].x, tab2[w].y).r+1, tab2[w].y);
         draw_billboard2(Tttt);
    }
    draw_obj1(Scale(3,3,3)*Translation(0,11,0));
    draw_obj2(Identity());
    draw_iman(Scale(5,7,5)*Translation(0,6,0));
    draw_iman2(Identity());
    dessineObjet(Scale(2,2,2)*Translation(-8,2,7)*RotationX(240));
    draw_obj_extru(Scale(1,1,1)*Translation(-20,2,-28)*RotationY(50));
    return 1;

}


/*
 * Fonction dans laquelle les mises a jours sont effectuees.
 */
int ViewerEtudiant::update( const float time, const float delta )
{
    // time est le temps ecoule depuis le demarrage de l'application, en millisecondes,
    // delta est le temps ecoule depuis l'affichage de la derniere image / le dernier appel a draw(), en millisecondes.
    float top=time/100;
    int ta=int(top);
    //Le compteur est un entier est va augmenter par rapport au temps
    compteur=ta;


    return 0;
}


/*
 * Constructeur.
 */

ViewerEtudiant::ViewerEtudiant() : Viewer()
{
}


/*
 * Programme principal.
 */
int main( int argc, char **argv )
{
    ViewerEtudiant v;
    v.run();
    return 0;
}
