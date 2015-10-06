#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
    #include "olb2D.hh" // include full template code
#endif

using namespace olb;

// Some C++ libraries wich are for the example and others
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <Magick++.h> 
#include <unistd.h>
#include <thread>  
#include <math.h>
#define PI 3.14159265

using namespace olb; // Pode ter erro aqui: declarando 2 vezes
using namespace olb::descriptors; // accessed in the examples
using namespace olb::graphics; 
using namespace std; // Namespace of standard C++ library

// Definindo a minha lattice
// Aqui eu posso definir tambem outros escalares naturais como
// forças. Para tais coisas memoria deve ser alocada na malha de lattice.
#define LATTICE D2Q9Descriptor
typedef double T;
int nx = 300;
int ny = 300;
int numIter = 250; // numero de iteracoes ou passos
T omega = 1.98; // viscosidade
T r = 30.; // raio do circulo


int main(int argc, char* argv[]){

    std::string ss;
    olbInit(&argc, &argv);
    //Ele pede para inserir a principal parte do codigo aqui, mas oq?
    BlockLattice2D<T, LATTICE> lattice(nx, ny); // Aqui o bloco de malha de lattice  nxXnyX9 é instanciado
    //Tipo de dinanimca, pode-se colocar perda de massa por exemplo
    BGKdynamics<T, LATTICE> bulkDynamics(omega, instances::getBulkMomenta<T, LATTICE>());
    //Deve-se indicar em quais lugares da malha lattice vai ocorrer a dinamica: nesse caso com todos os pontos
    lattice.defineDynamics(0,nx-1,0,ny-1,&bulkDynamics);
    
    //Setar a condicao inicial de equilibrio, espalhando as densidades na malha
    //Nesse caso tera uma distribuicao mais densa num circulo de raio r
    T rho = 1., u[2] = {0.,0.};
    for(int iX = 0; iX<nx; ++iX){
        for(int iY = 0; iY<ny; ++iY){
            lattice.get(iX, iY).iniEquilibrium(rho, u);
        }
    }

    //lattice.get(149, 149).iniEquilibrium(1.001, u);

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif("estado_inicial", lattice.getDataAnalysis().getDivRhoU());

    FILE * pFile_1;
    pFile_1 = fopen ("time_points.txt","w");
    // As colisoes definidas em bulkDynamics sao efetivadas aqui em cada celula
    T rho_varia;
    T lattice_speed_sound = 1/sqrt(3);
    T A_amplitude = 0.001;
    for(int iT = 0; iT <= numIter; ++iT){

        //lattice.collide();
        //Fluxo, true indica que as bordas sao periodicas
        //lattice.stream(true);
        lattice.collideAndStream(true);

        std::ostringstream s0;        
        s0 << lattice.get(nx/2, ny/2)[0];
        std::ostringstream s1;
        s1 << lattice.get(nx/2, ny/2)[1];
        std::ostringstream s2;
        s2 << lattice.get(nx/2, ny/2)[2];
        std::ostringstream s3;
        s3 << lattice.get(nx/2, ny/2)[3];
        std::ostringstream s4;
        s4 << lattice.get(nx/2, ny/2)[4];
        std::ostringstream s5;
        s5 << lattice.get(nx/2, ny/2)[5];
        std::ostringstream s6;
        s6 << lattice.get(nx/2, ny/2)[6];
        std::ostringstream s7;
        s7 << lattice.get(nx/2, ny/2)[7];
        std::ostringstream s8;
        s8 << lattice.get(nx/2, ny/2)[8];
        std::ostringstream sums;
        double sum = lattice.get(nx/2, ny/2)[0] + 
        lattice.get(nx/2, ny/2)[1] +
        lattice.get(nx/2, ny/2)[2] + 
        lattice.get(nx/2, ny/2)[3] +
        lattice.get(nx/2, ny/2)[4] +
        lattice.get(nx/2, ny/2)[5] +
        lattice.get(nx/2, ny/2)[6] +
        lattice.get(nx/2, ny/2)[7] +
        lattice.get(nx/2, ny/2)[8];
        sums << sum;

        fprintf (pFile_1, "%.10f\n", sum);

        /*string points_lattice = "1: " + s0.str() + 
        " -- 2:" + s1.str() + 
        " -- 3:" + s2.str() + 
        " -- 4:" + s3.str() + 
        " -- 5:" + s4.str() + 
        " -- 6:" + s5.str() + 
        " -- 7:" + s6.str() + 
        " -- 8:" + s7.str() + 
        " -- 9:" + s8.str()+ " -- soma:" + sums.str() + "\n\n";*/

        string points_lattice = " -- soma:" + sums.str() + "\n\n";

        //cout << points_lattice << endl;

       // std::string s = std::to_string(iT);
        std::ostringstream s;
        s << iT;
        ss = "images/" + s.str();
        ImageWriter<T> imageWriter("leeloo");
        imageWriter.writePpm(ss, lattice.getDataAnalysis().getDivRhoU(), 0, 0.0001);

       rho_varia = 1. + A_amplitude*sin(2*PI*(lattice_speed_sound/20)*iT);
        //printf("%f\n", rho_varia);
        lattice.get(149, 149).defineRho(rho_varia);
    }
    fclose (pFile_1);

    FILE * pFile;
    pFile = fopen ("space_points_1.txt","w");
    cout << "escrevendo resultado final" << endl;
    for (int i = 149; i <= 299; ++i){
        double sum = lattice.get(i, 149).computeRho();
        fprintf (pFile, "%.10f ", sum);
    }
    fclose (pFile);
}






















