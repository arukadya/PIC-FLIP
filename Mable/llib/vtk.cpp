//
//  vtk.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "vtk.hpp"

void outputVTK(std::string OutputFileName,myArray3d &Vertices,double dx,int nx,int ny,int nz){
    std::vector<float> origin = {0,0,0};
    std::ofstream writing_file;
    writing_file.open(OutputFileName, std::ios::out);
    std::string writing_text ="# vtk DataFile Version 2.0\nIsosurface\nASCII\nDATASET STRUCTURED_POINTS\n";
    writing_file << writing_text << std::endl;
    writing_file << "DIMENSIONS " << nx <<" "<< ny <<" "<< nz << std::endl;
    writing_file << "ORIGIN " << origin[0] <<" "<< origin[1] <<" "<< origin[2] << std::endl;
    writing_file << "SPACING " << dx <<" "<< dx <<" "<< dx << std::endl;
    writing_file << "POINT_DATA " << Vertices.size << std::endl;
    writing_file << "SCALARS value float 1" << std::endl;
    writing_file << "LOOKUP_TABLE default" << std::endl;
    for(int i=0;i<Vertices.nx;i++){
        for(int j=0;j<Vertices.ny;j++){
            for(int k=0;k<Vertices.nz;k++){
                writing_file << Vertices.value[i][j][k] << std::endl;
            }
        }
    }
    writing_file.close();
}

void outputVolume(std::string OutputFileName,std::vector<double> &volumes){
    //FILE *ofp = fopen(OutputFileName,"w");
    std::ofstream writing_file;
    writing_file.open(OutputFileName, std::ios::out);
    double max = volumes[0];
    double min = volumes[0];
    for(int i=0;i<volumes.size();i++){
        writing_file << i << " " << volumes[i] << std::endl;
        if(max < volumes[i])max = volumes[i];
        if(min > volumes[i])min = volumes[i];
    }
    std::cout << "(max,min) = " << std::endl;
    std::cout << std::setprecision(4) << max/volumes[0] << "," << min/volumes[0] << "," << (max/volumes[0] - min/volumes[0])*100 << "\%"<<std::endl;
    //std::cout << "min:" << min << " ration:" << min/volumes[0] << std::endl;
    //fclose(ofp);
    writing_file.close();
}
