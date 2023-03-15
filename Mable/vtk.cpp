//
//  vtk.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "vtk.hpp"
void inputVTK(const char* InputFileName,std::vector<Eigen::Vector3d> &Vertices,std::vector<std::vector<int>> &Faces){
    std::cout << "InputFileName:" << InputFileName << std::endl;
    //std::cout << "OutputFileName:" << OutputFileName << std::endl;
    FILE *ifp = fopen(InputFileName,"r");
    int num_vertices, num_faces,dummy;
    fscanf(ifp, "OFF %d %d %d", &num_vertices, &num_faces, &dummy);
    for(int i=0;i<num_vertices;i++){
        double x,y,z;//点の入力
        fscanf(ifp, "%lf %lf %lf", &x, &y, &z);
        Eigen::Vector3d v(x,y,z);
        Vertices.push_back(v);
    }
    for(int i=0;i<num_faces;i++){//面の入力
        int num_sides, v0, v1, v2;
        fscanf(ifp, "%d %d %d %d", &num_sides, &v0, &v1, &v2);
        Faces.push_back({v0,v1,v2});
    }
    fclose(ifp);
}
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

//    FILE *ofp = fopen(OutputFileName,"w");
//    fprintf(ofp,"# vtk DataFile Version 2.0\n");
//    fprintf(ofp,"Isosurface\n");
//    fprintf(ofp,"ASCII\n");
//    fprintf(ofp,"DATASET STRUCTURED_POINTS\n");
//    fprintf(ofp,"DIMENSIONS %d %d %d\n",nx,ny,nz);
//    fprintf(ofp,"ORIGIN %lf %lf %lf\n",origin[0] ,origin[1] ,origin[2] );
//    fprintf(ofp,"SPACING %lf %lf %lf\n",dx,dx,dx);
//    fprintf(ofp,"POINT_DATA %d\n",Vertices.size);
//    fprintf(ofp,"SCALARS value float 1\n");
//    fprintf(ofp,"LOOKUP_TABLE default\n");
    for(int i=0;i<Vertices.nx;i++){
        for(int j=0;j<Vertices.ny;j++){
            for(int k=0;k<Vertices.nz;k++){
                //fprintf(ofp,"%lf\n",Vertices.value[i][j][k]);
                writing_file << Vertices.value[i][j][k] << std::endl;
            }
        }
    }
    //fclose(ofp);
    writing_file.close();
}
void outputVTK_implicit(std::string OutputFileName,myArray3d &Vertices,double dx,double threshold){
    std::ofstream writing_file;
    writing_file.open(OutputFileName, std::ios::out);
    std::string writing_text ="# vtk DataFile Version 2.0\nIsosurface\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    writing_file << writing_text << std::endl;
    writing_file << "POINTS "<< (int)Vertices.nx*Vertices.ny*Vertices.nz << " float" << std::endl;
    for(int i=0;i<Vertices.nx;i++){
        for(int j=0;j<Vertices.ny;j++){
            for(int k=0;k<Vertices.nz;k++){
                //fprintf(ofp,"%lf %lf %lf\n",(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx);
                writing_file << (i+0.5)*dx <<" "<< (j+0.5)*dx <<" "<< (k+0.5)*dx << std::endl;
            }
        }
    }
    writing_file << "CELL_TYPES "<< (int)Vertices.nx*Vertices.ny*Vertices.nz << " float" << std::endl;
    
    writing_file << "POINT_DATA " << Vertices.size << std::endl;
    for(int i=0;i<(int)Vertices.nx*Vertices.ny*Vertices.nz;i++)writing_file << "1" <<std::endl;
    writing_file << std::endl;
    writing_file << "POINT_DATA " << (int)Vertices.nx*Vertices.ny*Vertices.nz << std::endl;
    writing_file << "SCALARS weight float" << std::endl;
    writing_file << "LOOKUP_TABLE default" << std::endl;
    for(int i=0;i<Vertices.nx;i++){
        for(int j=0;j<Vertices.ny;j++){
            for(int k=0;k<Vertices.nz;k++){
                //fprintf(ofp,"%lf\n",Vertices.value[i][j][k]);
                writing_file << Vertices.value[i][j][k] << std::endl;
            }
        }
    }
    writing_file << "SCALARS radius float" << std::endl;
    writing_file << "LOOKUP_TABLE default" << std::endl;
    for(int i=0;i<Vertices.nx;i++){
        for(int j=0;j<Vertices.ny;j++){
            for(int k=0;k<Vertices.nz;k++){
                if(Vertices.value[i][j][k] > threshold)writing_file << "0.5" <<std::endl;
                else writing_file << "0" <<std::endl;
            }
        }
    }
    writing_file.close();

//    FILE *ofp = fopen(OutputFileName,"w");
//    fprintf(ofp,"# vtk DataFile Version 2.0\n");
//    fprintf(ofp,"Title Data\n");
//    fprintf(ofp,"ASCII\n");
//    fprintf(ofp,"DATASET UNSTRUCTURED_GRID\n");
//    fprintf(ofp,"POINTS %d float\n",(int)Vertices.nx*Vertices.ny*Vertices.nz);
//    for(int i=0;i<Vertices.nx;i++){
//        for(int j=0;j<Vertices.ny;j++){
//            for(int k=0;k<Vertices.nz;k++){
//                fprintf(ofp,"%lf %lf %lf\n",(i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx);
//            }
//        }
//    }
//    fprintf(ofp,"\n");
//    fprintf(ofp,"CELL_TYPES %d\n",(int)Vertices.nx*Vertices.ny*Vertices.nz);
//    for(int i=0;i<(int)Vertices.nx*Vertices.ny*Vertices.nz;i++)fprintf(ofp,"1\n");
//    fprintf(ofp,"\n");
//    fprintf(ofp,"POINT_DATA %d\n",(int)Vertices.nx*Vertices.ny*Vertices.nz);
//    fprintf(ofp,"SCALARS weight float\n");
//    fprintf(ofp,"LOOKUP_TABLE default\n");
//    for(int i=0;i<Vertices.nx;i++){
//        for(int j=0;j<Vertices.ny;j++){
//            for(int k=0;k<Vertices.nz;k++){
//                fprintf(ofp,"%lf\n",Vertices.value[i][j][k]);
//            }
//        }
//    }
//    fprintf(ofp,"SCALARS radius float\n");
//    fprintf(ofp,"LOOKUP_TABLE default\n");
//    for(int i=0;i<Vertices.nx;i++){
//        for(int j=0;j<Vertices.ny;j++){
//            for(int k=0;k<Vertices.nz;k++){
//                if(Vertices.value[i][j][k] > threshold)fprintf(ofp,"0.5\n");
//                else fprintf(ofp,"0\n");
//            }
//        }
//    }
//    fclose(ofp);
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
