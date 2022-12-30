//
//  gnuplot.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/12/12.
//

#ifndef gnuplot_h
#define gnuplot_h
#include "particle.h"
void outputPLT_P(int nx,int ny,double dx,const char* OutputFileName,std::vector<std::vector<double>>&p){
    FILE *ofp = fopen(OutputFileName,"w");
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){
        fprintf(ofp,"%d %d %lf\n",i,j,p[i][j]);
    }
    fclose(ofp);
}
void outputPLT_M(int nx,int ny,const char* OutputFileName,std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<2>>&map){
    FILE *ofp = fopen(OutputFileName,"w");
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){
        std::vector<int>key = {i,j};
        if(map.find(key) == map.end()){
            fprintf(ofp,"%d %d 0\n",i,j);
            continue;
        }
        else fprintf(ofp,"%d %d %d\n",i,j,(int)map.at(key).size());
    }
    fclose(ofp);
}
void outputPLT_particles(const char* OutputFileName,std::vector<particle> &particles){
    FILE *ofp = fopen(OutputFileName,"w");
    for(auto p:particles){
        //double n = p.velocity.norm();
        double n = p.fixVector.norm();
        if(n > 1.0e-4)fprintf(ofp,"%lf %lf %lf %lf %lf\n ",p.position.x(),p.position.y(),p.fixVector.x(),p.fixVector.y(),n);
    }
    fclose(ofp);
}
#endif /* gnuplot_h */
