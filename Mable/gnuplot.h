//
//  gnuplot.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/12/12.
//

#ifndef gnuplot_h
#define gnuplot_h

void outputPLT_P(int nx,int ny,double dx,const char* OutputFileName,std::vector<std::vector<double>>&p){
    FILE *ofp = fopen(OutputFileName,"w");
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){
        fprintf(ofp,"%lf %lf %lf\n",i*dx,j*dx,p[i][j]);
    }
    fclose(ofp);
}
#endif /* gnuplot_h */
